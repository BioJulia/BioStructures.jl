export
    pdbentrylist,
    pdbstatuslist,
    pdbrecentchanges,
    pdbobsoletelist,
    downloadpdb,
    downloadentirepdb,
    updatelocalpdb,
    downloadallobsoletepdb,
    retrievepdb

const pdb_download_prefix = "https://files.wwpdb.org/pub/pdb"

"""
    pdbentrylist()

Obtain the list of all Protein Data Bank (PDB) entries from the RCSB server.

Requires an internet connection.
"""
function pdbentrylist()
    pdbidlist = String[]
    @info "Fetching the list of all PDB entries from the RCSB server"
    tempfilepath = tempname()
    try
        Downloads.download("$pdb_download_prefix/derived_data/index/entries.idx", tempfilepath)
        open(tempfilepath) do input
            reading = false
            for line in eachline(input)
                # Skip the header lines
                if !reading && !(startswith(line, "IDCODE") || startswith(line, "------") ||
                                 startswith(line, "\t"))
                    reading = true
                end
                if reading
                    # The first 4 characters in the line is the PDB ID
                    pdbid = uppercase(line[1:4])
                    # Check PDB ID is 4 characters long and only consits of alphanumeric characters
                    if !occursin(r"^[a-zA-Z0-9]{4}$", pdbid)
                        throw(ArgumentError("Not a valid PDB ID: \"$pdbid\""))
                    end
                    push!(pdbidlist, pdbid)
                end
            end
        end
    finally
        rm(tempfilepath, force=true)
    end
    return pdbidlist
end

"""
    pdbstatuslist(url::AbstractString)

Obtain the list of Protein Data Bank (PDB) entries from a RCSB weekly status
file by specifying its URL.

An example URL is $pdb_download_prefix/pub/pdb/data/status/latest/added.pdb.
Requires an internet connection.
"""
function pdbstatuslist(url::AbstractString)
    statuslist = String[]
    filename = split(url, "/")[end]
    @info "Fetching weekly status file $filename from the RCSB server"
    tempfilepath = tempname()
    try
        Downloads.download(url, tempfilepath)
        # Some operating systems don't create a file if the download is empty
        if isfile(tempfilepath)
            open(tempfilepath) do input
                for line in eachline(input)
                    # The first 4 characters in the line is the PDB ID
                    pdbid = uppercase(line[1:4])
                    # Check PDB ID is 4 characters long and only consits of alphanumeric characters
                    if !occursin(r"^[a-zA-Z0-9]{4}$", pdbid)
                        throw(ArgumentError("Not a valid PDB ID: \"$pdbid\""))
                    end
                    push!(statuslist, pdbid)
                end
            end
        end
    finally
        rm(tempfilepath, force=true)
    end
    return statuslist
end

"""
    pdbrecentchanges()

Obtain three lists giving the added, modified and obsolete Protein Data Bank
(PDB) entries from the recent RCSB weekly status files.

Requires an internet connection.
"""
function pdbrecentchanges()
    addedlist = pdbstatuslist("$pdb_download_prefix/data/status/latest/added.pdb")
    modifiedlist = pdbstatuslist("$pdb_download_prefix/data/status/latest/modified.pdb")
    obsoletelist = pdbstatuslist("$pdb_download_prefix/data/status/latest/obsolete.pdb")
    return addedlist, modifiedlist, obsoletelist
end

"""
    pdbobsoletelist()

Obtain the list of all obsolete Protein Data Bank (PDB) entries from the RCSB
server.

Requires an internet connection.
"""
function pdbobsoletelist()
    obsoletelist = String[]
    @info "Fetching the list of all obsolete PDB entries from the RCSB server"
    tempfilepath = tempname()
    try
        Downloads.download("$pdb_download_prefix/data/status/obsolete.dat", tempfilepath)
        open(tempfilepath) do input
            for line in eachline(input)
                # Check if its an obsolete pdb entry and not headers
                if line[1:6] == "OBSLTE"
                    # The 21st to 24th characters in obsolete pdb entry has the pdb id
                    pdbid = uppercase(line[21:24])
                    # Check PDB ID is 4 characters long and only consits of alphanumeric characters
                    if !occursin(r"^[a-zA-Z0-9]{4}$", pdbid)
                        throw(ArgumentError("Not a valid PDB ID: \"$pdbid\""))
                    end
                    push!(obsoletelist, pdbid)
                end
            end
        end
    finally
        rm(tempfilepath, force=true)
    end
    return obsoletelist
end

"""
    downloadpdb(pdbid::AbstractString; <keyword arguments>)
    downloadpdb(pdbid::AbstractArray{<:AbstractString, 1}; <keyword arguments>)
    downloadpdb(f::Function, args...)

Download files from the Protein Data Bank (PDB) via RCSB.

When given an `AbstractString`, e.g. `"1AKE"`, downloads the PDB file and
returns the path to the file.
When given an `Array{<:AbstractString, 1}`, downloads the PDB files in the array
and returns an array of the paths to the files.
When given a function as the first argument, runs the function with the
downloaded filepath(s) as an argument then removes the file(s).
Requires an internet connection.

# Arguments
- `dir::AbstractString=pwd()`: the directory to which the PDB file is
    downloaded; defaults to the current working directory.
- `format::Type=PDBFormat`: the format of the PDB file; options are PDBFormat,
    PDBXMLFormat and MMCIFFormat. MMTF files are no longer available to download.
- `obsolete::Bool=false`: if set `true`, the PDB file is downloaded in the
    auto-generated "obsolete" directory inside the specified `dir`.
- `overwrite::Bool=false`: if set `true`, overwrites the PDB file if it exists
    in `dir`; by default skips downloading the PDB file if it exists.
- `ba_number::Integer=0`: if set > 0 downloads the respective biological
    assembly; by default downloads the PDB file.
"""
function downloadpdb(pdbid::AbstractString;
    dir::AbstractString=pwd(),
    format::Type{<:Union{PDBFormat,PDBXMLFormat,MMCIFFormat,BCIFFormat}}=PDBFormat,
    obsolete::Bool=false,
    overwrite::Bool=false,
    ba_number::Integer=0)
    pdbid = uppercase(pdbid)
    # Check PDB ID is 4 characters long and only consits of alphanumeric characters
    if !occursin(r"^[a-zA-Z0-9]{4}$", pdbid)
        throw(ArgumentError("Not a valid PDB ID: \"$pdbid\""))
    end
    # Check if the PDB file is marked as obsolete
    if obsolete
        # Set the download path to obsolete directory inside dir
        dir = joinpath(dir, "obsolete")
    end
    # Check and create directory if it does not exists in filesystem
    if !isdir(dir)
        @info "Creating directory: $dir"
        mkpath(dir)
    end
    # Standard file name format for PDB and biological assembly
    if ba_number == 0
        pdbpath = joinpath(dir, "$pdbid.$(pdbextension[format])")
    else
        pdbpath = joinpath(dir, "$(pdbid)_ba$ba_number.$(pdbextension[format])")
    end
    # Download the PDB file only if it does not exist in the "dir" and when "overwrite" is true
    if isfile(pdbpath) && !overwrite
        @info "File exists: $pdbid"
    else
        # Temporary location to download compressed PDB file.
        archivefilepath = tempname()
        try
            # Download the compressed PDB file to the temporary location
            @info "Downloading file from PDB: $pdbid"
            if ba_number == 0
                if format == BCIFFormat
                    Downloads.download(
                        "https://models.rcsb.org/$pdbid.bcif",
                        pdbpath,
                    )
                    return pdbpath
                else
                    Downloads.download(
                        "http://files.rcsb.org/download/$pdbid.$(pdbextension[format]).gz",
                        archivefilepath,
                    )
                end
            else
                if format == PDBFormat
                    Downloads.download(
                        "http://files.rcsb.org/download/$pdbid.$(pdbextension[format])$ba_number.gz",
                        archivefilepath,
                    )
                elseif format == MMCIFFormat
                    Downloads.download(
                        "http://files.rcsb.org/download/$pdbid-assembly$ba_number.$(pdbextension[format]).gz",
                        archivefilepath,
                    )
                elseif format == BCIFFormat
                    Downloads.download(
                        "https://models.rcsb.org/$pdbid.bcif",
                        archivefilepath,
                    )
                else
                    throw(ArgumentError("Biological assemblies are available in the " *
                                        "PDB and mmCIF formats only"))
                end
            end
            # Verify if the compressed file is downloaded properly and extract it
            if isfile(archivefilepath) && filesize(archivefilepath) > 0
                stream = GzipDecompressorStream(open(archivefilepath))
                open(pdbpath, "w") do output
                    write(output, stream)
                end
                close(stream)
            end
            # Verify if the PDB file is downloaded and extracted without any error
            if !isfile(pdbpath) || filesize(pdbpath) == 0
                if format == PDBFormat
                    throw(ErrorException("Error downloading file: $pdbid; some PDB entries are " *
                                         "not available as PDB format files, consider downloading " *
                                         "the mmCIF file instead"))
                else
                    throw(ErrorException("Error downloading file: $pdbid"))
                end
            end
        finally
            # Remove the temporary compressd PDB file downloaded to clear up space
            rm(archivefilepath, force=true)
        end
    end

    return pdbpath
end

function downloadpdb(pdbidlist::AbstractArray{<:AbstractString,1}; kwargs...)
    pdbpaths = String[]
    failedlist = String[]
    for pdbid in pdbidlist
        try
            pdbpath = downloadpdb(pdbid; kwargs...)
            push!(pdbpaths, pdbpath)
        catch
            @warn "Error downloading PDB: $pdbid"
            push!(failedlist, pdbid)
        end
    end
    if length(failedlist) > 0
        @warn "$(length(failedlist)) PDB file(s) failed to download: $failedlist"
    end
    return pdbpaths
end

function downloadpdb(f::Function, args...; kwargs...)
    pdbpath = downloadpdb(args...; kwargs...)
    try
        f(pdbpath)
    finally
        isa(pdbpath, Array) ? rm.(pdbpath) : rm(pdbpath)
    end
end

"""
    downloadentirepdb(; <keyword arguments>)

Download the entire Protein Data Bank (PDB) from the RCSB server.

Returns the list of PDB IDs downloaded.
Ensure you have enough disk space and time before running.
The function can be stopped any time and called again to resume downloading.
Requires an internet connection.

# Arguments
- `dir::AbstractString=pwd()`: the directory to which the PDB files are
    downloaded; defaults to the current working directory.
- `format::Type=PDBFormat`: the format of the PDB file; options are PDBFormat,
    PDBXMLFormat and MMCIFFormat. MMTF files are no longer available to download.
- `overwrite::Bool=false`: if set `true`, overwrites the PDB file if it exists
    in `dir`; by default skips downloading the PDB file if it exists.
"""
function downloadentirepdb(; dir::AbstractString=pwd(),
    format::Type{<:Union{PDBFormat,PDBXMLFormat,MMCIFFormat}}=PDBFormat,
    overwrite::Bool=false)
    pdblist = pdbentrylist()
    @info "About to download $(length(pdblist)) PDB files, make sure you have enough disk space and time"
    @info "The function can be stopped any time and called again to resume downloading"
    downloadpdb(pdblist, dir=dir, overwrite=overwrite, format=format)
end

"""
    updatelocalpdb(; dir::AbstractString=pwd(), format::Type=PDBFormat)

Update a local copy of the Protein Data Bank (PDB).

Obtains the recent weekly lists of new, modified and obsolete PDB entries and
automatically updates the PDB files of the given `format` inside the local
`dir` directory.
Requires an internet connection.
"""
function updatelocalpdb(; dir::AbstractString=pwd(),
    format::Type{<:Union{PDBFormat,PDBXMLFormat,MMCIFFormat}}=PDBFormat)
    addedlist, modifiedlist, obsoletelist = pdbrecentchanges()
    # Download the newly added and modified pdb files
    downloadpdb(vcat(addedlist, modifiedlist), dir=dir, overwrite=true, format=format)
    # Set the obsolete directory to be inside dir
    obsolete_dir = joinpath(dir, "obsolete")
    for pdbid in obsoletelist
        oldfile = joinpath(dir, "$pdbid.$(pdbextension[format])")
        newfile = joinpath(obsolete_dir, "$pdbid.$(pdbextension[format])")
        # if obsolete pdb is in the "dir", move it to "obsolete" directory inside "dir"
        if isfile(oldfile)
            if !isdir(obsolete_dir)
                mkpath(obsolete_dir)
            end
            mv(oldfile, newfile)
            # If obsolete pdb is already in the obsolete directory, inform the user and skip
        elseif isfile(newfile)
            @info "PDB $pdbid is already moved to the obsolete directory"
            # If obsolete pdb not available in both dir and obsolete, inform the user and skip
        else
            @info "Obsolete PDB $pdbid is missing"
        end
    end
end

"""
    downloadallobsoletepdb(; <keyword arguments>)

Download all obsolete Protein Data Bank (PDB) files from the RCSB server.

Returns the list of PDB IDs downloaded.
Requires an internet connection.

# Arguments
- `obsolete_dir::AbstractString=pwd()`: the directory where the PDB files are
    downloaded; defaults to the current working directory.
- `format::Type=PDBFormat`: the format of the PDB file; options are PDBFormat,
    PDBXMLFormat and MMCIFFormat. MMTF files are no longer available to download.
- `overwrite::Bool=false`: if set `true`, overwrites the PDB file if it exists
    in `dir`; by default skips downloading the PDB file if it exists.
"""
function downloadallobsoletepdb(; obsolete_dir::AbstractString=pwd(),
    format::Type{<:Union{PDBFormat,PDBXMLFormat,MMCIFFormat}}=PDBFormat,
    overwrite::Bool=false)
    obsoletelist = pdbobsoletelist()
    downloadpdb(obsoletelist, dir=obsolete_dir, format=format, overwrite=overwrite)
end

"""
    retrievepdb(pdbid::AbstractString; <keyword arguments>)

Download and read a Protein Data Bank (PDB) file or biological assembly from the
RCSB server, returning a `MolecularStructure`.

Requires an internet connection.

# Arguments
- `pdbid::AbstractString`: the PDB ID to be downloaded and read.
- `dir::AbstractString=pwd()`: the directory to which the PDB file is
    downloaded; defaults to the current working directory.
- `format::Type=MMCIFFormat`: the format of the PDB file; options are PDBFormat,
    PDBXMLFormat and MMCIFFormat. MMTF files are no longer available to download.
- `obsolete::Bool=false`: if set `true`, the PDB file is downloaded in the
    auto-generated "obsolete" directory inside the specified `dir`.
- `overwrite::Bool=false`: if set `true`, overwrites the PDB file if it exists
    in `dir`; by default skips downloading the PDB file if it exists.
- `ba_number::Integer=0`: if set > 0 downloads the respective biological
    assembly; by default downloads the PDB file.
- `structure_name::AbstractString="\$pdbid.pdb"`: the name given to the returned
    `MolecularStructure`; defaults to the PDB ID.
- `remove_disorder::Bool=false`: whether to remove atoms with alt loc ID not ' '
    or 'A'.
- `read_std_atoms::Bool=true`: whether to read standard ATOM records.
- `read_het_atoms::Bool=true`: whether to read HETATOM records.
- `run_dssp::Bool=false`: whether to run DSSP to assign secondary structure.
    Requires the DSSP_jll.jl package to be imported if set to `true`.
- `run_stride::Bool=false`: whether to run STRIDE to assign secondary structure.
    Requires the STRIDE_jll.jl package to be imported if set to `true`.
"""
function retrievepdb(pdbid::AbstractString;
    dir::AbstractString=pwd(),
    format::Type{<:Union{PDBFormat,PDBXMLFormat,MMCIFFormat}}=MMCIFFormat,
    obsolete::Bool=false,
    overwrite::Bool=false,
    ba_number::Integer=0,
    structure_name::AbstractString="$(uppercase(pdbid)).pdb",
    kwargs...)
    downloadpdb(pdbid, dir=dir, format=format, obsolete=obsolete,
        overwrite=overwrite, ba_number=ba_number)
    if obsolete
        # If obsolete is set true, the PDB file is present in the obsolete directory inside dir
        dir = joinpath(dir, "obsolete")
    end
    pdbid_upper = uppercase(pdbid)
    if ba_number == 0
        pdbpath = joinpath(dir, "$pdbid_upper.$(pdbextension[format])")
    else
        pdbpath = joinpath(dir, "$(pdbid_upper)_ba$ba_number.$(pdbextension[format])")
    end
    read(pdbpath, format; structure_name=structure_name, kwargs...)
end
