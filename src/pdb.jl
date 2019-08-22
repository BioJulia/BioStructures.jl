export
    PDB,
    PDBXML,
    MMCIF,
    MMTF,
    pdbextension,
    PDBParseError,
    pdbentrylist,
    pdbstatuslist,
    pdbrecentchanges,
    pdbobsoletelist,
    downloadpdb,
    downloadentirepdb,
    updatelocalpdb,
    downloadallobsoletepdb,
    retrievepdb,
    readpdb,
    spaceatomname,
    pdbline,
    writepdb


# PDB file formats

"Protein Data Bank (PDB) file format."
struct PDB <: BioCore.IO.FileFormat end

"Protein Data Bank (PDB) XML file format."
struct PDBXML <: BioCore.IO.FileFormat end

"Protein Data Bank (PDB) mmCIF file format."
struct MMCIF <: BioCore.IO.FileFormat end

"Protein Data Bank (PDB) MMTF file format."
struct MMTF <: BioCore.IO.FileFormat end

"Mapping of Protein Data Bank (PDB) formats to their file extensions."
const pdbextension = Dict{Type, String}(PDB=> ".pdb", PDBXML=> ".xml",
                                        MMCIF=> ".cif", MMTF=> ".mmtf")


"Error arising from parsing a Protein Data Bank (PDB) file."
struct PDBParseError <: Exception
    message::String
    line_number::Int
    line::String
end

function Base.showerror(io::IO, e::PDBParseError)
    return print(io,
            e.message,
            " at line ",
            e.line_number,
            " of file:\n",
            e.line)
end


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
        download("ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx", tempfilepath)
        open(tempfilepath) do input
            # Skips the first two lines as it contains headers
            linecount = 1
            for line in eachline(input)
                if linecount > 2
                    # The first 4 characters in the line is the PDB ID
                    pdbid = uppercase(line[1:4])
                    # Check PDB ID is 4 characters long and only consits of alphanumeric characters
                    if !occursin(r"^[a-zA-Z0-9]{4}$", pdbid)
                        throw(ArgumentError("Not a valid PDB ID: \"$pdbid\""))
                    end
                    push!(pdbidlist, pdbid)
                end
                linecount +=1
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

An example URL is ftp://ftp.wwpdb.org/pub/pdb/data/status/latest/added.pdb.
Requires an internet connection.
"""
function pdbstatuslist(url::AbstractString)
    statuslist = String[]
    filename = split(url, "/")[end]
    @info "Fetching weekly status file $filename from the RCSB server"
    tempfilepath = tempname()
    try
        download(url, tempfilepath)
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
    addedlist = pdbstatuslist("ftp://ftp.wwpdb.org/pub/pdb/data/status/latest/added.pdb")
    modifiedlist = pdbstatuslist("ftp://ftp.wwpdb.org/pub/pdb/data/status/latest/modified.pdb")
    obsoletelist = pdbstatuslist("ftp://ftp.wwpdb.org/pub/pdb/data/status/latest/obsolete.pdb")
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
        download("ftp://ftp.wwpdb.org/pub/pdb/data/status/obsolete.dat", tempfilepath)
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
    downloadpdb(pdbid::Array{String, 1}; <keyword arguments>)
    downloadpdb(f::Function, args...)

Download files from the Protein Data Bank (PDB) via RCSB.

When given an `AbstractString`, e.g. `"1AKE"`, downloads the PDB file and
returns the path to the file.
When given an `Array{String, 1}`, downloads the PDB files in the array and
returns an array of the paths to the files.
When given a function as the first argument, runs the function with the
downloaded filepath(s) as an argument then removes the file(s).
Requires an internet connection.

# Arguments
- `pdb_dir::AbstractString=pwd()`: the directory to which the PDB file is
    downloaded; defaults to the current working directory.
- `file_format::Type=PDB`: the format of the PDB file; options are PDB, PDBXML,
    MMCIF and MMTF.
- `obsolete::Bool=false`: if set `true`, the PDB file is downloaded in the
    auto-generated "obsolete" directory inside the specified `pdb_dir`.
- `overwrite::Bool=false`: if set `true`, overwrites the PDB file if it exists
    in `pdb_dir`; by default skips downloading the PDB file if it exists.
- `ba_number::Integer=0`: if set > 0 downloads the respective biological
    assembly; by default downloads the PDB file.
"""
function downloadpdb(pdbid::AbstractString;
                pdb_dir::AbstractString=pwd(),
                file_format::Type=PDB,
                obsolete::Bool=false,
                overwrite::Bool=false,
                ba_number::Integer=0)
    pdbid = uppercase(pdbid)
    # Check PDB ID is 4 characters long and only consits of alphanumeric characters
    if !occursin(r"^[a-zA-Z0-9]{4}$", pdbid)
        throw(ArgumentError("Not a valid PDB ID: \"$pdbid\""))
    end
    # Check if PDB file format is valid
    if !haskey(pdbextension, file_format)
        throw(ArgumentError("Invalid PDB file format!"))
    end
    # Check if the PDB file is marked as obsolete
    if obsolete
        # Set the download path to obsolete directory inside the "pdb_dir"
        pdb_dir = joinpath(pdb_dir, "obsolete")
    end
    # Check and create directory if it does not exists in filesystem
    if !isdir(pdb_dir)
        @info "Creating directory: $pdb_dir"
        mkpath(pdb_dir)
    end
    # Standard file name format for PDB and biological assembly
    if ba_number==0
        pdbpath = joinpath(pdb_dir, "$pdbid$(pdbextension[file_format])")
    else
        pdbpath = joinpath(pdb_dir, "$(pdbid)_ba$ba_number$(pdbextension[file_format])")
    end
    # Download the PDB file only if it does not exist in the "pdb_dir" and when "overwrite" is true
    if isfile(pdbpath) && !overwrite
        @info "PDB exists: $pdbid"
    else
        # Temporary location to download compressed PDB file.
        archivefilepath = tempname()
        try
            # Download the compressed PDB file to the temporary location
            @info "Downloading PDB: $pdbid"
            if ba_number == 0
                if file_format == PDB || file_format == PDBXML || file_format == MMCIF
                    download("http://files.rcsb.org/download/$pdbid$(pdbextension[file_format]).gz", archivefilepath)
                else
                    download("http://mmtf.rcsb.org/v1.0/full/$pdbid.mmtf.gz", archivefilepath)
                end
            else
                if file_format == PDB
                    download("http://files.rcsb.org/download/$pdbid$(pdbextension[file_format])$ba_number.gz", archivefilepath)
                elseif file_format == MMCIF
                    download("http://files.rcsb.org/download/$pdbid-assembly$ba_number$(pdbextension[file_format]).gz", archivefilepath)
                else
                    throw(ArgumentError("Biological assemblies are available in the PDB and mmCIF formats only"))
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
                throw(ErrorException("Error downloading PDB: $pdbid"))
            end
        finally
            # Remove the temporary compressd PDB file downloaded to clear up space
            rm(archivefilepath, force=true)
        end
    end

    return pdbpath
end

function downloadpdb(pdbidlist::Array{String, 1}; kwargs...)
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
- `pdb_dir::AbstractString=pwd()`: the directory to which the PDB files are
    downloaded; defaults to the current working directory.
- `file_format::Type=PDB`: the format of the PDB file; options are PDB, PDBXML,
    MMCIF and MMTF.
- `overwrite::Bool=false`: if set `true`, overwrites the PDB file if it exists
    in `pdb_dir`; by default skips downloading the PDB file if it exists.
"""
function downloadentirepdb(; pdb_dir::AbstractString=pwd(), file_format::Type=PDB, overwrite::Bool=false)
    pdblist = pdbentrylist()
    @info "About to download $(length(pdblist)) PDB files, make sure you have enough disk space and time"
    @info "The function can be stopped any time and called again to resume downloading"
    downloadpdb(pdblist, pdb_dir=pdb_dir, overwrite=overwrite, file_format=file_format)
end

"""
    updatelocalpdb(; pdb_dir::AbstractString=pwd(), file_format::Type=PDB)

Update a local copy of the Protein Data Bank (PDB).

Obtains the recent weekly lists of new, modified and obsolete PDB entries and
automatically updates the PDB files of the given `file_format` inside the local
`pdb_dir` directory.
Requires an internet connection.
"""
function updatelocalpdb(; pdb_dir::AbstractString=pwd(), file_format::Type=PDB)
    addedlist, modifiedlist, obsoletelist = pdbrecentchanges()
    # Download the newly added and modified pdb files
    downloadpdb(vcat(addedlist, modifiedlist), pdb_dir=pdb_dir, overwrite=true, file_format=file_format)
    # Set the obsolete directory to be inside pdb_dir
    obsolete_dir=joinpath(pdb_dir, "obsolete")
    for pdbid in obsoletelist
        oldfile = joinpath(pdb_dir, "$pdbid$(pdbextension[file_format])")
        newfile = joinpath(obsolete_dir, "$pdbid$(pdbextension[file_format])")
        # if obsolete pdb is in the "pdb_dir", move it to "obsolete" directory inside "pdb_dir"
        if isfile(oldfile)
            if !isdir(obsolete_dir)
                mkpath(obsolete_dir)
            end
            mv(oldfile, newfile)
        # If obsolete pdb is already in the obsolete directory, inform the user and skip
        elseif isfile(newfile)
            @info "PDB $pdbid is already moved to the obsolete directory"
        # If obsolete pdb not available in both pdb_dir and obsolete, inform the user and skip
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
- `file_format::Type=PDB`: the format of the PDB file; options are PDB, PDBXML,
    MMCIF and MMTF.
- `overwrite::Bool=false`: if set `true`, overwrites the PDB file if it exists
    in `pdb_dir`; by default skips downloading the PDB file if it exists.
"""
function downloadallobsoletepdb(; obsolete_dir::AbstractString=pwd(), file_format::Type=PDB, overwrite::Bool=false)
    obsoletelist = pdbobsoletelist()
    downloadpdb(obsoletelist, pdb_dir=obsolete_dir, file_format=file_format, overwrite=overwrite)
end


"""
    retrievepdb(pdbid::AbstractString; <keyword arguments>)

Download and read a Protein Data Bank (PDB) file or biological assembly from the
RCSB server, returning a `ProteinStructure`.

Requires an internet connection.

# Arguments
- `pdbid::AbstractString`: the PDB ID to be downloaded and read.
- `pdb_dir::AbstractString=pwd()`: the directory to which the PDB file is
    downloaded; defaults to the current working directory.
- `obsolete::Bool=false`: if set `true`, the PDB file is downloaded in the
    auto-generated "obsolete" directory inside the specified `pdb_dir`.
- `overwrite::Bool=false`: if set `true`, overwrites the PDB file if it exists
    in `pdb_dir`; by default skips downloading the PDB file if it exists.
- `ba_number::Integer=0`: if set > 0 downloads the respective biological
    assembly; by default downloads the PDB file.
- `structure_name::AbstractString="\$pdbid.pdb"`: the name given to the returned
    `ProteinStructure`; defaults to the PDB ID.
- `remove_disorder::Bool=false`: whether to remove atoms with alt loc ID not ' '
    or 'A'.
- `read_std_atoms::Bool=true`: whether to read standard ATOM records.
- `read_het_atoms::Bool=true`: whether to read HETATOM records.
"""
function retrievepdb(pdbid::AbstractString;
            pdb_dir::AbstractString=pwd(),
            obsolete::Bool=false,
            overwrite::Bool=false,
            ba_number::Integer=0,
            structure_name::AbstractString="$(uppercase(pdbid)).pdb",
            kwargs...)
    downloadpdb(pdbid, pdb_dir=pdb_dir, obsolete=obsolete, overwrite=overwrite, ba_number=ba_number)
    if obsolete
        # If obsolete is set true, the PDB file is present in the obsolete directory inside "pdb_dir"
        pdb_dir = joinpath(pdb_dir, "obsolete")
    end
    readpdb(pdbid; pdb_dir=pdb_dir, ba_number=ba_number, structure_name=structure_name, kwargs...)
end

"""
    readpdb(pdbid::AbstractString; <keyword arguments>)

Read a Protein Data Bank (PDB) file and return a `ProteinStructure`.

This is an alternative to `read("\$pdb_dir/\$pdbid.pdb", PDB)` but is
effectively the same.

# Arguments
- `pdbid::AbstractString`: the PDB ID to be read.
- `pdb_dir::AbstractString=pwd()`: the directory to which the PDB file is
    downloaded; defaults to the current working directory.
- `ba_number::Integer=0`: if set > 0 downloads the respective biological
    assembly; by default downloads the PDB file.
- `structure_name::AbstractString="\$pdbid.pdb"`: the name given to the returned
    `ProteinStructure`; defaults to the PDB ID.
- `remove_disorder::Bool=false`: whether to remove atoms with alt loc ID not ' '
    or 'A'.
- `read_std_atoms::Bool=true`: whether to read standard ATOM records.
- `read_het_atoms::Bool=true`: whether to read HETATOM records.
"""
function readpdb(pdbid::AbstractString;
            pdb_dir::AbstractString=pwd(),
            ba_number::Integer=0,
            structure_name::AbstractString="$pdbid.pdb",
            kwargs...)
    pdbid = uppercase(pdbid)
    # Standard file name format for PDB and biological assembly
    if ba_number==0
        pdbpath = joinpath(pdb_dir, "$pdbid.pdb")
    else
        pdbpath = joinpath(pdb_dir, "$(pdbid)_ba$ba_number.pdb")
    end
    read(pdbpath, PDB; structure_name=structure_name, kwargs...)
end


function Base.read(input::IO,
            ::Type{PDB};
            structure_name::AbstractString="",
            remove_disorder::Bool=false,
            read_std_atoms::Bool=true,
            read_het_atoms::Bool=true)
    # Define ProteinStructure and add to it incrementally
    struc = ProteinStructure(structure_name)
    struc[1] = Model(1, struc)
    # Entries outside of a MODEL/ENDMDL block are added to model 1
    curr_model = 1
    line_n = 0
    for line in eachline(input)
        line_n += 1
        # Read ATOM and HETATM records as required
        if (read_std_atoms && startswith(line, "ATOM  ")) ||
                (read_het_atoms && startswith(line, "HETATM"))
            unsafe_addatomtomodel!(
                    struc[curr_model],
                    AtomRecord(line, line_n), remove_disorder=remove_disorder)
        # Read MODEL record
        elseif startswith(line, "MODEL ")
            try
                curr_model = parse(Int, line[11:min(14, end)])
            catch
                throw(PDBParseError(
                    "Could not read model serial number", line_n, line))
            end
            # Create model if required
            if !haskey(models(struc), curr_model)
                struc[curr_model] = Model(curr_model, struc)
            end
        # Read ENDMDL record
        elseif startswith(line, "ENDMDL")
            curr_model = 1
        end
    end
    # Remove any models that were not added to
    for model in struc
        if countchains(model) == 0
            delete!(models(struc), modelnumber(model))
        end
    end
    # Generate lists for iteration
    fixlists!(struc)
    return struc
end

function Base.read(filepath::AbstractString,
            t::Type{<:Union{PDB, MMCIF}};
            structure_name::AbstractString=splitdir(filepath)[2],
            kwargs...)
    open(filepath) do input
        read(input, t; structure_name=structure_name, kwargs...)
    end
end


# Constructor from PDB ATOM/HETATM line
AtomRecord(pdb_line::String, line_n::Integer=1) = AtomRecord(
    pdb_line[1] == 'H', # This assumes the line has already been checked as an ATOM/HETATM record
    parseserial(pdb_line, line_n),
    parseatomname(pdb_line, line_n),
    parsealtloc(pdb_line, line_n),
    parseresname(pdb_line, line_n),
    parsechainid(pdb_line, line_n),
    parseresnumber(pdb_line, line_n),
    parseinscode(pdb_line, line_n),
    [
        parsecoordx(pdb_line, line_n),
        parsecoordy(pdb_line, line_n),
        parsecoordz(pdb_line, line_n)
    ],
    parseoccupancy(pdb_line),
    parsetempfac(pdb_line),
    parseelement(pdb_line),
    parsecharge(pdb_line)
)


function parseserial(line::String, line_n::Integer=1)
    try
        return parse(Int, line[7:11])
    catch
        throw(PDBParseError("Could not read atom serial number", line_n, line))
    end
end

function parseatomname(line::String, line_n::Integer=1)
    try
        return line[13:16]
    catch
        throw(PDBParseError("Could not read atom name", line_n, line))
    end
end

function parsealtloc(line::String, line_n::Integer=1)
    try
        return line[17]
    catch
        throw(PDBParseError("Could not read alt loc identifier", line_n, line))
    end
end

function parseresname(line::String, line_n::Integer=1)
    try
        return line[18:20]
    catch
        throw(PDBParseError("Could not read residue name", line_n, line))
    end
end

function parsechainid(line::String, line_n::Integer=1)
    try
        return string(line[22])
    catch
        throw(PDBParseError("Could not read chain ID", line_n, line))
    end
end

function parseresnumber(line::String, line_n::Integer=1)
    try
        return parse(Int, line[23:26])
    catch
        throw(PDBParseError("Could not read residue number", line_n, line))
    end
end

function parseinscode(line::String, line_n::Integer=1)
    try
        return line[27]
    catch
        throw(PDBParseError("Could not read insertion code", line_n, line))
    end
end

function parsecoordx(line::String, line_n::Integer=1)
    try
        return parse(Float64, line[31:38])
    catch
        throw(PDBParseError("Could not read x coordinate", line_n, line))
    end
end

function parsecoordy(line::String, line_n::Integer=1)
    try
        return parse(Float64, line[39:46])
    catch
        throw(PDBParseError("Could not read y coordinate", line_n, line))
    end
end

function parsecoordz(line::String, line_n::Integer=1)
    try
        return parse(Float64, line[47:54])
    catch
        throw(PDBParseError("Could not read z coordinate", line_n, line))
    end
end

function parseoccupancy(line::String)
    try
        return parse(Float64, line[55:60])
    catch
        return 1.0
    end
end

function parsetempfac(line::String)
    try
        return parse(Float64, line[61:66])
    catch
        return 0.0
    end
end

function parseelement(line::String)
    try
        return line[77:78]
    catch
        return "  "
    end
end

function parsecharge(line::String)
    try
        return line[79:80]
    catch
        return "  "
    end
end


# Form a string of a certain length from a value by adding spaces to the left
# Throws an error if the value is too long
function spacestring(val_in, new_length::Integer)
    string_out = string(val_in)
    if length(string_out) > new_length
        throw(ArgumentError("Cannot fit value \"$string_out\" into $new_length space(s)"))
    end
    return lpad(string_out, new_length)
end

"""
    spaceatomname(at::Atom)

Space an `Atom` name such that the last element letter (generally) appears in
the second column.

If the `element` property of the `Atom` is set it is used to get the element,
otherwise the name starts from the second column where possible.
This function is generally not required as spacing is recorded when atom names
are read in from a Protein Data Bank (PDB) file.
However this spacing can be important, for example distinguising between C-alpha
and calcium atoms.
"""
function spaceatomname(at::Atom)
    at_name = atomname(at, strip=false)
    chars = length(at_name)
    if chars == 4
        return at_name
    end
    strip_el = element(at, strip=true)
    if chars > 4
        throw(ArgumentError("Atom name is greater than four characters: \"$at_name\""))
    end
    # In the absence of the element, the first index goes in column two
    if strip_el == "" || findfirst(isequal(strip_el[1]), at_name) == nothing
        cent_ind = 1
    # The last letter of the element goes in column two where possible
    else
        cent_ind = something(findfirst(isequal(strip_el[1]), at_name), 0) + length(strip_el) - 1
    end
    if cent_ind > 2
        throw(ArgumentError("Atom name is too long to space correctly: \"$at_name\""))
    end
    if cent_ind == 1 && chars < 4
        out_string = " $at_name"
    else
        out_string = "$at_name"
    end
    return rpad(out_string, 4)
end

# Decimal places to format output to
const coordspec = FormatSpec(".3f")
const floatspec = FormatSpec(".2f")

"""
    pdbline(at::Atom)
    pdbline(at::AtomRecord)

Form a Protein Data Bank (PDB) format ATOM or HETATM record as a `String` from
an `Atom` or `AtomRecord`.

This will throw an `ArgumentError` if a value cannot fit into the allocated
space, e.g. the chain ID is longer than one character or the atom serial is
greater than 99999.
In this case consider using `writemmcif` to write a mmCIF file.
"""
function pdbline(at::Atom)
    return (ishetero(at) ? "HETATM" : "ATOM  ") *
            # This will throw an error for serial numbers over 99999
            spacestring(serial(at), 5) *
            " " *
            spaceatomname(at) *
            string(altlocid(at)) *
            spacestring(resname(at), 3) *
            " " *
            chainid(at) * # This is checked during writing for being one character
            spacestring(resnumber(at), 4) *
            string(inscode(at)) *
            "   " *
            # This will throw an error for large coordinate values, e.g. -1000.123
            spacestring(pyfmt(coordspec, round(x(at), digits=3)), 8) *
            spacestring(pyfmt(coordspec, round(y(at), digits=3)), 8) *
            spacestring(pyfmt(coordspec, round(z(at), digits=3)), 8) *
            spacestring(pyfmt(floatspec, round(occupancy(at), digits=2)), 6) *
            # This will throw an error for large temp facs, e.g. 1000.12
            spacestring(pyfmt(floatspec, round(tempfactor(at), digits=2)), 6) *
            "          " *
            spacestring(element(at), 2) *
            spacestring(charge(at), 2)
end

function pdbline(at_rec::AtomRecord)
    return (at_rec.het_atom ? "HETATM" : "ATOM  ") *
            # This will throw an error for serial numbers over 99999
            spacestring(at_rec.serial, 5) *
            " " *
            spacestring(at_rec.atom_name, 4) *
            string(at_rec.alt_loc_id) *
            spacestring(at_rec.res_name, 3) *
            " " *
            at_rec.chain_id * # This is checked during writing for being one character
            spacestring(at_rec.res_number, 4) *
            string(at_rec.ins_code) *
            "   " *
            # This will throw an error for large coordinate values, e.g. -1000.123
            spacestring(pyfmt(coordspec, round(at_rec.coords[1], digits=3)), 8) *
            spacestring(pyfmt(coordspec, round(at_rec.coords[2], digits=3)), 8) *
            spacestring(pyfmt(coordspec, round(at_rec.coords[3], digits=3)), 8) *
            spacestring(pyfmt(floatspec, round(at_rec.occupancy, digits=2)), 6) *
            # This will throw an error for large temp facs, e.g. 1000.12
            spacestring(pyfmt(floatspec, round(at_rec.temp_factor, digits=2)), 6) *
            "          " *
            spacestring(at_rec.element, 2) *
            spacestring(at_rec.charge, 2)
end

function checkchainerror(el::Union{Chain, AbstractResidue, AbstractAtom})
    if length(chainid(el)) != 1
        throw(ArgumentError("Cannot write valid PDB file with non-single character chain ID: \"$(chainid(el))\""))
    end
end


"""
    writepdb(output, element, atom_selectors...)

Write a `StructuralElementOrList` to a Protein Data Bank (PDB) format file or
output stream.

Only ATOM, HETATM, MODEL and ENDMDL records are written - there is no header and
there are no TER records.
Atom selector functions can be given as additional arguments - only atoms that
return `true` from all the functions are retained.
"""
function writepdb(output::IO,
                el::Union{ProteinStructure, Vector{Model}},
                atom_selectors::Function...)
    # If there are multiple models, write out MODEL/ENDMDL lines
    if length(el) > 1
        for mod in sort(collect(el))
            println(output, "MODEL     ", spacestring(modelnumber(mod), 4), repeat(" ", 66))
            writepdb(output, mod, atom_selectors...)
            println(output, "ENDMDL$(repeat(" ", 74))")
        end
    # If there is only one model, do not write out MODEL/ENDMDL lines
    elseif isa(el, ProteinStructure)
        writepdb(output, defaultmodel(el), atom_selectors...)
    else
        writepdb(output, el[1], atom_selectors...)
    end
end

function writepdb(output::IO,
                el::Union{Model, Chain, AbstractResidue, Vector{Chain},
                    Vector{AbstractResidue}, Vector{Residue}, Vector{DisorderedResidue}},
                atom_selectors::Function...)
    # Collect residues then expand out disordered residues
    for res in collectresidues(el)
        checkchainerror(res)
        if isa(res, Residue)
            for at in collectatoms(res, atom_selectors...)
                writepdb(output, at)
            end
        else
            for res_name in resnames(res)
                for at in collectatoms(disorderedres(res, res_name), atom_selectors...)
                    writepdb(output, at)
                end
            end
        end
    end
end

function writepdb(output::IO, at::AbstractAtom, atom_selectors::Function...)
    checkchainerror(at)
    for atom_record in at
        println(output, pdbline(atom_record))
    end
end

function writepdb(output::IO,
                ats::Vector{<:AbstractAtom},
                atom_selectors::Function...)
    for at in collectatoms(ats, atom_selectors...)
        writepdb(output, at)
    end
end

function writepdb(filepath::AbstractString,
                el::StructuralElementOrList,
                atom_selectors::Function...)
    open(filepath, "w") do output
        writepdb(output, el, atom_selectors...)
    end
end
