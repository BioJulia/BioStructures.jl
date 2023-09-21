export rundssp!, rundssp, runstride!, runstride, sscode!, sscode, helixselector, sheetselector, coilselector, sscodeselector

const dssp_executable = `$(DSSP_jll.mkdssp()) --mmcif-dictionary $(DSSP_jll.mmcif_pdbx_dic) --output-format dssp`
const stride_executable = `$(STRIDE_jll.stride_exe())`
const dssp_init_line = "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA"
const header_line = "HEADER"
const helixsscodes = Set(["G", "H", "I", "P"])
const sheetsscodes = Set(["E", "B"])
const coilsscodes = Set(["C", "T", "S", " "])

"""
    rundssp(curr_model::Model, mod_n::Int)
"""
function rundssp(curr_model::Model, mod_n::Int)
    struc = ProteinStructure()
    struc.models[mod_n] = curr_model
    fixlists!(struc)

    # Write the structure to a temporary PDB file
    tmp_pdb_path = tempname() * ".pdb"
    open(tmp_pdb_path, "w") do io
        print(io, header_line)
        writepdb(io, struc)
    end

    # Run DSSP on the temporary PDB file
    dssp_output = try
        readchomp(pipeline(`$dssp_executable $tmp_pdb_path`))
    catch
        "Error: DSSP failed to run."
    end

    # Delete the temporary PDB file
    rm(tmp_pdb_path)

    data_begin = false
    end_of_line = Sys.iswindows() ? "\r\n" : "\n"
    for line in split(dssp_output, end_of_line)
        if line == dssp_init_line
            data_begin = true
            continue
        end
        !data_begin && continue
        line[14] == '!' && continue
        chainname = "$(line[12])"
        resnum = (parse(Int, line[6:10]))
        resnum = "$(resnum)"
        ss_code = "$(line[17])"
        sscode!(curr_model[chainname][resnum], ss_code)
    end
    return curr_model
end


"""
    rundssp(struc::ProteinStructure)

Runs DSSP (Dictionary of Secondary Structure of Proteins) on the given `struc` ProteinStructure object.
This function iterates over all models in the ProteinStructure object and runs DSSP on each model.
"""
function rundssp(struc::ProteinStructure)
    if countmodels(struc) > 1
        @info "Structure has $(countmodels(struc)) models, the keys are $(keys(struc.models))."
    end
    for (mod_n, curr_model) in models(struc)
        struc[mod_n] = rundssp(curr_model, mod_n)
    end
    return struc
end

"""
    rundssp!(struc::ProteinStructure)

Using DSSP to set the ss_code field of each residue in the given `struc` ProteinStructure object.
"""
rundssp!(struc::ProteinStructure) = (struc = rundssp(struc))

function runstride(curr_model::Model, mod_n::Int)
    struc = ProteinStructure()
    struc.models[mod_n] = curr_model
    fixlists!(struc)

    # Write the structure to a temporary PDB file
    tmp_pdb_path = tempname() * ".pdb"
    open(tmp_pdb_path, "w") do io
        print(io, header_line)
        writepdb(io, struc)
    end

    # Run STRIDE on the temporary PDB file
    stride_output = try
        readchomp(pipeline(`$stride_executable $tmp_pdb_path`))
    catch
        "Error: STRIDE failed to run."
    end

    # Delete the temporary PDB file
    rm(tmp_pdb_path)

    end_of_line = Sys.iswindows() ? "\r\n" : "\n"
    for line in split(stride_output, end_of_line)
        if startswith(line, "ASG")
            chainname = line[10:10]
            resnum = "$(parse(Int, line[12:15]))"
            ss_code = uppercase(line[25:25]) # stride sometimes returns "b" instead of "B"
            # @info "Assign $(sscode) to $(chainname) $(resnum)"
            sscode!(curr_model[chainname][resnum], ss_code)
        end
    end
    return curr_model
end

"""
    runstride(struc::ProteinStructure)

Runs the STRIDE algorithm to assign secondary structure to a protein structure.
This function iterates over all models in the ProteinStructure object and runs STRIDE on each model.
"""
function runstride(struc::ProteinStructure)
    if countmodels(struc) > 1
        @info "Structure has $(countmodels(struc)) models, the keys are $(keys(struc.models))"
    end

    for (mod_n, curr_model) in models(struc)
        struc[mod_n] = runstride(curr_model, mod_n)
    end
    return struc
end

"""
    runstride!(struc::ProteinStructure)

Using STRIDE to set the ss_code field of each residue in the given `struc` ProteinStructure object.
"""
runstride!(struc::ProteinStructure) = (struc = runstride(struc))


"""
    sscode!(res::Residue, ss_code::String)
    sscode!(at::Atom, ss_code::String)

Set the secondary structure code of a atom/residue.
"""
sscode!(res::Residue, ss_code::String) = (res.ss_code = ss_code)
sscode!(at::AbstractAtom, ss_code::String) = sscode!(residue(at), ss_code)

"""
    sscode(res::Residue)
    sscode(at::AbstractAtom)

Get the secondary structure code of a atom/residue.
"""
sscode(res::Residue) = res.ss_code
sscode(at::AbstractAtom) = sscode(residue(at))


"""
    helixselector(at)
    helixselector(res)

Returns `true` if the secondary structure code of the given `Atom` or `Residue` is "G", "H", "I" or "P", indicating that it is part of a helix.
"""
function helixselector(el::Union{Residue,Atom})
    return sscode(el) in helixsscodes
end

"""
    sheetselector(at)
    sheetselector(res)

Returns `true` if the secondary structure code of the given `Atom` or `Residue` is "E" (extended strand) or "B" (residue in isolated beta-bridge).
Otherwise, returns `false`.
"""
function sheetselector(el::Union{Residue,Atom})
    return sscode(el) in sheetsscodes
end

"""
    coilselector(at)
    coilselector(res)

Return `true` if the secondary structure code of the given `Atom` or `Residue` is "C" "T", "S", or " ", indicating that it is a coil or loop region.
"""
function coilselector(el::Union{Residue,Atom})
    return sscode(el) in coilsscodes
end

"""
    sscodeselector(at, ss_codes)
    sscodeselector(res, ss_codes)

Return `true` if the secondary structure code of `Atom` or `Residue` is in `ss_codes`, `false` otherwise.
"""
function sscodeselector(el::Union{AbstractResidue, AbstractAtom},
    ss_codes::Union{Set{String}, Vector{String}})
    return sscode(el) in ss_codes
end
