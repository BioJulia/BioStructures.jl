export run_dssp, run_stride, ss_code!, ss_code, helixselector, sheetselector, coilselector

dssp_executable = `$(DSSP_jll.mkdssp()) --mmcif-dictionary $(DSSP_jll.mmcif_pdbx_dic)`
const dssp_init_line = "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA"

stride_executable = STRIDE_jll.stride_exe()

# Write a PDB file with header and the ATOM records for DSSP/STRIDE to read
function writesubpdb(filepath::AbstractString, el::ProteinStructure)
    open(filepath, "w") do io
        print(
            io,
            "HEADER    XXXXXXXXXXXXXXXXXXXXXX                  25-DEC-23   XXXX              ",
        )
        for at in collectatoms(el)
            checkchainerror(at)
            println(io, pdbline(at))
        end
    end
end

function run_dssp(curr_model::Model, mod_n::Int)
    struc = ProteinStructure()
    struc.models[mod_n] = curr_model
    fixlists!(struc)

    # Write the structure to a temporary PDB file
    tmp_pdb_path = tempname() * ".pdb"
    writesubpdb(tmp_pdb_path, struc)

    # Run DSSP on the temporary PDB file
    dssp_output = try
        readchomp(pipeline(`$dssp_executable --output-format dssp $tmp_pdb_path`))
    catch
        "error"
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
        sscode = "$(line[17])"
        # @info "Assign $(sscode) to $(chainname) $(resnum)"
        ss_code!(curr_model[chainname][resnum], sscode)
    end
    return curr_model
end


"""
    run_dssp(struc::ProteinStructure)

Runs DSSP (Dictionary of Secondary Structure of Proteins) on the given `struc` ProteinStructure object.
This function iterates over all models in the ProteinStructure object and runs DSSP on each model.
"""
function run_dssp(struc::ProteinStructure)
    @info "Structure has $(countmodels(struc)) models, the keys are $(keys(struc.models))"
    for (mod_n, curr_model) in models(struc)
        struc[mod_n] = run_dssp(curr_model, mod_n)
    end
    return struc
end

function run_stride(curr_model::Model, mod_n::Int)
    struc = ProteinStructure()
    struc.models[mod_n] = curr_model
    fixlists!(struc)

    # Write the structure to a temporary PDB file
    tmp_pdb_path = tempname() * ".pdb"
    writesubpdb(tmp_pdb_path, struc)

    # Run STRIDE on the temporary PDB file
    stride_output = try
        readchomp(pipeline(`$stride_executable $tmp_pdb_path`))
    catch
        "error"
    end

    # Delete the temporary PDB file
    rm(tmp_pdb_path)

    end_of_line = Sys.iswindows() ? "\r\n" : "\n"
    for line in split(stride_output, end_of_line)
        if startswith(line, "ASG")
            chainname = line[10:10]
            resnum = "$(parse(Int, line[12:15]))"
            sscode = uppercase(line[25:25]) # stride sometimes returns "b" instead of "B"
            # @info "Assign $(sscode) to $(chainname) $(resnum)"
            ss_code!(curr_model[chainname][resnum], sscode)
        end
    end
    return curr_model
end

"""
    run_stride(struc::ProteinStructure)

Runs the STRIDE algorithm to assign secondary structure to a protein structure.
This function iterates over all models in the ProteinStructure object and runs STRIDE on each model.
"""
function run_stride(struc::ProteinStructure)
    @info "Structure has $(countmodels(struc)) models, the keys are $(keys(struc.models))"
    for (mod_n, curr_model) in models(struc)
        struc[mod_n] = run_stride(curr_model, mod_n)
    end
    return struc
end

"""
    ss_code!(res::Residue, ss_code::String)

Set the secondary structure code of a residue.
"""
ss_code!(res::Residue, ss_code::String) = (res.ss_code = ss_code)

"""
    ss_code(res::Residue)
    ss_code(at::Atom)

Get the secondary structure code of a residue/atom.
"""
ss_code(res::Residue) = res.ss_code
ss_code(at::AbstractAtom) = ss_code(residue(at))

"""
    helixselector(el)

Returns `true` if the secondary structure code of the given `el` is "H", "G", or "I", indicating that it is part of a helix.
"""
function helixselector(el::Union{Residue,Atom})
    sscode = ss_code(el)
    return sscode == "H" || sscode == "G" || sscode == "I"
end

"""
    sheetselector(el)

Returns `true` if the secondary structure code of the given `el` is "E" (extended strand) or "B" (residue in isolated beta-bridge).
Otherwise, returns `false`.
"""
function sheetselector(el::Union{Residue,Atom})
    sscode = ss_code(el)
    return sscode == "E" || sscode == "B"
end

"""
    coilselector(el)

Return `true` if the secondary structure code of the given `el` is "T", "S", or " ", indicating that it is a coil or loop region.
"""
function coilselector(el::Union{Residue,Atom})
    sscode = ss_code(el)
    return sscode == "T" || sscode == "S" || sscode == " "
end
