export
    rundssp!,
    rundssp,
    runstride!,
    runstride,
    sscode,
    sscode!,
    sscodeselector,
    helixsscodes,
    sheetsscodes,
    coilsscodes,
    helixselector,
    sheetselector,
    coilselector

const dssp_executable = `$(DSSP_jll.mkdssp()) --mmcif-dictionary $(DSSP_jll.mmcif_pdbx_dic) --output-format dssp`
const stride_executable = `$(STRIDE_jll.stride_exe())`
const ss_code_unassigned = '-'

"""
    rundssp!(struc)
    rundssp!(mod)

Run DSSP (Define Secondary Structure of Proteins) on the given structural element
to assign secondary structure.

A temporary PDB file is written, so this will fail if the structural element cannot
be written to a PDB file, for example if there are two-letter chain IDs.
"""
function rundssp!(mod::Model)
    # Write the structure to a temporary PDB file
    # Our mmCIF writer does not write out enough of the required entries for DSSP to read it
    tmp_pdb_path = tempname() * ".pdb"
    open(tmp_pdb_path, "w") do io
        println(io, "HEADER")
        writepdb(io, mod)
    end

    # Run DSSP on the temporary PDB file
    dssp_output_lines = readlines(pipeline(`$dssp_executable $tmp_pdb_path`))
    rm(tmp_pdb_path)

    data_begin = false
    for line in dssp_output_lines
        if startswith(line, "  #  RESIDUE AA STRUCTURE BP1")
            data_begin = true
            continue
        end
        data_begin || continue
        line[14] == '!' && continue
        res_id = strip(line[6:11]) # Insertion code is in column 11
        ch_id = line[12]
        ss_code = line[17]
        ch = mod[ch_id]
        # DSSP does not mark hetero atoms
        if haskey(ch.residues, res_id)
            sscode!(ch[res_id], ss_code)
        elseif haskey(ch.residues, "H_$res_id")
            sscode!(ch["H_$res_id"], ss_code)
        else
            error("Could not assign secondary structure to residue ID $res_id in chain $ch_id")
        end
    end
    return mod
end

function rundssp!(struc::ProteinStructure)
    for mod in values(models(struc))
        rundssp!(mod)
    end
    return struc
end

"""
    rundssp(struc)
    rundssp(mod)
    rundssp(filepath_in, dssp_filepath_out)

Return a copy of the structural element with DSSP (Define Secondary Structure of Proteins)
run to assign secondary structure, or run DSSP directly on a PDB or mmCIF file.
"""
rundssp(el::Union{ProteinStructure, Model}) = rundssp!(deepcopy(el))

function rundssp(filepath_in, dssp_filepath_out)
    run(`$dssp_executable $filepath_in $dssp_filepath_out`)
end

"""
    runstride!(struc)
    runstride!(mod)

Run STRIDE on the given structural element to assign secondary structure.

A temporary PDB file is written, so this will fail if the structural element cannot
be written to a PDB file, for example if there are two-letter chain IDs.
"""
function runstride!(mod::Model)
    # Write the structure to a temporary PDB file
    # STRIDE does not work with mmCIF files
    tmp_pdb_path = tempname() * ".pdb"
    open(tmp_pdb_path, "w") do io
        writepdb(io, mod)
    end

    # Run STRIDE on the temporary PDB file
    stride_output_lines = readlines(pipeline(`$stride_executable $tmp_pdb_path`))
    @assert length(stride_output_lines) > 1
    rm(tmp_pdb_path)

    for line in stride_output_lines
        if startswith(line, "ASG")
            ch_id = line[10]
            res_id = strip(line[11:15]) # Insertion code is directly after residue number
            ss_code = uppercase(line[25]) # STRIDE sometimes returns 'b' instead of 'B'
            sscode!(mod[ch_id][res_id], ss_code)
        end
    end
    return mod
end

function runstride!(struc::ProteinStructure)
    for mod in values(models(struc))
        runstride!(mod)
    end
    return struc
end

"""
    runstride(struc)
    runstride(mod)
    runstride(filepath_in, stride_filepath_out)

Return a copy of the structural element with STRIDE
run to assign secondary structure, or run STRIDE directly on a PDB file.
"""
runstride(el::Union{ProteinStructure, Model}) = runstride!(deepcopy(el))

function runstride(filepath_in, stride_filepath_out)
    run(`$stride_executable $filepath_in -f$stride_filepath_out`)
end

"""
    sscode(res)
    sscode(at)

Get the secondary structure code of an `AbstractResidue` or `AbstractAtom` as a `Char`.

`'$ss_code_unassigned'` represents unassigned secondary structure.
Secondary structure can be assigned using `rundssp!` or `runstride!`.
"""
sscode(res::Residue) = res.ss_code
sscode(dis_res::DisorderedResidue) = sscode(defaultresidue(dis_res))
sscode(at::Atom) = sscode(residue(at))
sscode(dis_at::DisorderedAtom) = sscode(defaultatom(dis_at))

"""
    sscode!(res, ss_code)

Set the secondary structure code of an `AbstractResidue` to a `Char`.
"""
function sscode!(res::Residue, ss_code)
    res.ss_code = ss_code
    return res
end

function sscode!(dis_res::DisorderedResidue, ss_code)
    for res_name in resnames(dis_res)
        sscode!(disorderedres(dis_res, res_name), ss_code)
    end
    return dis_res
end

"""
    sscodeselector(res, ss_codes)
    sscodeselector(at, ss_codes)

Determines if an `AbstractResidue` or `AbstractAtom` has its secondary structure code
in a list of secondary structure codes.
"""
function sscodeselector(el::Union{AbstractResidue, AbstractAtom}, ss_codes)
    return sscode(el) in ss_codes
end

"`Set` of secondary structure codes corresponding to an α-helix."
const helixsscodes = Set(['G', 'H', 'I', 'P'])

"`Set` of secondary structure codes corresponding to a β-sheet."
const sheetsscodes = Set(['E', 'B'])

"`Set` of secondary structure codes corresponding to a coil."
const coilsscodes = Set(['C', 'T', 'S', ' '])

"""
    helixselector(res)
    helixselector(at)

Determines if an `AbstractResidue` or `AbstractAtom` is part of an α-helix,
i.e. whether the secondary structure code is in `helixsscodes`.
"""
function helixselector(el::Union{AbstractResidue, AbstractAtom})
    return sscodeselector(el, helixsscodes)
end

"""
    sheetselector(res)
    sheetselector(at)

Determines if an `AbstractResidue` or `AbstractAtom` is part of a β-sheet,
i.e. whether the secondary structure code is in `sheetsscodes`.
"""
function sheetselector(el::Union{AbstractResidue, AbstractAtom})
    return sscodeselector(el, sheetsscodes)
end

"""
    coilselector(res)
    coilselector(at)

Determines if an `AbstractResidue` or `AbstractAtom` is part of a coil,
i.e. whether the secondary structure code is in `coilsscodes`.
"""
function coilselector(el::Union{AbstractResidue, AbstractAtom})
    return sscodeselector(el, coilsscodes)
end
