module BioStructuresDSSPExt

using BioStructures
using DSSP_jll

const dssp_executable = `$(DSSP_jll.mkdssp()) --mmcif-dictionary $(DSSP_jll.mmcif_pdbx_dic) --output-format dssp`

function BioStructures.rundssp!(mo::Model)
    # Write the structure to a temporary PDB file
    # Our mmCIF writer does not write out enough of the required entries for DSSP to read it
    tmp_pdb_path = tempname() * ".pdb"
    open(tmp_pdb_path, "w") do io
        println(io, "HEADER")
        writepdb(io, mo)
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
        ch = mo[ch_id]
        # DSSP does not mark hetero atoms
        if haskey(ch.residues, res_id)
            sscode!(ch[res_id], ss_code)
        elseif haskey(ch.residues, "H_$res_id")
            sscode!(ch["H_$res_id"], ss_code)
        else
            error("Could not assign secondary structure to residue ID $res_id in chain $ch_id")
        end
    end
    return mo
end

function BioStructures.rundssp!(struc::MolecularStructure)
    for mo in values(models(struc))
        rundssp!(mo)
    end
    return struc
end

BioStructures.rundssp(el::Union{MolecularStructure, Model}) = rundssp!(copy(el))

function BioStructures.rundssp(filepath_in, dssp_filepath_out)
    run(`$dssp_executable $filepath_in $dssp_filepath_out`)
end

end # BioStructuresDSSPExt
