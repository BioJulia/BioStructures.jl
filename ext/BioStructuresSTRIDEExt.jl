module BioStructuresSTRIDEExt

using BioStructures
using STRIDE_jll

const stride_executable = `$(STRIDE_jll.stride_exe())`

function BioStructures.runstride!(mo::Model)
    # Write the structure to a temporary PDB file
    # STRIDE does not work with mmCIF files
    tmp_pdb_path = tempname() * ".pdb"
    open(tmp_pdb_path, "w") do io
        writepdb(io, mo)
    end

    # Run STRIDE on the temporary PDB file
    stride_output_lines = readlines(pipeline(`$stride_executable $tmp_pdb_path`))
    rm(tmp_pdb_path)

    for line in stride_output_lines
        if startswith(line, "ASG")
            ch_id = line[10]
            res_id = strip(line[11:15]) # Insertion code is directly after residue number
            ss_code = uppercase(line[25]) # STRIDE sometimes returns 'b' instead of 'B'
            sscode!(mo[ch_id][res_id], ss_code)
        end
    end
    return mo
end

function BioStructures.runstride!(struc::MolecularStructure)
    for mo in values(models(struc))
        runstride!(mo)
    end
    return struc
end

BioStructures.runstride(el::Union{MolecularStructure, Model}) = runstride!(copy(el))

function BioStructures.runstride(filepath_in, stride_filepath_out)
    run(`$stride_executable $filepath_in -f$stride_filepath_out`)
end

end # BioStructuresSTRIDEExt
