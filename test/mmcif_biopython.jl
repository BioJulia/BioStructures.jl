# Script to check that all mmCIF dictionary entries are the same as those from Biopython
# Requires a Python installation with Biopython and NumPy installed

using BioStructures
using PyCall

version = pyimport("Bio").__version__
@info "Biopython version is $version"

MMCIF2Dict_python = pyimport("Bio.PDB.MMCIF2Dict").MMCIF2Dict

pdbids = pdbentrylist()

for pdbid in pdbids
    filepath = downloadpdb(pdbid, format=MMCIF)
    if isfile(filepath)
        dict_julia = MMCIFDict(filepath).dict
        dict_python = MMCIF2Dict_python(filepath)
        rm(filepath)
        status = "identical"
        if convert(Array{String, 1}, sort(collect(keys(dict_python)))) == sort(collect(keys(dict_julia)))
            for key in keys(dict_julia)
                if key == "data_"
                    value_python = [dict_python[key]]
                else
                    value_python = dict_python[key]
                end
                if dict_julia[key] != value_python
                    status = "different $key"
                    break
                end
            end
        else
            status = "different keys"
        end
    else
        status = "download_failed"
    end
    println(pdbid, " ", status)
end
