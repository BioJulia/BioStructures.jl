# Benchmark suite for BioStructures
# Run with something like:
#   using BioStructures
#   using PkgBenchmark
#   results = benchmarkpkg(".")
#   export_markdown(out_file, results)

using BioStructures
using BenchmarkTools
import BioCore # Imported to avoid clash with BioGenerics distance
using Random

# Use files in BioFmtSpecimens
fmtdir = BioCore.Testing.get_bio_fmt_specimens()
testfilepath(path::AbstractString...) = joinpath(fmtdir, path...)

# All writing is done to one temporary file
temp_filename, io = mktemp()
close(io)

pdbids = ["1AKE", "1EN2", "1SSU"]
formats = Dict("PDB"=> PDB, "mmCIF"=> MMCIF, "MMTF"=> MMTF)
writefunctions = Dict("PDB"=> writepdb, "mmCIF"=> writemmcif, "MMTF"=> writemmtf)

const SUITE = BenchmarkGroup(
        [],
        "read"   => BenchmarkGroup([], [f=> BenchmarkGroup() for f in keys(formats)]...),
        "write"  => BenchmarkGroup([], [f=> BenchmarkGroup() for f in keys(formats)]...),
        "dict"   => BenchmarkGroup(),
        "model"  => BenchmarkGroup(),
        "collect"=> BenchmarkGroup(),
        "spatial"=> BenchmarkGroup(),
)

struc = Dict{String, MolecularStructure}()
for pdbid in pdbids
    struc[pdbid] = read(testfilepath("PDB", "$pdbid.pdb"), PDB)
end

for pdbid in pdbids
    for f in keys(formats)
        SUITE["read"][f][pdbid] = @benchmarkable read(
                $(testfilepath(f, "$pdbid.$(pdbextension[formats[f]])")), $(formats[f]))
        SUITE["write"][f][pdbid] = @benchmarkable $(writefunctions[f])(
                $temp_filename, $(struc[pdbid])) teardown=(rm(temp_filename, force=true))
    end
end

SUITE["dict"]["mmCIF"] = @benchmarkable MMCIFDict($(testfilepath("mmCIF", "1AKE.cif" )))
SUITE["dict"]["MMTF" ] = @benchmarkable MMTFDict( $(testfilepath("MMTF" , "1AKE.mmtf")))

SUITE["model"]["atomaccess"  ] = @benchmarkable $(struc["1EN2"])["A"][20]["CA"]
SUITE["model"]["sortatoms"   ] = @benchmarkable sort(ats) setup=(ats = shuffle(collectatoms(   struc["1EN2"])))
SUITE["model"]["sortresidues"] = @benchmarkable sort(res) setup=(res = shuffle(collectresidues(struc["1EN2"])))
SUITE["model"]["iterate"     ] = @benchmarkable for mo in $(struc["1AKE"]) for ch in mo for res in ch for at in res end end end end
SUITE["model"]["sequence"    ] = @benchmarkable LongAA($(collectresidues(struc["1EN2"])))
SUITE["model"]["pdbline"     ] = @benchmarkable pdbline($(struc["1EN2"]["A"][20]["CA"]))

SUITE["collect"]["atoms"       ] = @benchmarkable collectatoms(   $(struc["1EN2"]))
SUITE["collect"]["atomssel"    ] = @benchmarkable collectatoms(   $(struc["1EN2"]), calphaselector)
SUITE["collect"]["atomsdis"    ] = @benchmarkable collectatoms(   $(struc["1EN2"]), expand_disordered=true)
SUITE["collect"]["residues"    ] = @benchmarkable collectresidues($(struc["1EN2"]))
SUITE["collect"]["residuessel" ] = @benchmarkable collectresidues($(struc["1EN2"]), standardselector)
SUITE["collect"]["residuesdis" ] = @benchmarkable collectresidues($(struc["1EN2"]), expand_disordered=true)
SUITE["collect"]["chains"      ] = @benchmarkable collectchains(  $(struc["1AKE"]))
SUITE["collect"]["models"      ] = @benchmarkable collectmodels(  $(struc["1SSU"]))

SUITE["spatial"]["coordarray"    ] = @benchmarkable coordarray($(collectatoms(struc["1AKE"])))
SUITE["spatial"]["transformation"] = @benchmarkable Transformation($(collectresidues(struc["1SSU"][5], standardselector)), $(collectresidues(struc["1SSU"][10], standardselector)))
SUITE["spatial"]["rmsd"          ] = @benchmarkable rmsd($(coordarray(struc["1SSU"][5], heavyatomselector)), $(coordarray(struc["1SSU"][10], heavyatomselector)))
SUITE["spatial"]["distance"      ] = @benchmarkable distance($(struc["1AKE"]["A"][50]), $(struc["1AKE"]["A"][60]))
SUITE["spatial"]["bondangle"     ] = @benchmarkable bondangle($(struc["1AKE"]["A"][20]["N"]), $(struc["1AKE"]["A"][20]["CA"]), $(struc["1AKE"]["A"][20]["C"]))
SUITE["spatial"]["dihedralangle" ] = @benchmarkable dihedralangle($(struc["1AKE"]["A"][20]["N"]), $(struc["1AKE"]["A"][20]["CA"]), $(struc["1AKE"]["A"][20]["C"]), $(struc["1AKE"]["A"][21]["N"]))
SUITE["spatial"]["ramachandran"  ] = @benchmarkable ramachandranangles($(collectresidues(struc["1AKE"]["A"], standardselector)))
SUITE["spatial"]["contactmap"    ] = @benchmarkable ContactMap($(collectatoms(struc["1AKE"]["A"], cbetaselector)), 8.0)
SUITE["spatial"]["distancemap"   ] = @benchmarkable DistanceMap($(collectresidues(struc["1AKE"]["A"], standardselector)))
SUITE["spatial"]["contactgraph"  ] = @benchmarkable MetaGraph($(collectatoms(struc["1AKE"]["A"], cbetaselector)), 8.0)
