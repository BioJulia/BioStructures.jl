# BioStructures examples

The best way to learn how to use the package is to read the [BioStructures documentation](@ref).
Here we give further examples, showing what you can do with the package.

**A)** Print the centroid coordinate of the sidechain heavy atoms for each residue in a protein:

```julia
using BioStructures
using Statistics

struc = read("1AKE.pdb", PDBFormat)
res_list = collectresidues(struc, standardselector)

function sidechainheavyselector(a::AbstractAtom)
    return !hydrogenselector(a) && !atomnameselector(a, backboneatomnames)
end

for res in res_list
    print(resid(res, full=true), "  ")
    if resname(res) == "GLY"
        println("no sidechain")
    else
        coord_array = coordarray(res, sidechainheavyselector)
        println(join(mean(coord_array, dims=2), "  "))
    end
end
```

**B)** Plot the temperature factors of a protein:

```julia
using Plots
calphas = collectatoms(struc, calphaselector)
plot(
    resnumber.(calphas),
    tempfactor.(calphas),
    xlabel="Residue number",
    ylabel="Temperature factor",
    label="",
)
```

**C)** Print the PDB records for all Cα atoms within 5 Å of residue 38:

```julia
for at in calphas
    if distance(struc['A'][38], at) < 5.0 && resnumber(at) != 38
        println(pdbline(at))
    end
end
```

**D)** Find the residues at the interface of a protein-protein interaction:

```julia
for res_a in collectresidues(struc["A"], standardselector)
    for res_b in collectresidues(struc["B"], standardselector)
        if distance(res_a, res_b) < 5.0
            println(resnumber(res_a), "A ", resnumber(res_b), "B")
        end
    end
end
```

**E)** Show the Ramachandran phi/psi angle plot of a protein structure:

```julia
using Plots
phi_angles, psi_angles = ramachandranangles(struc, standardselector)
scatter(
    rad2deg.(phi_angles),
    rad2deg.(psi_angles),
    title="Ramachandran plot",
    xlabel="Phi / degrees",
    ylabel="Psi / degrees",
    label="",
    xticks=[-180, -90, 0, 90, 180],
    yticks=[-180, -90, 0, 90, 180],
    xlims=(-180, 180),
    ylims=(-180, 180),
)
```

**F)** Calculate the RMSD and displacements between the heavy (non-hydrogen) atoms of two models in an NMR structure:

```julia
downloadpdb("1SSU")
struc_nmr = read("1SSU.pdb", PDBFormat)
rmsd(struc_nmr[5], struc_nmr[10], superimpose=false, rmsdatoms=heavyatomselector)
displacements(struc_nmr[5], struc_nmr[10], superimpose=false, rmsdatoms=heavyatomselector)
```

**G)** Calculate the cysteine fraction of every protein in the PDB:

```julia
l = pdbentrylist()
for p in l
    downloadpdb(p, format=MMCIFFormat) do fp
        s = read(fp, MMCIFFormat)
        nres = countresidues(s, standardselector)
        if nres > 0
            frac = countresidues(s, standardselector, x -> resname(x) == "CYS") / nres
            println(p, "  ", round(frac, digits=2))
        end
    end
end
```

**H)** Interoperability is possible with other packages in the [Julia ecosystem](https://pkg.julialang.org/docs).
For example, use [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) to find the 10 nearest residues to each residue:

```julia
using NearestNeighbors
struc = retrievepdb("1AKE")
ca = coordarray(struc["A"], cbetaselector)
kdtree = KDTree(ca; leafsize=10)
idxs, dists = knn(kdtree, ca, 10, true)
```

**I)** Interoperability with [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) gives access to filtering, sorting, summary statistics and other writing options:

```julia
using DataFrames
using CSV
using Statistics
struc = retrievepdb("1ALW")
df = DataFrame(collectatoms(struc))
describe(df) # Show summary
mean(df.tempfactor) # Column-wise operations
sort(df, :x) # Sorting
CSV.write("1ALW.csv", df) # CSV file writing
```

**J)** Download a model from the [AlphaFold database](https://alphafold.ebi.ac.uk) and show the secondary structure assignment of each residue.

```julia
using DSSP_jll
fp = download("https://alphafold.ebi.ac.uk/files/AF-P24941-F1-model_v4.pdb")
struc = read(fp, PDBFormat; run_dssp=true)
sscodes = sscode.(collectresidues(struc))
println(join(sscodes))
```
