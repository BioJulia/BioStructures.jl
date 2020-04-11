# BioStructures examples

Further examples of BioStructures usage are given below.

**A)** Plot the temperature factors of a protein:

```julia
using Plots
calphas = collectatoms(struc, calphaselector)
plot(resnumber.(calphas),
     tempfactor.(calphas),
     xlabel="Residue number",
     ylabel="Temperature factor",
     label="")
```

**B)** Print the PDB records for all Cα atoms within 5 Å of residue 38:

```julia
for at in calphas
    if distance(struc['A'][38], at) < 5.0 && resnumber(at) != 38
        println(pdbline(at))
    end
end
```

**C)** Find the residues at the interface of a protein-protein interaction:

```julia
for res_a in collectresidues(struc["A"], standardselector)
    for res_b in collectresidues(struc["B"], standardselector)
        if distance(res_a, res_b) < 5.0
            println(resnumber(res_a), "A ", resnumber(res_b), "B")
        end
    end
end
```

**D)** Show the Ramachandran phi/psi angle plot of a structure:

```julia
using Plots
phi_angles, psi_angles = ramachandranangles(struc, standardselector)
scatter(rad2deg.(phi_angles),
     rad2deg.(psi_angles),
     title="Ramachandran plot",
     xlabel="Phi / degrees",
     ylabel="Psi / degrees",
     label="",
     xticks=[-180, -90, 0, 90, 180],
     yticks=[-180, -90, 0, 90, 180],
     xlims=(-180, 180),
     ylims=(-180, 180))
```

**E)** Calculate the RMSD and displacements between the heavy (non-hydrogen) atoms of two models in an NMR structure:

```julia
downloadpdb("1SSU")
struc_nmr = read("1SSU.pdb", PDB)
rmsd(struc_nmr[5], struc_nmr[10], superimpose=false, rmsdatoms=heavyatomselector)
displacements(struc_nmr[5], struc_nmr[10], superimpose=false, rmsdatoms=heavyatomselector)
```

**F)** Calculate the cysteine fraction of every structure in the PDB:

```julia
l = pdbentrylist()
for p in l
    downloadpdb(p, format=MMCIF) do fp
        s = read(fp, MMCIF)
        nres = countresidues(s, standardselector)
        if nres > 0
            frac = countresidues(s, standardselector, x -> resname(x) == "CYS") / nres
            println(p, "  ", round(frac, digits=2))
        end
    end
end
```

**G)** Interoperability is possible with other packages in the [Julia ecosystem](https://pkg.julialang.org/docs).
For example, use [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) to find the 10 nearest residues to each residue:

```julia
using NearestNeighbors
struc = retrievepdb("1AKE")
ca = coordarray(struc["A"], cbetaselector)
kdtree = KDTree(ca; leafsize=10)
idxs, dists = knn(kdtree, ca, 10, true)
```

**H)** Interoperability with [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) gives access to filtering, sorting, summary statistics and other writing options:

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
