# Running this file generates src/bonding.jl, which contains the atomtypes and residuedata dictionaries.
# Thanks to OpenMM for the ff14SB force field XML file.

using Downloads

if !isfile(joinpath(@__DIR__, "protein.ff14SB.xml"))
    Downloads.download("https://raw.githubusercontent.com/openmm/openmm/refs/heads/master/wrappers/python/openmm/app/data/amber14/protein.ff14SB.xml", "protein.ff14SB.xml")
end

function parsexmblock(f, io::IO, key)
    while !eof(io)
        line = strip(readline(io))
        line == key && return nothing
        f(line)
    end
end

function parsestring(str)
    @assert startswith(str, '"')
    @assert endswith(str, '"')
    return String(str[2:end-1])
end

function parsexmlline(f, line, tag, keynames...; skip=())
    @assert startswith(line, "<$tag ")
    @assert endswith(line, "/>")
    kv = split(strip(line[length(tag)+2:end-2]), ' ')
    key = String[]
    vals = Pair{Symbol,Any}[]
    for kvp in kv
        k, v = split(kvp, '=')
        k ∈ skip && continue
        if k ∈ keynames
            push!(key, parsestring(v))
        else
            push!(vals, Symbol(k) => f(k, v))
        end
    end
    @assert length(key) == length(keynames)
    key = length(key) == 1 ? only(key) : (key...,)
    return key => (; vals...)
end

rawatomtypes, residues, bondlengths, bondangles, nonbonded, scale14 = open("protein.ff14SB.xml", "r") do io
    line = readline(io)
    @assert line == "<ForceField>"
    atomtypes = Dict{String, @NamedTuple{element::String, mass::Float32, name::String}}()
    residues = Dict{String, @NamedTuple{atoms::Dict{String, @NamedTuple{charge::Float32, type::String}}, bonds::Vector{Tuple{String,String}}, externalbonds::Vector{String}}}()
    harmonicbonds = Dict{Tuple{String,String}, @NamedTuple{length::Float32}}()
    harmonicangles = Dict{Tuple{String,String,String}, @NamedTuple{angle::Float32}}()
    # The XML is OpenMM's, so lengths are nm and energies kJ/mol; convert to Å and kcal/mol.
    nonbonded = Dict{String, @NamedTuple{epsilon::Float64, sigma::Float64}}()
    scale14 = nothing
    parsexmblock(io, "</ForceField>") do line
        if line == "<AtomTypes>"
            parsexmblock(io, "</AtomTypes>") do line
                push!(atomtypes, parsexmlline(line, "Type", "class") do k, v
                    if k == "element"
                        return parsestring(v)
                    elseif k == "mass"
                        return parse(Float32, v[2:end-1])  # strip the quotes
                    elseif k == "name"
                        return parsestring(v)
                    else
                        error("Unknown AtomType key $k")
                    end
                end)
            end
        elseif line == "<Residues>"
            parsexmblock(io, "</Residues>") do line
                if startswith(line, "<Residue name=")
                    resname = parsestring(line[15:end-1])
                    atoms = Dict{String, @NamedTuple{charge::Float32, type::String}}()
                    bonds = Vector{Tuple{String,String}}()
                    externalbonds = Vector{String}()
                    parsexmblock(io, "</Residue>") do line
                        if startswith(line, "<Atom")
                            push!(atoms, parsexmlline(line, "Atom", "name") do k, v
                                if k == "charge"
                                    return parse(Float32, v[2:end-1])  # strip the quotes
                                elseif k == "type"
                                    return parsestring(v)
                                else
                                    error("Unknown Atom key $k")
                                end
                            end)
                        elseif startswith(line, "<Bond")
                            line = line[6:end-2]
                            a1, a2 = split(strip(line), ' ')
                            push!(bonds, (only(match(r"atomName1=\"(.*)\"", a1).captures), only(match(r"atomName2=\"(.*)\"", a2).captures)))
                        elseif startswith(line, "<ExternalBond")
                            line = line[14:end-2]
                            push!(externalbonds, only(match(r"atomName=\"(.*)\"", line).captures))
                        else
                            error("Unknown Residue line $line")
                        end
                    end
                    residues[resname] = (; atoms, bonds, externalbonds)
                else
                    error("Unknown Residues line $line")
                end
            end
        elseif line == "<HarmonicBondForce>"
            parsexmblock(io, "</HarmonicBondForce>") do line
                push!(harmonicbonds, parsexmlline(line, "Bond", "type1", "type2"; skip=("k",)) do k, v
                    if k == "length"
                        return parse(Float32, v[2:end-1])
                    else
                        error("Unknown Bond key $k")
                    end
                end)
            end
        elseif line == "<HarmonicAngleForce>"
            parsexmblock(io, "</HarmonicAngleForce>") do line
                push!(harmonicangles, parsexmlline(line, "Angle", "type1", "type2", "type3"; skip=("k",)) do k, v
                    if k == "k"
                        return parse(Float32, v[2:end-1])  # strip the quotes
                    elseif k == "angle"
                        return parse(Float32, v[2:end-1])  # strip the quotes
                    else
                        error("Unknown Angle key $k")
                    end
                end)
            end
        elseif startswith(line, "<NonbondedForce ")
            m = match(r"^<NonbondedForce coulomb14scale=\"([^\"]*)\" lj14scale=\"([^\"]*)\">$", line)
            m === nothing && error("Unrecognized NonbondedForce header $line")
            scale14 = (coulomb = parse(Float32, m.captures[1]), lj = parse(Float32, m.captures[2]))
            parsexmblock(io, "</NonbondedForce>") do line
                # <UseAttributeFromResidue name="charge"/> declares that charges come from
                # the residue templates, which is where residuedata already takes them from.
                startswith(line, "<UseAttributeFromResidue ") && return nothing
                push!(nonbonded, parsexmlline(line, "Atom", "type") do k, v
                    if k == "sigma"
                        return 10 * parse(Float64, v[2:end-1])  # nm -> Å
                    elseif k == "epsilon"
                        return parse(Float64, v[2:end-1]) / 4.184  # kJ/mol -> kcal/mol
                    else
                        error("Unknown nonbonded Atom key $k")
                    end
                end)
            end
        end
    end
    atomtypes, residues, Dict(k => 10 * v.length for (k, v) in harmonicbonds), Dict(k => v.angle for (k, v) in harmonicangles), nonbonded, scale14
end

scale14 === nothing && error("no NonbondedForce block found")
let extra = setdiff(keys(nonbonded), (v.name for v in values(rawatomtypes)))
    isempty(extra) || error("nonbonded parameters given for unknown types: $(sort!(collect(extra)))")
end
atomtypes = Dict{String, @NamedTuple{element::String, mass::Float32, name::String, sigma::Float32, epsilon::Float32}}(
    class => begin
        lj = get(() -> error("no nonbonded parameters for atom type $class ($(v.name))"), nonbonded, v.name)
        (; v.element, v.mass, v.name, sigma = Float32(lj.sigma), epsilon = Float32(lj.epsilon))
    end for (class, v) in rawatomtypes)

const atomtypesdoc = raw"""
    atomtypes

The atom types of the Amber ff14SB force field, keyed by type class (e.g., `atomtypes["CT"]`).
Each value is a `NamedTuple` with fields:

- `element`: the chemical element symbol.
- `mass`: the atomic mass in amu.
- `name`: the force field's full name for the type, e.g., `"protein-CT"`. This is the value
  stored in the `type` field of `residuedata`'s per-atom entries.
- `sigma`: the Lennard-Jones distance parameter in Å, the separation at which the pair
  potential crosses zero. Amber's published parameter tables instead list `Rmin/2`, where
  `Rmin = 2^(1/6) * sigma` is the separation of minimum energy.
- `epsilon`: the Lennard-Jones well depth in kcal/mol. Type `"HO"` has `epsilon = 0` and a
  placeholder `sigma`: it carries no Lennard-Jones interaction.

Partial charges are not properties of the type, as they vary between residues; they are in
`residuedata`.

See also `ff14SB_scale14` for the scaling of nonbonded interactions between close neighbors.
"""

const scale14doc = raw"""
    ff14SB_scale14

The factors by which the ff14SB force field scales nonbonded interactions between atoms
separated by exactly three bonds (1-4 pairs): `coulomb` for electrostatics and `lj` for the
Lennard-Jones terms of `atomtypes`.
"""

printdoc(io, doc) = (println(io, "\"\"\""); print(io, doc); println(io, "\"\"\""))

open(joinpath(dirname(@__DIR__), "src", "bonding.jl"), "w") do io
    println(io, "# This file is auto-generated by extractdata/bonding.jl; do not edit directly.")
    println(io, "# It defines only data tables (atomtypes, residuedata, bondlengths, bondangles);")
    println(io, "# code that uses them belongs in a hand-written file such as src/atombonds.jl.\n")
    printdoc(io, atomtypesdoc)
    println(io, "const atomtypes = Dict{String, @NamedTuple{element::String, mass::Float32, name::String, sigma::Float32, epsilon::Float32}}(")
    at = sort!(collect(atomtypes); by=first)
    for pr in at
        println(io, "    ", pr, ',')
    end
    println(io, ")\n")

    printdoc(io, scale14doc)
    println(io, "const ff14SB_scale14 = ", scale14, "\n")

    println(io, "const RDADict = Dict{String, @NamedTuple{charge::Float32, type::String}}")

    println(io, "const residuedata = Dict{String, @NamedTuple{atoms::RDADict, bonds::Vector{Tuple{String,String}}, externalbonds::Vector{String}}}(")
    rd = sort!(collect(residues); by=first)
    for (k, v) in rd
        print(io, "    "*" "^(4-length(k)))
        show(io, k)
        println(io, " => (atoms = ", replace(sprint(show, v.atoms), "Dict{String, @NamedTuple{charge::Float32, type::String}}" => "RDADict"), ',')
        println(io, "               bonds = ", v.bonds, ",")
        println(io, "               externalbonds = ", v.externalbonds, "),")
    end
    println(io, ")")

    println(io, "\nconst bondlengths = Dict{Tuple{String,String}, Float32}(")
    hb = sort!(collect(bondlengths); by=first)
    for pr in hb
        println(io, "    ", pr, ',')
    end
    println(io, ")")

    println(io, "\nconst bondangles = Dict{Tuple{String,String,String}, Float32}(")
    ha = sort!(collect(bondangles); by=first)
    for pr in ha
        println(io, "    ", pr, ',')
    end
    println(io, ")")
end
