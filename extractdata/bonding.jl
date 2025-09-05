# Running this file generates src/bonding.jl, which contains the atomtypes and residuedata dictionaries.
# Thanks to OpenMM for the ff14SB force field XML file.
# van der Waals radii are from
#    Mantina, Manjeera, et al. "Consistent van der Waals radii for the whole main group." The Journal of Physical Chemistry A 113.19 (2009): 5806-5812. (Table 12)

using Downloads

if !isfile(joinpath(@__DIR__, "protein.ff14SB.xml"))
    Downloads.download("https://raw.githubusercontent.com/openmm/openmm/refs/heads/master/wrappers/python/openmm/app/data/amber14/protein.ff14SB.xml", "protein.ff14SB.xml")
end
if !isfile(joinpath(@__DIR__, "vanderWaals_radii.csv"))
    Downloads.download("https://raw.githubusercontent.com/mojaie/MolecularGraph.jl/refs/heads/master/assets/const/vanderWaals_radii.csv", "vanderWaals_radii.csv")
    Downloads.download("https://raw.githubusercontent.com/mojaie/MolecularGraph.jl/refs/heads/master/assets/const/symboltonumber.yaml", "symboltonumber.yaml")
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

atomtypes, residues, bondlengths, bondangles = open("protein.ff14SB.xml", "r") do io
    line = readline(io)
    @assert line == "<ForceField>"
    atomtypes = Dict{String, @NamedTuple{element::String, mass::Float32, name::String}}()
    residues = Dict{String, @NamedTuple{atoms::Dict{String, @NamedTuple{charge::Float32, type::String}}, bonds::Vector{Tuple{String,String}}, externalbonds::Vector{String}}}()
    harmonicbonds = Dict{Tuple{String,String}, @NamedTuple{length::Float32}}()
    harmonicangles = Dict{Tuple{String,String,String}, @NamedTuple{angle::Float32}}()
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
        end
    end
    atomtypes, residues, Dict(k => 10 * v.length for (k, v) in harmonicbonds), Dict(k => v.angle for (k, v) in harmonicangles)
end

num2sym = open("symboltonumber.yaml", "r") do io
    d = Dict{Int, String}()
    while !eof(io)
        line = strip(readline(io))
        el, n = chomp.(split(line, ':'))
        n = parse(Int, n)
        n == 1 && el != "H" && continue # skip Deuterium and Tritium
        startswith(el, "Uu") && continue # skip three-letter element codes
        d[n] = el
    end
    return d
end

dequote(s) = startswith(s, '"') ? String(s[2:end-1]) : s

vdw = open("vanderWaals_radii.csv", "r") do io
    line = strip(readline(io))  # skip comments and header
    while line[1] == '#'
        line = strip(readline(io))
    end
    d = Dict{String, Tuple{Int, Float32}}()
    while !eof(io)
        line = strip(readline(io))
        el, r = dequote.(split(line, '\t'))
        el = parse(Int, el)
        d[num2sym[el]] = (el, parse(Float32, r))
    end
    return d
end

open(joinpath(dirname(@__DIR__), "src", "bonding.jl"), "w") do io
    println(io, "# This file is autogenerated by extractdata/bonding.jl")
    println(io, "# Do not edit this file directly.")
    println(io)

    println(io, "const atomtypes = Dict{String, @NamedTuple{element::String, mass::Float32, name::String}}(")
    at = sort!(collect(atomtypes); by=first)
    for pr in at
        println(io, "    ", pr, ',')
    end
    println(io, ")\n")

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

    println(io, "\nconst vanderWaals_radii = Dict{String, Float32}(")
    vdw_sorted = sort!(collect(vdw); by=pr->pr.second[1])
    for (k, v) in vdw_sorted
        println(io, "    ", k => v[2], ',')
    end
    println(io, ")")
end
