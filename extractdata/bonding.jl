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

function parsexmlline(f, line, tag, keyname)
    @assert startswith(line, "<$tag ")
    @assert endswith(line, "/>")
    kv = split(strip(line[length(tag)+2:end-2]), ' ')
    key = ""
    vals = Pair{Symbol,Any}[]
    for kvp in kv
        k, v = split(kvp, '=')
        if k == keyname
            key = parsestring(v)
        else
            push!(vals, Symbol(k) => f(k, v))
        end
    end
    return key => (; vals...)
end


atomtypes, residues = open("protein.ff14SB.xml", "r") do io
    line = readline(io)
    @assert line == "<ForceField>"
    atomtypes = Dict{String, @NamedTuple{element::String, mass::Float32, name::String}}()
    residues = Dict{String, @NamedTuple{atoms::Dict{String, @NamedTuple{charge::Float32, type::String}}, bonds::Vector{Tuple{String,String}}, externalbonds::Vector{String}}}()
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
        end
    end
    atomtypes, residues
end

open(joinpath(dirname(@__DIR__), "src", "bonding.jl"), "w") do io
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
end
