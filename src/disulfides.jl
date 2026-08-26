export
    disulfidebonds,
    renamedisulfides!

"""
    disulfidebonds(el; cutoff=2.5)

Return pairs of cysteine SG atoms within `cutoff` Å in `el`.

Residues named `CYS`, `NCYS`, or `CCYS` are considered. An error is thrown
if an SG atom has more than one possible partner.

`el` can be a `Chain`, `Model`, `MolecularStructure`, or any other
`StructuralElementOrList`. Pass a `Model` or `MolecularStructure` to find
disulfides between chains.
"""
function disulfidebonds(el::StructuralElementOrList; cutoff::Real=2.5)
    sgs = Tuple{Residue,AbstractAtom}[]
    for r in collectresidues(el; expand_disordered=true)
        endswith(resname(r), "CYS") || continue
        at = findatombyname(r, "SG"; strict=false)
        at === nothing || push!(sgs, (r, at))
    end
    pairs = Tuple{AbstractAtom,AbstractAtom}[]
    for i in eachindex(sgs)
        partners = [j for j in eachindex(sgs) if j != i && distance(sgs[i][2], sgs[j][2]) <= cutoff]
        if length(partners) > 1
            r = sgs[i][1]
            error("SG of $(resname(r)) $(resnumber(r)) in chain $(chainid(r)) is within $cutoff Å of $(length(partners)) other SG atoms; disulfide pairing is ambiguous")
        elseif length(partners) == 1
            j = only(partners)
            i < j && push!(pairs, (sgs[i][2], sgs[j][2]))
        end
    end
    return pairs
end

cyxname(rname::AbstractString) = rname[1:end-3] * "CYX"

function renamemap(bonds)
    torename = Dict{Residue,String}()
    for (ai, aj) in bonds, at in (ai, aj)
        r = residue(at)
        torename[r] = cyxname(resname(r))
    end
    return torename
end

# Keep DisorderedResidue name lookups consistent with renamed residues.
function applyrenames!(c::Chain, torename::Dict{Residue,String})
    for (key, res) in c.residues
        if isa(res, DisorderedResidue)
            for r in values(res.names)
                haskey(torename, r) && (r.name = torename[r])
            end
            names = Dict(resname(r, strip=false) => r for r in values(res.names))
            c.residues[key] = DisorderedResidue(names, resname(res.names[defaultresname(res)]))
        elseif haskey(torename, res)
            res.name = torename[res]
        end
    end
    return c
end

"""
    renamedisulfides!(el; cutoff=2.5)

Rename disulfide-bonded cysteines to the corresponding Amber `CYX` name.

Partners are identified by [`disulfidebonds`](@ref). `CYS`, `NCYS`, and
`CCYS` become `CYX`, `NCYX`, and `CCYX`, respectively.

`el` can be a `Chain`, `Model` or `MolecularStructure`. A `Chain` argument
only pairs cysteines within that chain. Each model in a
`MolecularStructure` is processed separately.

Call this before [`specializeresnames!`](@ref) when also renaming terminal
residues.
"""
function renamedisulfides!(c::Chain; cutoff::Real=2.5)
    applyrenames!(c, renamemap(disulfidebonds(c; cutoff)))
    return c
end

function renamedisulfides!(m::Model; cutoff::Real=2.5)
    torename = renamemap(disulfidebonds(m; cutoff))
    foreach(c -> applyrenames!(c, torename), m)
    return m
end

renamedisulfides!(s::MolecularStructure; cutoff::Real=2.5) = (foreach(m -> renamedisulfides!(m; cutoff), s); s)
