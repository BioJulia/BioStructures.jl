module BioStructuresGraphsExt

using BioStructures
using Graphs
using MetaGraphs
using BioStructures: findatombyname

"""
    MetaGraph(element, contact_distance)

Construct a graph of elements where edges are contacts separated by less
than `contact_distance`.

See Graphs.jl and MetaGraphs.jl for more on how to use graphs.
"""
function MetaGraphs.MetaGraph(el::StructuralElementOrList, contact_dist::Real)
    sq_contact_dist = contact_dist ^ 2
    el_list = collect(el)
    mg = MetaGraph(length(el_list))
    set_prop!(mg, :contactdist, Float64(contact_dist))
    for (i, el) in enumerate(el_list)
        set_prop!(mg, i, :element, el)
        for j in 1:(i - 1)
            if sqdistance(el, el_list[j]) <= sq_contact_dist
                add_edge!(mg, j, i)
            end
        end
    end
    set_indexing_prop!(mg, :element)
    return mg
end

"""
    MetaGraph(chain::Chain; strict::Bool=true)

Construct a graph of atoms where edges are determined by the known bonds
of residues in the chain.

By default, the graph is constructed in `strict` mode, which means that:

- residue and atom names must be standard
- all hydrogens are present
- HIS must be disambiguated as HIE, HID, or HIP

These constraints can be relaxed by setting `strict = false`, at some
risk to accuracy.

See Graphs.jl and MetaGraphs.jl for more on how to use graphs.
"""
function MetaGraphs.MetaGraph(chain::Chain; strict::Bool=true)
    el_list = collectatoms(chain; expand_disordered=true)
    mg = MetaGraph(length(el_list))
    for (i, el) in enumerate(el_list)
        set_prop!(mg, i, :element, el)
    end
    set_indexing_prop!(mg, :element)

    prev = nothing
    for r in chain
        ishetero(r) && continue
        if prev !== nothing
            if resnumber(r) == resnumber(prev) + 1
                # Add the peptide bond(s) (disordered residues may need multiples)
                for _r in collectresidues(r; expand_disordered=true), _rp in collectresidues(prev; expand_disordered=true)
                    for a in _r["N"], ap in _rp["C"]
                        add_edge!(mg, mg[a, :element], mg[ap, :element])
                    end
                end
            else
                prev = nothing
            end
        end
        # Add the residue bonds
        addresiduebonds!(mg, r, strict)
        prev = r
        # The "OXT" (C-terminus oxygen) is not connected
        for _r in collectresidues(r; expand_disordered=true)
            aj = findatombyname(_r, "OXT"; strict=false)
            if aj !== nothing
                ai = findatombyname(_r, "C"; strict)
                connect_atoms!(mg, ai, aj)
            end
        end
    end

    return mg
end

function addresiduebonds!(mg, r::Residue, strict::Bool)
    rname = residuekey(r, strict)
    strict || haskey(BioStructures.residuedata, rname) || return
    rd = BioStructures.residuedata[rname]
    for (ni, nj) in rd.bonds
        connect_atoms!(mg, r, ni, nj, strict)
    end
end

connect_atoms!(mg, r, ni::AbstractString, nj::AbstractString, strict::Bool) =
    connect_atoms!(mg, findatombyname(r, ni; strict), findatombyname(r, nj; strict))

function connect_atoms!(mg, ai::Union{AbstractAtom,Nothing}, aj::Union{AbstractAtom,Nothing})
    (ai === nothing || aj === nothing) && return
    for _ai in ai, _aj in aj  # handle disordered atoms
        i = mg[_ai, :element]
        j = mg[_aj, :element]
        add_edge!(mg, i, j)
    end
end

addresiduebonds!(mg, dr::DisorderedResidue, strict::Bool) = for r in values(dr.names)
    addresiduebonds!(mg, r, strict)
end

function residuekey(r::Residue, strict::Bool)
    rname = resname(r)
    if rname == "HIS"
        strict && error("HIS is not a valid residue name in strict mode; use HIE, HID, or HIP")
        rname = "HIE"
    end
    return rname
end

end # BioStructuresGraphsExt
