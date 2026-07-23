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
of residues in the chain, as given by [`atombonds`](@ref).

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

    for (ai, aj) in atombonds(chain; strict)
        add_edge!(mg, mg[ai, :element], mg[aj, :element])
    end

    return mg
end

end # BioStructuresGraphsExt
