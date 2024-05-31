module BioStructuresDataFramesExt

using BioStructures
using DataFrames

"""
    DataFrame(atom_list, atom_selectors...)
    DataFrame(residue_list, residue_selectors...)

Construct a `DataFrame` from a list of atoms or residues.

Additional arguments are selector functions - only atoms or residues that return
`true` from all the functions are retained.
The keyword argument `expand_disordered` (default `true`) determines whether to
return all copies of disordered atoms or residues separately.
See DataFrames.jl for more on how to use `DataFrame`s.
"""
function DataFrames.DataFrame(ats::AbstractVector{<:AbstractAtom},
                    atom_selectors::Function...;
                    expand_disordered::Bool=true)
    df = DataFrame(ishetero=Bool[],
                    serial=Int[],
                    atomname=String[],
                    altlocid=Char[],
                    resname=String[],
                    chainid=String[],
                    resnumber=Int[],
                    inscode=Char[],
                    x=Float64[],
                    y=Float64[],
                    z=Float64[],
                    occupancy=Float64[],
                    tempfactor=Float64[],
                    element=String[],
                    charge=String[],
                    modelnumber=Int[],
                    isdisorderedatom=Bool[])
    for a in collectatoms(ats, atom_selectors...;
                            expand_disordered=expand_disordered)
        push!(df, (ishetero(a), serial(a), atomname(a), altlocid(a), resname(a),
                    chainid(a), resnumber(a), inscode(a), BioStructures.x(a),
                    BioStructures.y(a), BioStructures.z(a), occupancy(a), tempfactor(a),
                    element(a), charge(a), modelnumber(a), isdisorderedatom(a)))
    end
    return df
end

function DataFrames.DataFrame(res::AbstractVector{<:AbstractResidue},
                    residue_selectors::Function...;
                    expand_disordered::Bool=true)
    df = DataFrame(ishetero=Bool[],
                    resname=String[],
                    chainid=String[],
                    resnumber=Int[],
                    inscode=Char[],
                    countatoms=Int[],
                    modelnumber=Int[],
                    isdisorderedres=Bool[])
    for r in collectresidues(res, residue_selectors...;
                                expand_disordered=expand_disordered)
        push!(df, (ishetero(r), resname(r), chainid(r), resnumber(r),
                    inscode(r), countatoms(r; expand_disordered=expand_disordered),
                    modelnumber(r), isdisorderedres(r)))
    end
    return df
end

end # BioStructuresDataFramesExt
