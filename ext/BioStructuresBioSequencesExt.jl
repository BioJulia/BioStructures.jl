module BioStructuresBioSequencesExt

using BioStructures
using BioSequences

"""
    LongAA(el)

Return the amino acid sequence of a structural element.

Additional arguments are residue selector functions - only residues that
return `true` from all the functions are retained.
The `gaps` keyword argument determines whether to add gaps to the sequence
based on missing residue numbers (default `true`).
See BioSequences.jl for more on how to use sequences.
`LongAA` is an alias for `LongSequence{AminoAcidAlphabet}`.
"""
function BioSequences.LongAA(el::Union{StructuralElement, Vector{Model},
                        Vector{Chain}, Vector{<:AbstractAtom}},
                        residue_selectors::Function...;
                        gaps::Bool=true)
    return LongAA(collectresidues(el, residue_selectors...); gaps=gaps)
end

function BioSequences.LongAA(res::AbstractVector{<:AbstractResidue}; gaps::Bool=true)
    seq = AminoAcid[]
    for i in 1:length(res)
        res_symbol = get(threeletter_to_aa, resname(res[i], strip=false), AA_X)
        push!(seq, res_symbol)
        # Add gaps based on missing residue numbers
        if gaps &&
                i + 1 <= length(res) &&
                resnumber(res[i + 1]) - resnumber(res[i]) > 1 &&
                chainid(res[i]) == chainid(res[i + 1])
            append!(seq, [AA_Gap for _ in 1:(resnumber(res[i + 1]) - resnumber(res[i]) - 1)])
        end
    end
    return LongAA(seq)
end

end # BioStructuresBioSequencesExt
