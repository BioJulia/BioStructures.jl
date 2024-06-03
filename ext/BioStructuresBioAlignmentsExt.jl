module BioStructuresBioAlignmentsExt

using BioStructures
using BioAlignments
using BioSequences

"""
    pairalign(el1, el2, residue_selectors...)

Carries out a pairwise sequence alignment between the sequences of two
structural elements.

Additional arguments are residue selector functions - only residues that return
`true` from all the functions are retained.
The keyword arguments `scoremodel` (default
`AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1)`) and `aligntype`
(default `GlobalAlignment()`) determine the properties of the alignment.
See BioAlignments.jl for more on how to use alignments.
"""
function BioAlignments.pairalign(el1::StructuralElementOrList,
                            el2::StructuralElementOrList,
                            residue_selectors::Function...;
                            scoremodel::AbstractScoreModel=AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1),
                            aligntype::BioAlignments.AbstractAlignment=GlobalAlignment())
    seq1 = LongAA(el1, residue_selectors...; gaps=false)
    seq2 = LongAA(el2, residue_selectors...; gaps=false)
    return pairalign(aligntype, seq1, seq2, scoremodel)
end

function BioStructures.Transformation(el1::StructuralElementOrList,
                el2::StructuralElementOrList,
                residue_selectors::Function...;
                scoremodel::AbstractScoreModel=AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1),
                aligntype::BioAlignments.AbstractAlignment=LocalAlignment(),
                alignatoms::Function=calphaselector)
    res1 = collectresidues(el1, residue_selectors...)
    res2 = collectresidues(el2, residue_selectors...)
    # Shortcut if the sequences are the same
    if LongAA(res1; gaps=false) == LongAA(res2; gaps=false)
        inds1 = collect(1:length(res1))
        inds2 = collect(1:length(res2))
    else
        alres = pairalign(res1, res2; scoremodel=scoremodel, aligntype=aligntype)
        al = alignment(alres)
        inds1, inds2 = Int[], Int[]
        # Offset residue counter based on start of aligned region
        first_anchor = first(al.a.aln.anchors)
        counter1 = first_anchor.seqpos
        counter2 = first_anchor.refpos
        # Obtain indices of residues used in alignment
        for (v1, v2) in al
            if v1 != AA_Gap
                counter1 += 1
            end
            if v2 != AA_Gap
                counter2 += 1
            end
            if v1 != AA_Gap && v2 != AA_Gap
                push!(inds1, counter1)
                push!(inds2, counter2)
            end
        end
    end
    @info "Superimposing based on a sequence alignment between $(length(inds1)) residues"
    atoms1, atoms2 = AbstractAtom[], AbstractAtom[]
    for (i1, i2) in zip(inds1, inds2)
        sel_ats1 = collectatoms(res1[i1], alignatoms)
        sel_ats2 = collectatoms(res2[i2], alignatoms)
        if length(atoms1) == length(atoms2)
            append!(atoms1, sel_ats1)
            append!(atoms2, sel_ats2)
        end
    end
    if length(atoms1) == 0
        throw(ArgumentError("No atoms found to superimpose"))
    end
    @info "Superimposing based on $(length(atoms1)) atoms"
    return Transformation(coordarray(atoms1), coordarray(atoms2), inds1, inds2)
end

end # BioStructuresBioAlignmentsExt
