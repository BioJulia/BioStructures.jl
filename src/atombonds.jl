export atombonds

# Bond derivation
# ===============
#
# Which atoms are bonded is determined by `residuedata`, and so depends only on
# residue and atom names — not on coordinates. Deriving bonds from interatomic
# distances instead would need element-dependent cutoffs (S-H is ~1.34 Å against
# C-H's ~1.09 Å), would still have to special-case geminal hydrogen pairs, and
# would give a different answer for a distorted or clashing conformation.
#
# This lives outside the Graphs extension so that callers who want connectivity
# but not a graph type — force fields, internal-coordinate parametrizations —
# need not depend on Graphs.jl and MetaGraphs.jl. `residuedata` itself is in the
# generated src/bonding.jl; the derivation logic belongs in a hand-written file.

"""
    atombonds(chain::Chain; strict::Bool=true) -> Vector{Tuple{AbstractAtom,AbstractAtom}}

The covalently bonded atom pairs of `chain`, from the known bonds of its
residues (`residuedata`).

Three sources are combined: each residue's internal bonds, the peptide bond
joining residues consecutive in residue number, and the terminal bonds `C`-`OXT`
and `N`-`H1`/`H2`/`H3`, which the standard residue entries do not record. The
terminal bonds follow from the atom names, so terminal residues bond the same
way whether or not [`specializeresnames!`](@ref) has renamed them to `NALA`,
`CALA`, and the like. Hetero residues are skipped. Disordered atoms and residues
are expanded, so every alternative location appears.

By default the chain is required to be `strict`, meaning:

- residue and atom names must be standard
- all hydrogens must be present
- `HIS` must be disambiguated as `HIE`, `HID`, or `HIP`

A missing atom then throws a `KeyError`. With `strict = false` unrecognized
residues are skipped and missing atoms leave their bonds out, at some risk to
accuracy: absent bonds are silent, so a structure that is quietly under-bonded
looks the same as one that is correct.

Bonding is a chemical property, not a geometric one, so the result is
unchanged by conformation.
"""
function atombonds(chain::Chain; strict::Bool=true)
    bonds = Tuple{AbstractAtom,AbstractAtom}[]
    prev = nothing
    for r in chain
        ishetero(r) && continue
        if prev !== nothing && resnumber(r) == resnumber(prev) + 1
            # Peptide bond(s); disordered residues may need more than one
            for _r in collectresidues(r; expand_disordered=true),
                _rp in collectresidues(prev; expand_disordered=true)
                for a in _r["N"], ap in _rp["C"]
                    push!(bonds, (a, ap))
                end
            end
        end
        addresiduebonds!(bonds, r, strict)
        addterminalbonds!(bonds, r, strict)
        prev = r
    end
    return bonds
end

# Atoms that a chain's terminal residues carry but that the standard residue
# entries omit, paired with the heavy atom they bond to. The terminal variants
# (`CALA`, `NALA`, …) do list them, so a residue is only repaired when its own
# entry lacks the bond; otherwise it would be bonded twice.
const terminalbonds = (("C", "OXT"), ("N", "H1"), ("N", "H2"), ("N", "H3"))

function addterminalbonds!(bonds, r::AbstractResidue, strict::Bool)
    rd = get(residuedata, resname(r), nothing)
    for (heavy, terminal) in terminalbonds
        rd !== nothing && any(((a, b),) -> a == terminal || b == terminal, rd.bonds) && continue
        for _r in collectresidues(r; expand_disordered=true)
            aj = findatombyname(_r, terminal; strict=false)
            aj === nothing && continue
            pushbond!(bonds, findatombyname(_r, heavy; strict), aj)
        end
    end
    return bonds
end

function addresiduebonds!(bonds, r::Residue, strict::Bool)
    rname = residuekey(r, strict)
    strict || haskey(residuedata, rname) || return bonds
    # An N-terminal residue named without the `N` prefix carries "H1"/"H2"/"H3"
    # in place of the amide "H", leaving the entry's ("N", "H") bond with no
    # atom to refer to.
    ntermH = findatombyname(r, "H"; strict=false) === nothing &&
             any(n -> findatombyname(r, n; strict=false) !== nothing, ("H1", "H2", "H3"))
    for (ni, nj) in residuedata[rname].bonds
        ntermH && (ni == "H" || nj == "H") && continue
        pushbond!(bonds, findatombyname(r, ni; strict), findatombyname(r, nj; strict))
    end
    return bonds
end

addresiduebonds!(bonds, dr::DisorderedResidue, strict::Bool) = for r in values(dr.names)
    addresiduebonds!(bonds, r, strict)
end

function pushbond!(bonds, ai::Union{AbstractAtom,Nothing}, aj::Union{AbstractAtom,Nothing})
    (ai === nothing || aj === nothing) && return bonds
    for _ai in ai, _aj in aj  # handle disordered atoms
        push!(bonds, (_ai, _aj))
    end
    return bonds
end

function residuekey(r::Residue, strict::Bool)
    rname = resname(r)
    if rname == "HIS"
        strict && error("HIS is not a valid residue name in strict mode; use HIE, HID, or HIP")
        rname = "HIE"
    end
    return rname
end
