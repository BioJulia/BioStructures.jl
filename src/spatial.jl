export
    Transformation,
    coordarray,
    applytransform!,
    applytransform,
    superimpose!,
    rmsd,
    displacements,
    sqdistance,
    distance,
    bondangle,
    dihedralangle,
    omegaangle,
    phiangle,
    psiangle,
    omegaangles,
    phiangles,
    psiangles,
    chiangle,
    chiangles,
    ramachandranangles,
    SpatialMap,
    ContactMap,
    DistanceMap,
    showcontactmap

"""
    Transformation(el1, el2, residue_selectors...)
    Transformation(coords1, coords2)
    Transformation(trans1, trans2, rot)

A 3D transformation to map one set of coordinates onto another, found using the
Kabsch algorithm.

When called with structural elements, carries out a pairwise alignment and
superimposes on atoms from aligned residues.
In this case the BioSequences.jl and BioAlignments.jl packages should be imported.
Keyword arguments for pairwise alignment can be given, see `pairalign`.
The residue selectors determine which residues to do the pairwise alignment on.
The keyword argument `alignatoms` is an atom selector that selects the atoms to
calculate the superimposition on (default `calphaselector`).
Can also be called with two sets of coordinates of the same size, with the
number of dimensions in the first axis and the number of points in the second
axis.

The returned `Transformation` object consists of the mean coordinates of the
first set, the mean coordinates of the second set, the rotation to map the
first centred set onto the second centred set, and the indices of the aligned
residues in the first and second elements if relevant.
"""
struct Transformation
    trans1::Vector{Float64}
    trans2::Vector{Float64}
    rot::Matrix{Float64}
    inds1::Vector{Int}
    inds2::Vector{Int}
end

function Transformation(trans1::AbstractVector{<:Real},
                        trans2::AbstractVector{<:Real},
                        rot::AbstractMatrix{<:Real})
    return Transformation(trans1, trans2, rot, Int[], Int[])
end

Base.show(io::IO, trans::Transformation) = print("3D transformation with ",
    "translation 1 ", trans.trans1, ", translation 2 ", trans.trans2,
    ", rotation ", trans.rot)

"""
    coordarray(element, atom_selectors...)

Get the atomic coordinates in Å of a `StructuralElementOrList` as a 2D `Array`.

Each column corresponds to one atom, so the size is (3, n_atoms).
Additional arguments are atom selector functions - only atoms that return
`true` from all the functions are retained.
The keyword argument `expand_disordered` (default `false`) determines whether to
return coordinates for all copies of disordered atoms separately.
"""
function coordarray(el::StructuralElementOrList,
                    atom_selectors::Function...;
                    expand_disordered::Bool=false)
    at_list = collectatoms(el, atom_selectors...; expand_disordered=expand_disordered)
    coords_out = zeros(3, length(at_list))
    for j in eachindex(at_list)
        coords_out[1, j] = x(at_list[j])
        coords_out[2, j] = y(at_list[j])
        coords_out[3, j] = z(at_list[j])
    end
    return coords_out
end

# Selector functions ignored
coordarray(coords_in::AbstractArray{<:Real}, atom_selectors::Function...) = coords_in

"""
    applytransform!(el, transformation)

Modify all coordinates in an element according to a transformation.
"""
function applytransform!(el::StructuralElementOrList,
                        transformation::Transformation)
    ats = collectatoms(el; expand_disordered=true)
    cs = coordarray(ats)
    new_coords = applytransform(cs, transformation)
    for (i, at) in enumerate(ats)
        coords!(at, new_coords[:, i])
    end
    return el
end

"""
    applytransform(coords, transformation)

Modify coordinates according to a transformation.
"""
function applytransform(cs::AbstractArray{<:Real, 2},
                        t::Transformation)
    return t.rot * (cs .- t.trans1) .+ t.trans2
end

"""
    superimpose!(el1, el2, residue_selectors...)

Calculate the `Transformation` that maps the first element onto the second,
and modify all coordinates in the first element according to the transformation.

Requires the BioSequences.jl and BioAlignments.jl packages to be imported.
See `Transformation` for keyword arguments.
"""
function superimpose!(el1::StructuralElementOrList,
                        el2::StructuralElementOrList,
                        residue_selectors::Function...;
                        kwargs...)
    transformation = Transformation(el1, el2, residue_selectors...; kwargs...)
    applytransform!(el1, transformation)
    return el1
end

function Transformation(coords1::AbstractMatrix{<:Real},
                        coords2::AbstractMatrix{<:Real},
                        inds1::AbstractVector{<:Integer}=Int[],
                        inds2::AbstractVector{<:Integer}=Int[])
    if size(coords1) != size(coords2)
        throw(ArgumentError("Size of coordinate arrays differ: $(size(coords1)) and $(size(coords2))"))
    end
    trans1 = mean(coords1, dims=2)
    trans2 = mean(coords2, dims=2)
    p = coords1 .- trans1
    q = coords2 .- trans2
    # Find the rotation that maps the coordinates
    cov = p * transpose(q)
    svd_res = svd(cov)
    Ut = transpose(svd_res.U)
    # Check sign of determinant
    d = sign(det(svd_res.V * Ut))
    @view(svd_res.V[:, end]) .*= d
    rot = svd_res.V * Ut
    return Transformation(vec(trans1), vec(trans2), rot, inds1, inds2)
end

"""
    rmsd(element_one, element_two, residue_selectors...)
    rmsd(element_one, element_two, superimpose=false)
    rmsd(coords_one, coords_two)

Get the root-mean-square deviation (RMSD) in Å between two
`StructuralElementOrList`s or two coordinate `Array`s.

If `superimpose` is `true` (the default), the elements are superimposed before
RMSD calculation and the RMSD is calculated on the superimposed residues.
In this case the BioSequences.jl and BioAlignments.jl packages should be imported.
See `Transformation` for keyword arguments.
If `superimpose` is `false` the elements are assumed to be superimposed and must
be of the same length.
The keyword argument `rmsdatoms` is an atom selector that selects the atoms to
calculate RMSD on (default `calphaselector`).
"""
function rmsd(coords_one::AbstractArray{<:Real}, coords_two::AbstractArray{<:Real})
    if size(coords_one) != size(coords_two)
        throw(ArgumentError("Coordinate arrays have size $(size(coords_one)) " *
            "and $(size(coords_two)) but must be the same to calculate RMSD"))
    end
    diff = coords_one - coords_two
    return sqrt.(dot(diff, diff) / size(coords_one, 2))
end

function rmsd(el1::StructuralElementOrList,
            el2::StructuralElementOrList,
            residue_selectors::Function...;
            superimpose::Bool=true,
            rmsdatoms::Function=calphaselector,
            kwargs...)
    if superimpose
        res1 = collectresidues(el1, residue_selectors...)
        res2 = collectresidues(el2, residue_selectors...)
        trans = Transformation(res1, res2; kwargs...)
        return rmsd(applytransform(coordarray(res1[trans.inds1], rmsdatoms), trans),
                    coordarray(res2[trans.inds2], rmsdatoms))
    else
        return rmsd(coordarray(el1, rmsdatoms),
                    coordarray(el2, rmsdatoms))
    end
end

"""
    displacements(element_one, element_two, residue_selectors...)
    displacements(element_one, element_two, superimpose=false)
    displacements(coords_one, coords_two)

Get the displacements in Å between atomic coordinates from two
`StructuralElementOrList`s or two coordinate `Array`s.

If `superimpose` is `true` (the default), the elements are superimposed before
calculation and the displacements are calculated on the superimposed residues.
In this case the BioSequences.jl and BioAlignments.jl packages should be imported.
See `Transformation` for keyword arguments.
If `superimpose` is `false` the elements are assumed to be superimposed and must
be of the same length.
The keyword argument `dispatoms` is an atom selector that selects the atoms to
calculate displacements on (default `calphaselector`).
"""
function displacements(coords_one::AbstractArray{<:Real}, coords_two::AbstractArray{<:Real})
    if size(coords_one) != size(coords_two)
        throw(ArgumentError("Coordinate arrays have size $(size(coords_one)) " *
            "and $(size(coords_two)) but must be the same to calculate displacements"))
    end
    diff = coords_one - coords_two
    return sqrt.(sum(diff .* diff, dims=1))[:]
end

function displacements(el1::StructuralElementOrList,
            el2::StructuralElementOrList,
            residue_selectors::Function...;
            superimpose::Bool=true,
            dispatoms::Function=calphaselector,
            kwargs...)
    if superimpose
        res1 = collectresidues(el1, residue_selectors...)
        res2 = collectresidues(el2, residue_selectors...)
        trans = Transformation(res1, res2; kwargs...)
        return displacements(applytransform(coordarray(res1[trans.inds1], dispatoms), trans),
                    coordarray(res2[trans.inds2], dispatoms))
    else
        return displacements(coordarray(el1, dispatoms),
                            coordarray(el2, dispatoms))
    end
end

"""
    sqdistance(element_one, element_two, atom_selectors...)

Get the minimum square distance in Å between two
`StructuralElementOrList`s.

Additional arguments are atom selector functions - only atoms that return
`true` from the functions are retained.
"""
function sqdistance(el1::StructuralElementOrList,
                    el2::StructuralElementOrList,
                    atom_selectors::Function...)
    coords_one = coordarray(el1, atom_selectors...)
    coords_two = coordarray(el2, atom_selectors...)
    min_sq_dist = Inf
    for i in 1:size(coords_one, 2)
        for j in 1:size(coords_two, 2)
            @inbounds sq_dist = (coords_one[1, i] - coords_two[1, j]) ^ 2 + (coords_one[2, i] - coords_two[2, j]) ^ 2 + (coords_one[3, i] - coords_two[3, j]) ^ 2
            if sq_dist < min_sq_dist
                min_sq_dist = sq_dist
            end
        end
    end
    return min_sq_dist
end

function sqdistance(at_one::AbstractAtom, at_two::AbstractAtom)
    return (x(at_one) - x(at_two)) ^ 2 +
           (y(at_one) - y(at_two)) ^ 2 +
           (z(at_one) - z(at_two)) ^ 2
end

"""
    distance(element_one, element_two, atom_selectors...)

Get the minimum distance in Å between two `StructuralElementOrList`s.

Additional arguments are atom selector functions - only atoms that return
`true` from the functions are retained.
"""
function BioGenerics.distance(el1::StructuralElementOrList,
                      el2::StructuralElementOrList,
                      atom_selectors::Function...)
    return sqrt(sqdistance(el1, el2, atom_selectors...))
end

function BioGenerics.distance(at_one::AbstractAtom, at_two::AbstractAtom)
    return sqrt(sqdistance(at_one, at_two))
end

"""
    bondangle(atom_a, atom_b, atom_c)
    bondangle(vec_ba, vec_bc)

Calculate the bond or pseudo-bond angle in radians between three
`AbstractAtom`s or two vectors.

The angle between B→A and B→C is returned in the range 0 to π.
"""
function bondangle(at_a::AbstractAtom,
                at_b::AbstractAtom,
                at_c::AbstractAtom)
    return bondangle(
        coords(at_a) - coords(at_b),
        coords(at_c) - coords(at_b)
    )
end

function bondangle(vec_a::AbstractVector{<:Real}, vec_b::AbstractVector{<:Real})
    return acos(dot(vec_a, vec_b) / (norm(vec_a) * norm(vec_b)))
end

"""
    dihedralangle(atom_a, atom_b, atom_c, atom_d)
    dihedralangle(vec_ab, vec_bc, vec_cd)

Calculate the dihedral angle in radians defined by four `AbstractAtom`s or
three vectors.

The angle between the planes defined by atoms (A, B, C) and (B, C, D) is
returned in the range -π to π.
"""
function dihedralangle(at_a::AbstractAtom,
            at_b::AbstractAtom,
            at_c::AbstractAtom,
            at_d::AbstractAtom)
    return dihedralangle(
        coords(at_b) - coords(at_a),
        coords(at_c) - coords(at_b),
        coords(at_d) - coords(at_c))
end

function dihedralangle(vec_a::AbstractVector{<:Real},
                    vec_b::AbstractVector{<:Real},
                    vec_c::AbstractVector{<:Real})
    return atan(
        norm(vec_b) * dot(vec_a, cross(vec_b, vec_c)),
        dot(cross(vec_a, vec_b), cross(vec_b, vec_c)))
end

"""
    omegaangle(res, res_previous)
    omegaangle(chain, res_id)

Calculate the omega angle in radians for an `AbstractResidue`.

Arguments can either be a residue and the previous residue or a chain and
the position as a residue ID.
The first residue (or one at the given index) requires the atoms "N" and
"CA" and the previous residue requires the atoms "CA" and "C".
The angle is in the range -π to π.
"""
function omegaangle(res::AbstractResidue, res_prev::AbstractResidue)
    at_names = atomnames(res, strip=false)
    at_names_prev = atomnames(res_prev, strip=false)
    if !("CA" in at_names_prev || " CA " in at_names_prev)
        throw(ArgumentError("Atom with atom name \"CA\" not found in previous residue"))
    elseif !("C" in at_names_prev || " C  " in at_names_prev)
        throw(ArgumentError("Atom with atom name \"C\" not found in previous residue"))
    elseif !("N" in at_names || " N  " in at_names)
        throw(ArgumentError("Atom with atom name \"N\" not found in residue"))
    elseif !("CA" in at_names || " CA " in at_names)
        throw(ArgumentError("Atom with atom name \"CA\" not found in residue"))
    end
    return dihedralangle(res_prev["CA"], res_prev["C"], res["N"], res["CA"])
end

function omegaangle(chain::Chain, res_id::Union{Integer, AbstractString})
    inds = findall(r -> r == string(res_id), resids(chain))
    if length(inds) != 1
        throw(ArgumentError("\"$res_id\" is an invalid residue ID"))
    end
    i = inds[1]
    if i == 1 || !sequentialresidues(chain[resids(chain)[i - 1]], chain[resids(chain)[i]])
        throw(ArgumentError("Cannot calculate omega angle for residue \"$res_id\" due to a lack of connected residues"))
    end
    return omegaangle(chain[resids(chain)[i]], chain[resids(chain)[i - 1]])
end

"""
    phiangle(res, res_previous)
    phiangle(chain, res_id)

Calculate the phi angle in radians for an `AbstractResidue`.

Arguments can either be a residue and the previous residue or a chain and the
position as a residue ID.
The first residue (or one at the given index) requires the atoms "N", "CA" and
"C" and the previous residue requires the atom "C".
The angle is in the range -π to π.
"""
function phiangle(res::AbstractResidue, res_prev::AbstractResidue)
    at_names = atomnames(res, strip=false)
    at_names_prev = atomnames(res_prev, strip=false)
    if !("C" in at_names_prev || " C  " in at_names_prev)
        throw(ArgumentError("Atom with atom name \"C\" not found in previous residue"))
    elseif !("N" in at_names || " N  " in at_names)
        throw(ArgumentError("Atom with atom name \"N\" not found in residue"))
    elseif !("CA" in at_names || " CA " in at_names)
        throw(ArgumentError("Atom with atom name \"CA\" not found in residue"))
    elseif !("C" in at_names || " C  " in at_names)
        throw(ArgumentError("Atom with atom name \"C\" not found in residue"))
    end
    return dihedralangle(res_prev["C"], res["N"], res["CA"], res["C"])
end

function phiangle(chain::Chain, res_id::Union{Integer, AbstractString})
    inds = findall(r -> r == string(res_id), resids(chain))
    if length(inds) != 1
        throw(ArgumentError("\"$res_id\" is an invalid residue ID"))
    end
    i = inds[1]
    if i == 1 || !sequentialresidues(chain[resids(chain)[i - 1]], chain[resids(chain)[i]])
        throw(ArgumentError("Cannot calculate phi angle for residue \"$res_id\" due to a lack of connected residues"))
    end
    return phiangle(chain[resids(chain)[i]], chain[resids(chain)[i - 1]])
end

"""
    psiangle(res, res_next)
    psiangle(chain, res_id)

Calculate the psi angle in radians for an `AbstractResidue`.

Arguments can either be a residue and the next residue or a chain and the
position as a residue ID.
The first residue (or one at the given index) requires the atoms "N", "CA" and
"C" and the next residue requires the atom "N".
The angle is in the range -π to π.
"""
function psiangle(res::AbstractResidue, res_next::AbstractResidue)
    at_names = atomnames(res, strip=false)
    at_names_next = atomnames(res_next, strip=false)
    if !("N" in at_names || " N  " in at_names)
        throw(ArgumentError("Atom with atom name \"N\" not found in residue"))
    elseif !("CA" in at_names || " CA " in at_names)
        throw(ArgumentError("Atom with atom name \"CA\" not found in residue"))
    elseif !("C" in at_names || " C  " in at_names)
        throw(ArgumentError("Atom with atom name \"C\" not found in residue"))
    elseif !("N" in at_names_next || " N  " in at_names_next)
        throw(ArgumentError("Atom with atom name \"N\" not found in next residue"))
    end
    return dihedralangle(res["N"], res["CA"], res["C"], res_next["N"])
end

function psiangle(chain::Chain, res_id::Union{Integer, AbstractString})
    inds = findall(r -> r == string(res_id), resids(chain))
    if length(inds) != 1
        throw(ArgumentError("\"$res_id\" is an invalid residue ID"))
    end
    i = inds[1]
    if i == length(chain) || !sequentialresidues(chain[resids(chain)[i]], chain[resids(chain)[i + 1]])
        throw(ArgumentError("Cannot calculate psi angle for residue \"$res_id\" due to a lack of connected residues"))
    end
    return psiangle(chain[resids(chain)[i]], chain[resids(chain)[i + 1]])
end

# source: http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
const chitables = [
    Dict{String,NTuple{4,String}}(
        "ARG" => ("N", "CA", "CB", "CG"),
        "ASN" => ("N", "CA", "CB", "CG"),
        "ASP" => ("N", "CA", "CB", "CG"),
        "CYS" => ("N", "CA", "CB", "SG"),
        "GLN" => ("N", "CA", "CB", "CG"),
        "GLU" => ("N", "CA", "CB", "CG"),
        "HIS" => ("N", "CA", "CB", "CG"),
        "ILE" => ("N", "CA", "CB", "CG1"),
        "LEU" => ("N", "CA", "CB", "CG"),
        "LYS" => ("N", "CA", "CB", "CG"),
        "MET" => ("N", "CA", "CB", "CG"),
        "PHE" => ("N", "CA", "CB", "CG"),
        "PRO" => ("N", "CA", "CB", "CG"),
        "SER" => ("N", "CA", "CB", "OG"),
        "THR" => ("N", "CA", "CB", "OG1"),
        "TRP" => ("N", "CA", "CB", "CG"),
        "TYR" => ("N", "CA", "CB", "CG"),
        "VAL" => ("N", "CA", "CB", "CG1"),
    ),
    Dict{String,NTuple{4,String}}(
        "ARG" => ("CA", "CB", "CG", "CD"),
        "ASN" => ("CA", "CB", "CG", "OD1"),
        "ASP" => ("CA", "CB", "CG", "OD1"),
        "GLN" => ("CA", "CB", "CG", "CD"),
        "GLU" => ("CA", "CB", "CG", "CD"),
        "HIS" => ("CA", "CB", "CG", "ND1"),
        "ILE" => ("CA", "CB", "CG1", "CD1"),
        "LEU" => ("CA", "CB", "CG", "CD1"),
        "LYS" => ("CA", "CB", "CG", "CD"),
        "MET" => ("CA", "CB", "CG", "SD"),
        "PHE" => ("CA", "CB", "CG", "CD1"),
        "PRO" => ("CA", "CB", "CG", "CD"),
        "TRP" => ("CA", "CB", "CG", "CD1"),
        "TYR" => ("CA", "CB", "CG", "CD1"),
    ),
    Dict{String,NTuple{4,String}}(
        "ARG" => ("CB", "CG", "CD", "NE"),
        "GLN" => ("CB", "CG", "CD", "OE1"),
        "GLU" => ("CB", "CG", "CD", "OE1"),
        "LYS" => ("CB", "CG", "CD", "CE"),
        "MET" => ("CB", "CG", "SD", "CE"),
    ),
    Dict{String,NTuple{4,String}}(
        "ARG" => ("CG", "CD", "NE", "CZ"),
        "LYS" => ("CG", "CD", "CE", "NZ"),
    ),
    Dict{String,NTuple{4,String}}(
        "ARG" => ("CD", "NE", "CZ", "NH1"),
    ),
]

# source: Jumper et al 2021, Supp. Table 2
const chisymmetries = [
    ("ASP", 2),
    ("GLU", 3),
    ("PHE", 2),
    ("TYR", 2),
]

function chiangle(a1, a2, a3, a4, sym::Bool)
    χ = dihedralangle(a1, a2, a3, a4)
    if sym && χ < 0
        χ += π
    end
    return χ
end

"""
    chiangle(res, i)

Calculate the χᵢ angle in radians for a standard `AbstractResidue` with standard atom names.
The angle is in the range -π to π.
"""
function chiangle(res::AbstractResidue, i::Integer)
    1 <= i <= 5 || throw(ArgumentError("χᵢ index `i` must be between 1 and 5"))
    ct = chitables[i]
    rname = resname(res)
    t = get(ct, rname, nothing)
    if t === nothing
        # is it because this isn't a standard residue?
        (rname == "GLY" || rname == "ALA") && throw(ArgumentError("$rname does not have any χ angles"))
        haskey(chitables[1], rname) || throw(ArgumentError("no χ angles are defined for residues with name $rname"))
        # it must be missing specifically for `i`
        j = 1
        while haskey(chitables[j+1], rname)
            j += 1
        end
        throw(ArgumentError("χ angle with index $i does not exist for residue $rname (max is $j)"))
    end
    return chiangle(res[t[1]], res[t[2]], res[t[3]], res[t[4]], (rname, i) in chisymmetries)
end

"""
    chiangles(res)

Calculate the `Vector` of standard χ angles for a standard `AbstractResidue` with standard atom names.
The length of the vector ranges from 0 (GLY, ALA) to 5 (ARG).
The angles are each in the range -π to π.
"""
function chiangles(res::AbstractResidue)
    chi = Float64[]
    rname = resname(res)
    for i = 1:5
        ct = chitables[i]
        t = get(ct, rname, nothing)
        if t === nothing
            if i == 1
                (rname == "GLY" || rname == "ALA") && break
                throw(ArgumentError("no χ angles are defined for residues with name $rname"))
            end
            break
        end
        push!(chi, chiangle(res[t[1]], res[t[2]], res[t[3]], res[t[4]], (rname, i) in chisymmetries))
    end
    return chi
end

"""
    omegaangles(element, residue_selectors...)

Calculate the `Vector` of omega angles of a `StructuralElementOrList`.

The vectors have `NaN` for residues where an angle cannot be calculated, e.g.
due to missing atoms or lack of an adjacent residue.
The angle is in the range -π to π.
Additional arguments are residue selector functions - only residues that return
`true` from the functions are retained.
"""
function omegaangles(el::StructuralElementOrList,
                    residue_selectors::Function...)
    res_list = collectresidues(el, residue_selectors...)
    if length(res_list) < 2
        throw(ArgumentError("At least 2 residues required to calculate dihedral angles"))
    end

    omega_angles = Float64[NaN]

    for i in 2:length(res_list)
        res = res_list[i]
        res_prev = res_list[i - 1]
        if sequentialresidues(res_prev, res)
            try
                omega_angle = omegaangle(res, res_prev)
                push!(omega_angles, omega_angle)
            catch ex
                isa(ex, ArgumentError) || rethrow()
                push!(omega_angles, NaN)
            end
        else
            push!(omega_angles, NaN)
        end
    end

    return omega_angles
end

"""
    phiangles(element, residue_selectors...)

Calculate the `Vector` of phi angles of a `StructuralElementOrList`.

The vectors have `NaN` for residues where an angle cannot be calculated, e.g.
due to missing atoms or lack of an adjacent residue.
The angle is in the range -π to π.
Additional arguments are residue selector functions - only residues that return
`true` from the functions are retained.
"""
function phiangles(el::StructuralElementOrList,
                    residue_selectors::Function...)
    res_list = collectresidues(el, residue_selectors...)
    if length(res_list) < 2
        throw(ArgumentError("At least 2 residues required to calculate dihedral angles"))
    end

    phi_angles = Float64[NaN]

    for i in 2:length(res_list)
        res = res_list[i]
        res_prev = res_list[i - 1]
        if sequentialresidues(res_prev, res)
            try
                phi_angle = phiangle(res, res_prev)
                push!(phi_angles, phi_angle)
            catch ex
                isa(ex, ArgumentError) || rethrow()
                push!(phi_angles, NaN)
            end
        else
            push!(phi_angles, NaN)
        end
    end

    return phi_angles
end

"""
    psiangles(element, residue_selectors...)

Calculate the `Vector` of psi angles of a `StructuralElementOrList`.

The vectors have `NaN` for residues where an angle cannot be calculated, e.g.
due to missing atoms or lack of an adjacent residue.
The angle is in the range -π to π.
Additional arguments are residue selector functions - only residues that return
`true` from the functions are retained.
"""
function psiangles(el::StructuralElementOrList,
                    residue_selectors::Function...)
    res_list = collectresidues(el, residue_selectors...)

    if length(res_list) < 2
        throw(ArgumentError("At least 2 residues required to calculate dihedral angles"))
    end

    psi_angles = Float64[]

    for i in 1:(length(res_list) - 1)
        res = res_list[i]
        res_next = res_list[i + 1]
        if sequentialresidues(res, res_next)
            try
                psi_angle = psiangle(res, res_next)
                push!(psi_angles, psi_angle)
            catch ex
                isa(ex, ArgumentError) || rethrow()
                push!(psi_angles, NaN)
            end
        else
            push!(psi_angles, NaN)
        end
    end

    push!(psi_angles, NaN)

    return psi_angles
end

"""
    ramachandranangles(element, residue_selectors...)

Calculate the `Vector`s of phi and psi angles of a `StructuralElementOrList`.

The vectors have `NaN` for residues where an angle cannot be calculated, e.g.
due to missing atoms or lack of an adjacent residue.
The angles are in the range -π to π.
Additional arguments are residue selector functions - only residues that return
`true` from the functions are retained.
"""
function ramachandranangles(el::StructuralElementOrList,
                    residue_selectors::Function...)
    return phiangles(el, residue_selectors...), psiangles(el, residue_selectors...)
end

"""
    chiangles(element, residue_selectors...)

Calculate the `Vector` of standard χ angles for each residue in a `StructuralElementOrList`.
This returns a `Vector` of `Vector`s, where all angles are in the range -π to π.
Additional arguments are residue selector functions - only residues that return
`true` from the functions are retained.
"""
chiangles(el::StructuralElementOrList, residue_selectors::Function...) =
    [chiangles(r) for r in collectresidues(el, residue_selectors...)]

"A map of a structural property, e.g. a `ContactMap` or a `DistanceMap`."
abstract type SpatialMap end

"""
    ContactMap(element, contact_distance)
    ContactMap(element_one, element_two, contact_distance)
    ContactMap(bit_array_2D)

Calculate the contact map for a `StructuralElementOrList`, or between two
`StructuralElementOrList`s.

This returns a `ContactMap` type containing a `BitArray{2}` with `true` where
the sub-elements are no further than the contact distance and `false` otherwise.
When one element is given as input this returns a symmetric square matrix.
To directly access the underlying data of `ContactMap` `cm`, use `cm.data`.

# Examples
```julia
cbetas_A = collectatoms(struc["A"], cbetaselector)
cbetas_B = collectatoms(struc["B"], cbetaselector)

# Contact map of chain A using conventional Cβ and 8 Å definitions
cm = ContactMap(cbetas_A, 8.0)

# Returns true if a contact is present between the tenth and twentieth element
cm[10, 20]

# Rectangular contact map of chains A and B
cm = ContactMap(cbetas_A, cbetas_B, 8.0)

# Write the contact map to file
using DelimitedFiles
writedlm("contacts.out", Int64.(cm.data), " ")
```
"""
struct ContactMap <: SpatialMap
    data::BitArray{2}
end

"""
    DistanceMap(element)
    DistanceMap(element_one, element_two)
    DistanceMap(float_array_2D)

Calculate the distance map for a `StructuralElementOrList`, or between two
`StructuralElementOrList`s.

This returns a `DistanceMap` type containing a `Array{Float64, 2}` with minimum
distances between the sub-elements.
When one element is given as input this returns a symmetric square matrix.
To directly access the underlying data of `DistanceMap` `dm`, use `dm.data`.

# Examples
```julia
cbetas_A = collectatoms(struc["A"], cbetaselector)
cbetas_B = collectatoms(struc["B"], cbetaselector)

# Distance map of chain A showing how far each Cβ atom is from the others
dm = DistanceMap(cbetas_A)

# Returns the distance between the tenth and twentieth element
dm[10, 20]

# Rectangular distance map of chains A and B
dm = DistanceMap(cbetas_A, cbetas_B)

# Write the distance map to file
using DelimitedFiles
writedlm("distances.out", dm.data, " ")
```
"""
struct DistanceMap <: SpatialMap
    data::Array{Float64, 2}
end

Base.size(m::SpatialMap) = size(m.data)
Base.size(m::SpatialMap, dim::Integer) = size(m.data, dim)

Base.getindex(m::SpatialMap, args::Integer...) = m.data[args...]

function Base.setindex!(m::SpatialMap, v, args::Integer...)
    m.data[args...] = v
    return m
end

Base.firstindex(::SpatialMap) = 1
Base.firstindex(::SpatialMap, dim) = 1
Base.lastindex(m::SpatialMap) = reduce(*, size(m))
Base.lastindex(m::SpatialMap, dim) = size(m, dim)

function Base.show(io::IO, cm::ContactMap)
    print(io, "Contact map of size $(size(cm))")
end

function Base.show(io::IO, dm::DistanceMap)
    print(io, "Distance map of size $(size(dm))")
end

function ContactMap(el1::StructuralElementOrList,
                el2::StructuralElementOrList,
                contact_dist::Real)
    sq_contact_dist = contact_dist ^ 2
    contacts = falses(length(el1), length(el2))
    for (i, subel1) in enumerate(el1)
        for (j, subel2) in enumerate(el2)
            if sqdistance(subel1, subel2) <= sq_contact_dist
                contacts[i, j] = true
            end
        end
    end
    return ContactMap(contacts)
end

function ContactMap(el::StructuralElementOrList, contact_dist::Real)
    sq_contact_dist = contact_dist ^ 2
    contacts = falses(length(el), length(el))
    el_list = collect(el)
    for i in 1:length(el)
        contacts[i, i] = true
        for j in 1:(i - 1)
            if sqdistance(el_list[i], el_list[j]) <= sq_contact_dist
                contacts[i, j] = true
                contacts[j, i] = true
            end
        end
    end
    return ContactMap(contacts)
end

function DistanceMap(el1::StructuralElementOrList,
                el2::StructuralElementOrList)
    dists = zeros(length(el1), length(el2))
    for (i, subel1) in enumerate(el1)
        for (j, subel2) in enumerate(el2)
            dists[i, j] = distance(subel1, subel2)
        end
    end
    return DistanceMap(dists)
end

function DistanceMap(el::StructuralElementOrList)
    dists = zeros(length(el), length(el))
    el_list = collect(el)
    for i in 1:length(el)
        for j in 1:(i - 1)
            dist = distance(el_list[i], el_list[j])
            dists[i, j] = dist
            dists[j, i] = dist
        end
    end
    return DistanceMap(dists)
end

# Plot recipe to show a ContactMap
@recipe function plot(cm::ContactMap)
    seriestype := :heatmap
    fillcolor --> :dense
    aspect_ratio --> 1
    xmirror --> true
    colorbar --> false
    xs = string.(1:size(cm, 2))
    ys = string.(size(cm, 1):-1:1)
    xs, ys, reverse(cm.data, dims=1)
end

# Plot recipe to show a DistanceMap
@recipe function plot(dm::DistanceMap)
    seriestype := :heatmap
    fillcolor --> :inferno
    aspect_ratio --> 1
    xmirror --> true
    colorbar --> true
    xs = string.(1:size(dm, 2))
    ys = string.(size(dm, 1):-1:1)
    xs, ys, reverse(dm.data, dims=1)
end

"""
    showcontactmap(contact_map)
    showcontactmap(io, contact_map)

Print a representation of a `ContactMap` to `stdout`, or a specified `IO`
instance.

A fully plotted version can be obtained with `plot(contact_map)` but that
requires Plots.jl; `showcontactmap` works without that dependency.
"""
function showcontactmap(io::IO, cm::ContactMap)
    size1 = size(cm, 1)
    # Print two y values to each line for a nicer output
    for i in 1:2:size1
        for j in 1:size(cm, 2)
            cont_one = cm[i, j]
            # Check we aren't going over the end of the map
            cont_two = i + 1 <= size1 && cm[i + 1, j]
            if cont_one && cont_two
                char = "█"
            elseif cont_one
                char = "▀"
            elseif cont_two
                char = "▄"
            else
                char = " "
            end
            print(io, char)
        end
        println(io)
    end
end

showcontactmap(cm::ContactMap) = showcontactmap(stdout, cm)
