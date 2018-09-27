export
    coordarray,
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
    ramachandranangles,
    SpatialMap,
    ContactMap,
    DistanceMap,
    showcontactmap


"""
    coordarray(element, atom_selectors...)

Get the atomic coordinates in Å of a `StructuralElementOrList` as a 2D
`Array` with each column corresponding to one atom.

Additional arguments are atom selector functions - only atoms that return
`true` from all the functions are retained.
"""
function coordarray(el::StructuralElementOrList, atom_selectors::Function...)
    at_list = collectatoms(el, atom_selectors...)
    coords_out = zeros(3, length(at_list))
    for j in eachindex(at_list)
        coords_out[1, j] = x(at_list[j])
        coords_out[2, j] = y(at_list[j])
        coords_out[3, j] = z(at_list[j])
    end
    return coords_out
end

# Selector functions ignored
coordarray(coords_in::Array{Float64}, atom_selectors::Function...) = coords_in


"""
    rmsd(element_one, element_two, atom_selectors...)
    rmsd(coords_one, coords_two)

Get the root-mean-square deviation (RMSD) in Å between two
`StructuralElementOrList`s or two coordinate `Array`s of the same size.

Assumes they are already superimposed.
Additional arguments are atom selector functions - only atoms that return
`true` from the functions are retained.
"""
function rmsd(coords_one::Array{Float64}, coords_two::Array{Float64})
    if size(coords_one) != size(coords_two)
        throw(ArgumentError("Sizes of coordinate arrays are different - cannot calculate RMSD"))
    end
    diff = coords_one - coords_two
    return sqrt.(dot(diff, diff) / size(coords_one, 2))
end

function rmsd(el_one::StructuralElementOrList,
            el_two::StructuralElementOrList,
            atom_selectors::Function...)
    return rmsd(
        coordarray(el_one, atom_selectors...),
        coordarray(el_two, atom_selectors...))
end


"""
    displacements(element_one, element_two, atom_selectors...)
    displacements(coords_one, coords_two)

Get the displacements in Å between atomic coordinates from two
`StructuralElementOrList`s or two coordinate `Array`s of the same size.

Assumes they are already superimposed.
Additional arguments are atom selector functions - only atoms that return
`true` from the functions are retained.
"""
function displacements(coords_one::Array{Float64}, coords_two::Array{Float64})
    if size(coords_one) != size(coords_two)
        throw(ArgumentError("Sizes of coordinate arrays are different - cannot calculate displacements"))
    end
    diff = coords_one - coords_two
    return sqrt.(sum(diff .* diff, dims=1))[:]
end

function displacements(el_one::StructuralElementOrList,
                    el_two::StructuralElementOrList,
                    atom_selectors::Function...)
    return displacements(
        coordarray(el_one, atom_selectors...),
        coordarray(el_two, atom_selectors...))
end


"""
    sqdistance(element_one, element_two, atom_selectors...)

Get the minimum square distance in Å between two
`StructuralElementOrList`s.

Additional arguments are atom selector functions - only atoms that return
`true` from the functions are retained.
"""
function sqdistance(el_one::StructuralElementOrList,
                    el_two::StructuralElementOrList,
                    atom_selectors::Function...)
    coords_one = coordarray(el_one, atom_selectors...)
    coords_two = coordarray(el_two, atom_selectors...)
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
function distance(el_one::StructuralElementOrList,
                      el_two::StructuralElementOrList,
                      atom_selectors::Function...)
    return sqrt(sqdistance(el_one, el_two, atom_selectors...))
end

function distance(at_one::AbstractAtom, at_two::AbstractAtom)
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

function bondangle(vec_a::Vector{Float64}, vec_b::Vector{Float64})
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

function dihedralangle(vec_a::Vector{Float64},
                    vec_b::Vector{Float64},
                    vec_c::Vector{Float64})
    return atan(
        dot(cross(cross(vec_a, vec_b), cross(vec_b, vec_c)), vec_b / norm(vec_b)),
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
    if !("CA" in atomnames(res_prev))
        throw(ArgumentError("Atom with atom name \"CA\" not found in previous residue"))
    elseif !("C"  in atomnames(res_prev))
        throw(ArgumentError("Atom with atom name \"C\" not found in previous residue"))
    elseif !("N"  in atomnames(res))
        throw(ArgumentError("Atom with atom name \"N\" not found in residue"))
    elseif !("CA" in atomnames(res))
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
    if i == 1 || !sequentialresidues(chain[resids(chain)[i-1]], chain[resids(chain)[i]])
        throw(ArgumentError("Cannot calculate omega angle for residue \"$res_id\" due to a lack of connected residues"))
    end
    return omegaangle(chain[resids(chain)[i]], chain[resids(chain)[i-1]])
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
    if !("C"  in atomnames(res_prev))
        throw(ArgumentError("Atom with atom name \"C\" not found in previous residue"))
    elseif !("N"  in atomnames(res))
        throw(ArgumentError("Atom with atom name \"N\" not found in residue"))
    elseif !("CA" in atomnames(res))
        throw(ArgumentError("Atom with atom name \"CA\" not found in residue"))
    elseif !("C"  in atomnames(res))
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
    if i == 1 || !sequentialresidues(chain[resids(chain)[i-1]], chain[resids(chain)[i]])
        throw(ArgumentError("Cannot calculate phi angle for residue \"$res_id\" due to a lack of connected residues"))
    end
    return phiangle(chain[resids(chain)[i]], chain[resids(chain)[i-1]])
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
    if !("N"  in atomnames(res))
        throw(ArgumentError("Atom with atom name \"N\" not found in residue"))
    elseif !("CA" in atomnames(res))
        throw(ArgumentError("Atom with atom name \"CA\" not found in residue"))
    elseif !("C"  in atomnames(res))
        throw(ArgumentError("Atom with atom name \"C\" not found in residue"))
    elseif !("N"  in atomnames(res_next))
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
    if i == length(chain) || !sequentialresidues(chain[resids(chain)[i]], chain[resids(chain)[i+1]])
        throw(ArgumentError("Cannot calculate psi angle for residue \"$res_id\" due to a lack of connected residues"))
    end
    return psiangle(chain[resids(chain)[i]], chain[resids(chain)[i+1]])
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
        res_prev = res_list[i-1]
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
        res_prev = res_list[i-1]
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

    for i in 1:length(res_list)-1
        res = res_list[i]
        res_next = res_list[i+1]
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

# Examples
```julia
cbetas_A = collectatoms(struc["A"], cbetaselector)
cbetas_B = collectatoms(struc["B"], cbetaselector)

# Contact map of chain A using standard C-beta and 8.0 Å definitions
ContactMap(cbetas_A, 8.0)

# Rectangular contact map of chains A and B
ContactMap(cbetas_A, cbetas_B, 8.0)
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

# Examples
```julia
cbetas_A = collectatoms(struc["A"], cbetaselector)
cbetas_B = collectatoms(struc["B"], cbetaselector)

# Distance map of chain A showing how far each residue is from the others
DistanceMap(cbetas_A)

# Rectangular distance map of chains A and B
DistanceMap(cbetas_A, cbetas_B)
```
"""
struct DistanceMap <: SpatialMap
    data::Array{Float64, 2}
end


Base.getindex(m::SpatialMap, args::Integer...) = m.data[args...]

function Base.setindex!(m::SpatialMap, v, args::Integer...)
    m.data[args...] = v
    return m
end

Base.size(m::SpatialMap) = size(m.data)
Base.size(m::SpatialMap, dim::Integer) = size(m.data, dim)

function Base.show(io::IO, cm::ContactMap)
    print(io, "Contact map of size $(size(cm))")
end

function Base.show(io::IO, dm::DistanceMap)
    print(io, "Distance map of size $(size(dm))")
end


function ContactMap(el_one::StructuralElementOrList,
                el_two::StructuralElementOrList,
                contact_dist::Real)
    sq_contact_dist = contact_dist ^ 2
    contacts = falses(length(el_one), length(el_two))
    for (i, subel_one) in enumerate(el_one)
        for (j, subel_two) in enumerate(el_two)
            if sqdistance(subel_one, subel_two) <= sq_contact_dist
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
        for j in 1:i-1
            if sqdistance(el_list[i], el_list[j]) <= sq_contact_dist
                contacts[i, j] = true
                contacts[j, i] = true
            end
        end
    end
    return ContactMap(contacts)
end

function DistanceMap(el_one::StructuralElementOrList,
                el_two::StructuralElementOrList)
    dists = zeros(length(el_one), length(el_two))
    for (i, subel_one) in enumerate(el_one)
        for (j, subel_two) in enumerate(el_two)
            dists[i, j] = distance(subel_one, subel_two)
        end
    end
    return DistanceMap(dists)
end

function DistanceMap(el::StructuralElementOrList)
    dists = zeros(length(el), length(el))
    el_list = collect(el)
    for i in 1:length(el)
        for j in 1:i-1
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
    aspectratio --> 1
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
    aspectratio --> 1
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
