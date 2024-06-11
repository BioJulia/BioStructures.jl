export
    StructuralElement,
    AbstractAtom,
    PDBConsistencyError,
    Atom,
    DisorderedAtom,
    AbstractResidue,
    Residue,
    DisorderedResidue,
    Chain,
    Model,
    MolecularStructure,
    AtomRecord,
    StructuralElementOrList,
    serial,
    atomname,
    altlocid,
    coords,
    coords!,
    occupancy,
    tempfactor,
    element,
    charge,
    residue,
    ishetero,
    isdisorderedatom,
    defaultaltlocid,
    defaultatom,
    altlocids,
    atomid,
    resname,
    resnumber,
    inscode,
    resid,
    atomnames,
    atoms,
    isdisorderedres,
    disorderedres,
    defaultresname,
    defaultresidue,
    resnames,
    chain,
    chainid,
    chainid!,
    resids,
    residues,
    model,
    modelnumber,
    chainids,
    chains,
    structure,
    structurename,
    modelnumbers,
    models,
    defaultmodel,
    sequentialresidues,
    collect,
    collectmodels,
    countmodels,
    collectchains,
    countchains,
    collectresidues,
    countresidues,
    collectatoms,
    countatoms,
    choosedefaultaltlocid,
    threeletter_to_aa,
    PDB,
    PDBXML,
    MMCIF,
    MMTF,
    pdbextension,
    generatechainid,
    MMTFDict,
    writemmtf

"A macromolecular structural element."
abstract type StructuralElement end

"""
An atom that is part of a macromolecule - either an `Atom` or a
`DisorderedAtom`.
"""
abstract type AbstractAtom <: StructuralElement end

"""
Error arising from an attempt to make an inconsistent structural state.
"""
struct PDBConsistencyError <: Exception
    message::String
end

Base.showerror(io::IO, e::PDBConsistencyError) = print(io, "PDBConsistencyError: ", e.message)

"An atom that is part of a macromolecule."
struct Atom <: AbstractAtom
    serial::Int
    name::String
    alt_loc_id::Char
    coords::Vector{Float64}
    occupancy::Float64
    temp_factor::Float64
    element::String
    charge::String
    residue::StructuralElement
end

"A container to hold different locations of the same atom."
struct DisorderedAtom <: AbstractAtom
    alt_loc_ids::Dict{Char, Atom}
    default::Char
end

"""
A residue (amino acid) or other molecule - either a `Residue` or a
`DisorderedResidue`.
"""
abstract type AbstractResidue <: StructuralElement end

"A residue (amino acid) or other molecule."
mutable struct Residue <: AbstractResidue
    name::String
    number::Int
    ins_code::Char
    het_res::Bool # Does the residue consist of hetatoms?
    atom_list::Vector{String}
    atoms::Dict{String, AbstractAtom}
    chain::StructuralElement
    ss_code::Char
end

"""
A container to hold different versions of the same residue (point
mutations).
"""
struct DisorderedResidue <: AbstractResidue
    names::Dict{String, Residue}
    default::String
end

"A chain (molecule) from a macromolecular structure."
mutable struct Chain <: StructuralElement
    id::String # mmCIF files can have multi-character chain IDs
    res_list::Vector{String}
    residues::Dict{String, AbstractResidue}
    model::StructuralElement
end

"A conformation of a macromolecular structure."
struct Model <: StructuralElement
    number::Int
    chains::Dict{String, Chain}
    structure::StructuralElement
end

"""
A container for multiple `Model`s that represents a Protein Data Bank (PDB)
entry.
"""
struct MolecularStructure <: StructuralElement
    name::String
    models::Dict{Int, Model}
end

"""
A record for a single atom, e.g. as represented in a Protein Data Bank
(PDB) file.
"""
struct AtomRecord
    het_atom::Bool
    serial::Int
    atom_name::String
    alt_loc_id::Char
    res_name::String
    chain_id::String
    res_number::Int
    ins_code::Char
    coords::Vector{Float64}
    occupancy::Float64
    temp_factor::Float64
    element::String
    charge::String
end

# Constructors without sub-elements and without super-elements

MolecularStructure(name::AbstractString) = MolecularStructure(name, Dict())

MolecularStructure() = MolecularStructure("")

Model(number::Integer, struc::MolecularStructure) = Model(number, Dict(), struc)

Model(number::Integer) = Model(number, MolecularStructure())

Model() = Model(1)

Chain(id::AbstractString, mo::Model) = Chain(id, [], Dict(), mo)
Chain(id::Char, mo::Model) = Chain(string(id), [], Dict(), mo)

Chain(id::Union{AbstractString, Char}) = Chain(id, Model())

function Residue(name::AbstractString,
                number::Integer,
                ins_code::Char,
                het_res::Bool,
                ch::Chain,
                ss_code=ss_code_unassigned)
    return Residue(name, number, ins_code, het_res, [], Dict(), ch, ss_code_unassigned)
end

"""
A `StructuralElement` or `Vector` of `StructuralElement`s up to
a `Vector{Model}`.
"""
const StructuralElementOrList = Union{
        StructuralElement,
        Vector{Model},
        Vector{Chain},
        Vector{<:AbstractResidue},
        Vector{<:AbstractAtom},
    }

# Allow accessing sub elements contained in an element like a dictionary
# e.g. allows you to do res[atom_name] rather than res.atoms[atom_name]
# setindex! should be used with caution as it can lead to inconsistencies
# e.g. adding an atom to a residue atom dict does not update the atom list
# Hence it is used internally but not documented for general use

# Accessing a DisorderedAtom with a character returns the Atom with that alt
#   loc ID
Base.getindex(dis_at::DisorderedAtom, alt_loc_id::Char) = dis_at.alt_loc_ids[alt_loc_id]

function Base.setindex!(dis_at::DisorderedAtom, at::Atom, alt_loc_id::Char)
    dis_at.alt_loc_ids[alt_loc_id] = at
    return dis_at
end

Base.firstindex(dis_at::DisorderedAtom) = first(altlocids(dis_at))
Base.lastindex(dis_at::DisorderedAtom) = last(altlocids(dis_at))

# Accessing a Residue with an AbstractString returns the AbstractAtom with that
#   atom name
Base.getindex(res::Residue, atom_name::AbstractString) = findatombyname(res, atom_name)

function Base.setindex!(res::Residue, at::AbstractAtom, atom_name::AbstractString)
    res.atoms[atom_name] = at
    return res
end

Base.firstindex(res::Residue) = first(atomnames(res, strip=false))
Base.lastindex(res::Residue) = last(atomnames(res, strip=false))

# Accessing a DisorderedResidue with an AbstractString returns the AbstractAtom
#   in the default Residue with that atom name
# This is not necessarily intuitive, it may be expected to return the Residue
#   with that residue name
# However this way accessing an AbstractResidue always returns an AbstractAtom
Base.getindex(dis_res::DisorderedResidue, atom_name::AbstractString) = findatombyname(defaultresidue(dis_res), atom_name)

function Base.setindex!(dis_res::DisorderedResidue, at::AbstractAtom, atom_name::AbstractString)
    defaultresidue(dis_res)[atom_name] = at
    return dis_res
end

Base.firstindex(dis_res::DisorderedResidue) = first(atomnames(defaultresidue(dis_res), strip=false))
Base.lastindex(dis_res::DisorderedResidue) = last(atomnames(defaultresidue(dis_res), strip=false))

# Accessing a Chain with an AbstractString returns the AbstractResidue with that
#   residue ID
Base.getindex(ch::Chain, res_id::AbstractString) = ch.residues[res_id]

function Base.setindex!(ch::Chain, res::AbstractResidue, res_id::AbstractString)
    ch.residues[res_id] = res
    return ch
end

# Accessing a Chain with an Integer returns the AbstractResidue with that residue ID
#   converted to a string
Base.getindex(ch::Chain, res_n::Integer) = ch.residues[string(res_n)]

function Base.setindex!(ch::Chain, res::AbstractResidue, res_n::Integer)
    ch.residues[string(res_n)] = res
    return ch
end

Base.firstindex(ch::Chain) = first(resids(ch))
Base.lastindex(ch::Chain) = last(resids(ch))

# Accessing a Model with a Char or AbstractString returns the Chain with that
#   chain ID
Base.getindex(mo::Model, ch_id::AbstractString) = mo.chains[ch_id]
Base.getindex(mo::Model, ch_id::Char) = mo.chains[string(ch_id)]

function Base.setindex!(mo::Model, ch::Chain, ch_id::AbstractString)
    mo.chains[ch_id] = ch
    return mo
end

function Base.setindex!(mo::Model, ch::Chain, ch_id::Char)
    return setindex!(mo, ch, string(ch_id))
end

Base.firstindex(mo::Model) = first(chainids(mo))
Base.lastindex(mo::Model) = last(chainids(mo))

# Accessing a MolecularStructure with an Integer returns the Model with that model
#   number
Base.getindex(struc::MolecularStructure, mo_n::Integer) = struc.models[mo_n]

function Base.setindex!(struc::MolecularStructure, mo::Model, mo_n::Integer)
    struc.models[mo_n] = mo
    return struc
end

# Accessing a MolecularStructure with a Char returns the Chain with that chain ID
#   on the default model
Base.getindex(struc::MolecularStructure, ch_id::AbstractString) = defaultmodel(struc)[ch_id]
Base.getindex(struc::MolecularStructure, ch_id::Char) = defaultmodel(struc)[string(ch_id)]

function Base.setindex!(struc::MolecularStructure, ch::Chain, ch_id::AbstractString)
    defaultmodel(struc)[ch_id] = ch
    return struc
end

function Base.setindex!(struc::MolecularStructure, ch::Chain, ch_id::Char)
    return setindex!(struc, ch, string(ch_id))
end

Base.firstindex(struc::MolecularStructure) = first(modelnumbers(struc))
Base.lastindex(struc::MolecularStructure) = last(modelnumbers(struc))

# Check if an atom name exists in a residue as a whitespace-padded version
function findatombyname(res::Residue, atom_name::AbstractString)
    # Look for atom name directly
    if haskey(res.atoms, atom_name)
        return res.atoms[atom_name]
    # Pad out name to 4 characters to read PDB atom names with whitespace
    elseif length(atom_name) == 3
        if haskey(res.atoms, " $atom_name")
            return res.atoms[" $atom_name"]
        elseif haskey(res.atoms, "$atom_name ")
            return res.atoms["$atom_name "]
        end
    elseif length(atom_name) == 2
        if haskey(res.atoms, " $atom_name ")
            return res.atoms[" $atom_name "]
        elseif haskey(res.atoms, "  $atom_name")
            return res.atoms["  $atom_name"]
        elseif haskey(res.atoms, "$atom_name  ")
            return res.atoms["$atom_name  "]
        end
    elseif length(atom_name) == 1
        if haskey(res.atoms, " $atom_name  ")
            return res.atoms[" $atom_name  "]
        elseif haskey(res.atoms, "  $atom_name ")
            return res.atoms["  $atom_name "]
        elseif haskey(res.atoms, "   $atom_name")
            return res.atoms["   $atom_name"]
        elseif haskey(res.atoms, "$atom_name   ")
            return res.atoms["$atom_name   "]
        end
    end
    # Could not find atom name
    throw(KeyError(atom_name))
end

# Getters and setters for structural elements

"""
    serial(at)

Get the serial number of an `AbstractAtom` as an `Int`.
"""
serial(at::Atom) = at.serial
serial(dis_at::DisorderedAtom) = serial(defaultatom(dis_at))

"""
    atomname(at; strip=true)

Get the atom name of an `AbstractAtom` as a `String`.
`strip` determines whether surrounding whitespace is stripped.
"""
function atomname(at::Atom; strip::Bool=true)
    if strip
        return Base.strip(at.name)
    else
        return at.name
    end
end

atomname(dis_at::DisorderedAtom; strip::Bool=true) = atomname(defaultatom(dis_at), strip=strip)

"""
    altlocid(at)

Get the alternative location ID of an `AbstractAtom` as a `Char`.
"""
altlocid(at::Atom) = at.alt_loc_id
altlocid(dis_at::DisorderedAtom) = defaultaltlocid(dis_at)

"""
    x(at)

Get the x coordinate in Å of an `AbstractAtom` as a `Float64`.
"""
x(at::Atom) = at.coords[1]
x(dis_at::DisorderedAtom) = x(defaultatom(dis_at))

"""
    x!(at, val)

Set the x coordinate in Å of an `AbstractAtom` to `val`.

For `DisorderedAtom`s only the default atom is updated.
"""
x!(at::Atom, x::Real) = (at.coords[1] = x; at)
x!(dis_at::DisorderedAtom, x::Real) = x!(defaultatom(dis_at), x)

"""
    y(at)

Get the y coordinate in Å of an `AbstractAtom` as a `Float64`.
"""
y(at::Atom) = at.coords[2]
y(dis_at::DisorderedAtom) = y(defaultatom(dis_at))

"""
    y!(at, val)

Set the y coordinate in Å of an `AbstractAtom` to `val`.

For `DisorderedAtom`s only the default atom is updated.
"""
y!(at::Atom, y::Real) = (at.coords[2] = y; at)
y!(dis_at::DisorderedAtom, y::Real) = y!(defaultatom(dis_at), y)

"""
    z(at)

Get the z coordinate in Å of an `AbstractAtom` as a `Float64`.
"""
z(at::Atom) = at.coords[3]
z(dis_at::DisorderedAtom) = z(defaultatom(dis_at))

"""
    z!(at, val)

Set the z coordinate in Å of an `AbstractAtom` to `val`.

For `DisorderedAtom`s only the default atom is updated.
"""
z!(at::Atom, z::Real) = (at.coords[3] = z; at)
z!(dis_at::DisorderedAtom, z::Real) = z!(defaultatom(dis_at), z)

"""
    coords(at)

Get the coordinates in Å of an `AbstractAtom` as a `Vector{Float64}`.
"""
coords(at::Atom) = at.coords
coords(dis_at::DisorderedAtom) = coords(defaultatom(dis_at))

"""
    coords!(at, new_coords)

Set the coordinates in Å of an `AbstractAtom` to a `Vector` of 3 numbers.

For `DisorderedAtom`s only the default atom is updated.
"""
function coords!(at::Atom, new_coords)
    if length(new_coords) != 3
        throw(ArgumentError("3 coordinates must be given"))
    end
    x!(at, new_coords[1])
    y!(at, new_coords[2])
    z!(at, new_coords[3])
    return at
end

function coords!(dis_at::DisorderedAtom, new_coords)
    coords!(defaultatom(dis_at), new_coords)
end

"""
    occupancy(at)

Get the occupancy of an `AbstractAtom` as a `Float64`.

The occupancy is set to `1.0` if not specified during atom creation.
"""
occupancy(at::Atom) = at.occupancy
occupancy(dis_at::DisorderedAtom) = occupancy(defaultatom(dis_at))

"""
    tempfactor(at)

Get the temperature factor of an `AbstractAtom` as a `Float64`.

The temperature factor is set to `0.0` if not specified during atom creation.
"""
tempfactor(at::Atom) = at.temp_factor
tempfactor(dis_at::DisorderedAtom) = tempfactor(defaultatom(dis_at))

"""
    element(at; strip=true)

Get the element of an `AbstractAtom` as a `String`.

The element is set to `"  "` if not specified during atom creation.
`strip` determines whether surrounding whitespace is stripped.
"""
function element(at::Atom; strip::Bool=true)
    if strip
        return Base.strip(at.element)
    else
        return at.element
    end
end

element(dis_at::DisorderedAtom; strip::Bool=true) = element(defaultatom(dis_at), strip=strip)

"""
    charge(at; strip=true)

Get the charge on an `AbstractAtom` as a `String`.
The charge is set to `"  "` if not specified during atom creation.
`strip` determines whether surrounding whitespace is stripped.
"""
function charge(at::Atom; strip::Bool=true)
    if strip
        return Base.strip(at.charge)
    else
        return at.charge
    end
end

charge(dis_at::DisorderedAtom; strip::Bool=true) = charge(defaultatom(dis_at), strip=strip)

"""
    residue(at)

Get the `Residue` that an `AbstractAtom` belongs to.
"""
residue(at::Atom) = at.residue
residue(dis_at::DisorderedAtom) = residue(defaultatom(dis_at))
residue(res::AbstractResidue) = res

"""
    ishetero(at)
    ishetero(res)

Determines if an `AbstractAtom` represents a hetero atom, e.g. came from a
HETATM record in a Protein Data Bank (PDB) file, or if an `AbstractResidue`
represents a hetero molecule, e.g. consists of HETATM records from a PDB file.
"""
ishetero(at::AbstractAtom) = ishetero(residue(at))
ishetero(res::Residue) = res.het_res
ishetero(dis_res::DisorderedResidue) = ishetero(defaultresidue(dis_res))

"""
    isdisorderedatom(at)

Determines if an `AbstractAtom` is a `DisorderedAtom`, i.e. if there are
multiple locations present for an atom.
"""
isdisorderedatom(::Atom) = false
isdisorderedatom(::DisorderedAtom) = true

"""
    defaultaltlocid(dis_at)

Get the alternative location ID of the default `Atom` in a `DisorderedAtom` as a
`Char`.

The default is the highest occupancy, or lowest character alternative location
ID for ties (i.e. 'A' beats 'B').
"""
defaultaltlocid(dis_at::DisorderedAtom) = dis_at.default

# Constructor acts as a setter for the default alt loc ID
function DisorderedAtom(dis_at::DisorderedAtom, default::Char)
    if !(default in altlocids(dis_at))
        throw(ArgumentError("The new default alternative location ID \"$default\" must be present in the atom"))
    end
    return DisorderedAtom(dis_at.alt_loc_ids, default)
end

"""
    defaultatom(dis_at)

Return the default `Atom` in a `DisorderedAtom`.

The default is the highest occupancy, or lowest character alternative location
ID for ties (i.e. 'A' beats 'B').
"""
defaultatom(dis_at::DisorderedAtom) = dis_at[defaultaltlocid(dis_at)]

"""
    altlocids(dis_at)

Get the list of alternative location IDs in an `AbstractAtom` as a
`Vector{Char}`, sorted by atom serial.
"""
function altlocids(dis_at::DisorderedAtom)
    return sort(collect(keys(dis_at.alt_loc_ids)),
        by=alt_loc_id -> dis_at[alt_loc_id])
end

altlocids(at::Atom) = [altlocid(at)]

"""
    atomid(at)

Get a descriptive atom ID for an `AbstractAtom` as a `Tuple` of the form
(full residue ID, residue name, atom name).
"""
atomid(at::Atom) = (resid(at, full=true), resname(at), atomname(at))
atomid(dis_at::DisorderedAtom) = atomid(defaultatom(dis_at))

"""
    resname(at; strip=true)
    resname(res; strip=true)

Get the residue name of an `AbstractAtom` or `AbstractResidue` as a `String`.

`strip` determines whether surrounding whitespace is stripped.
"""
function resname(res::Residue; strip::Bool=true)
    if strip
        return Base.strip(res.name)
    else
        return res.name
    end
end

resname(at::Atom; strip::Bool=true) = resname(residue(at), strip=strip)
resname(dis_at::DisorderedAtom; strip::Bool=true) = resname(defaultatom(dis_at), strip=strip)
resname(dis_res::DisorderedResidue; strip::Bool=true) = resname(defaultresidue(dis_res), strip=strip)

"""
    resnumber(at)
    resnumber(res)

Get the residue number of an `AbstractAtom` or `AbstractResidue` as an `Int`.
"""
resnumber(at::Atom) = resnumber(residue(at))
resnumber(dis_at::DisorderedAtom) = resnumber(defaultatom(dis_at))
resnumber(res::Residue) = res.number
resnumber(dis_res::DisorderedResidue) = resnumber(defaultresidue(dis_res))

"""
    inscode(at)
    inscode(res)

Get the insertion code of an `AbstractAtom` or `AbstractResidue` as a `Char`.
"""
inscode(at::Atom) = inscode(residue(at))
inscode(dis_at::DisorderedAtom) = inscode(defaultatom(dis_at))
inscode(res::Residue) = res.ins_code
inscode(dis_res::DisorderedResidue) = inscode(defaultresidue(dis_res))

"""
    resid(res; full=true)

Get a descriptive residue ID `String` for an `AbstractAtom` or
`AbstractResidue`.

Format is residue number then insertion code with "H_" in front for hetero
residues.
If `full` equals `true` the chain ID is also added after a colon.
Examples are "50A", "H_20" and "10:A".
"""
function resid(res::AbstractResidue; full::Bool=false)
    if ishetero(res)
        if full
            if inscode(res) == ' '
                return "H_$(resnumber(res)):$(chainid(res))"
            else
                return "H_$(resnumber(res))$(inscode(res)):$(chainid(res))"
            end
        else
            if inscode(res) == ' '
                return "H_$(resnumber(res))"
            else
                return "H_$(resnumber(res))$(inscode(res))"
            end
        end
    else
        if full
            if inscode(res) == ' '
                return "$(resnumber(res)):$(chainid(res))"
            else
                return "$(resnumber(res))$(inscode(res)):$(chainid(res))"
            end
        else
            if inscode(res) == ' '
                return "$(resnumber(res))"
            else
                return "$(resnumber(res))$(inscode(res))"
            end
        end
    end
end

resid(at::AbstractAtom; full::Bool=false) = resid(residue(at), full=full)

# Shortcut from atomic properties
function resid(hetatm::Bool, resnum::Int, inscode::Char)
    if hetatm
        if inscode == ' '
            return "H_$resnum"
        else
            return "H_$resnum$inscode"
        end
    else
        if inscode == ' '
            return "$resnum"
        else
            return "$resnum$inscode"
        end
    end
end

"""
    atomnames(res; strip=true)

Get the sorted list of `AbstractAtom`s in an `AbstractResidue`.
`strip` determines whether surrounding whitespace is stripped.
"""
function atomnames(res::Residue; strip::Bool=true)
    if strip
        return Base.strip.(res.atom_list)
    else
        return res.atom_list
    end
end

atomnames(dis_res::DisorderedResidue; strip::Bool=true) = atomnames(defaultresidue(dis_res), strip=strip)

"""
    atoms(res)

Return the dictionary of `AbstractAtom`s in an `AbstractResidue`.
"""
atoms(res::Residue) = res.atoms
atoms(dis_res::DisorderedResidue) = atoms(defaultresidue(dis_res))

"""
    isdisorderedres(res)

Determine if an `AbstractResidue` is a `DisorderedResidue`, i.e. there are
multiple residue names with the same residue ID.
"""
isdisorderedres(::Residue) = false
isdisorderedres(::DisorderedResidue) = true

"""
    disorderedres(dis_res, res_name)

Return the `Residue` in a `DisorderedResidue` with a given residue name.
"""
disorderedres(dis_res::DisorderedResidue, res_name::AbstractString) = dis_res.names[res_name]

"""
    defaultresname(dis_res)

Get the name of the default `Residue` in a `DisorderedResidue` as a `String`.

The default is the first name read in.
"""
defaultresname(dis_res::DisorderedResidue) = dis_res.default

"""
    defaultresidue(dis_res)

Return the default `Residue` in a `DisorderedResidue`.

The default is the first name read in.
"""
defaultresidue(dis_res::DisorderedResidue) = dis_res.names[defaultresname(dis_res)]

"""
    resnames(dis_res)

Get the residue names in an `AbstractResidue` as a `Vector{String}`.

For a `DisorderedResidue` there will be multiple residue names - in this case
the default residue name is placed first, then the others are ordered
alphabetically.
"""
function resnames(dis_res::DisorderedResidue)
    return sort(collect(keys(dis_res.names)),
        lt= (res_name_one, res_name_two) ->
            (isless(res_name_one, res_name_two) && res_name_two != defaultresname(dis_res)) ||
            res_name_one == defaultresname(dis_res)
    )
end

resnames(res::Residue) = [resname(res, strip=false)]

# Constructor acts as a setter for the default residue name
function DisorderedResidue(dis_res::DisorderedResidue, default::AbstractString)
    if !(default in resnames(dis_res))
        throw(ArgumentError("The new default residue name \"$default\" must be present in the residue"))
    end
    return DisorderedResidue(dis_res.names, default)
end

"""
    chain(at)
    chain(res)

Return the `Chain` that an `AbstractAtom` or `AbstractResidue` belongs to.
"""
chain(at::Atom) = chain(residue(at))
chain(dis_at::DisorderedAtom) = chain(defaultatom(dis_at))
chain(res::Residue) = res.chain
chain(dis_res::DisorderedResidue) = chain(defaultresidue(dis_res))
chain(ch::Chain) = ch

"""
    chainid(el)

Get the chain ID of an `AbstractAtom`, `AbstractResidue` or `Chain` as a
`String`.
"""
chainid(el::Union{AbstractResidue, AbstractAtom}) = chainid(chain(el))
chainid(ch::Chain) = ch.id

"""
    chainid!(ch, id)
    chainid!(res, id)

Set the chain ID of a `Chain` or an `AbstractResidue` to a new `String`.

If a chain with this ID already exists, it will be removed from its current
chain and added to that chain. If a chain with this ID does not exist, a new
chain will be added to the model and this residue will be added to it. If
moving this residue from a chain to a new chain leaves the old chain without
residues, the old chain will be removed from the `Model`.
"""
function chainid!(ch::Chain, id::String)
    if haskey(ch.model.chains, id)
        throw(PDBConsistencyError("invalid ID $id, the model already has a chain with this ID"))
    end

    old_id = ch.id
    ch.id = id
    ch.model.chains[id] = ch
    delete!(ch.model.chains, old_id)
end

function chainid!(res::AbstractResidue, id::String)
    current_chain = res.chain
    model_chains = current_chain.model.chains

    # Find the currently-assigned resid, which may not have been created from the resid function
    current_resid = findfirst(isequal(res), current_chain.residues)

    if id in keys(model_chains)
        if haskey(model_chains[id].residues, current_resid) && model_chains[id].residues[current_resid] != res
            throw(PDBConsistencyError("a residue with ID $current_resid already exists in chain $id, cannot copy this residue there"))
        end
        model_chains[id].residues[resid(res)] = res
    else
        model_chains[id] = Chain(id, [], Dict(current_resid => res), current_chain.model)
    end
    res.chain = model_chains[id]

    # Remove the residue from its current chain
    delete!(current_chain.residues, current_resid)
    if isempty(current_chain.residues)
        delete!(model_chains, current_chain.id)
    end

    fixlists!(structure(res))
end

"""
    resids(ch)

Get the sorted list of `AbstractResidue`s in a `Chain`.
"""
resids(ch::Chain) = ch.res_list

"""
    residues(ch)

Return the dictionary of `AbstractResidue`s in a `Chain`.
"""
residues(ch::Chain) = ch.residues

"""
    model(el)

Return the `Model` that an `AbstractAtom`, `AbstractResidue` or `Chain` belongs
to.
"""
model(at::Atom) = model(chain(at))
model(dis_at::DisorderedAtom) = model(defaultatom(dis_at))
model(res::Residue) = model(chain(res))
model(dis_res::DisorderedResidue) = model(defaultresidue(dis_res))
model(ch::Chain) = ch.model
model(mo::Model) = mo

"""
    modelnumber(el)

Get the model number of a `Model`, `Chain`, `AbstractResidue` or `AbstractAtom`
as an `Int`.
"""
modelnumber(mo::Model) = mo.number
modelnumber(el::Union{Chain, AbstractResidue, AbstractAtom}) = modelnumber(model(el))

"""
    chainids(model)
    chainids(struc)

Get the sorted chain IDs of the chains in a `Model`, or the default `Model` of a
`MolecularStructure`, as a `Vector{String}`.
"""
chainids(mo::Model) = chainid.(sort(collect(values(chains(mo)))))

function chainids(struc::MolecularStructure)
    if countmodels(struc) > 0
        return chainids(defaultmodel(struc))
    else
        return String[]
    end
end

"""
    chains(model)
    chains(struc)

Return the dictionary of `Chain`s in a `Model`, or the default `Model` of a
`MolecularStructure`.
"""
chains(mo::Model) = mo.chains

function chains(struc::MolecularStructure)
    if countmodels(struc) > 0
        return chains(defaultmodel(struc))
    else
        return Dict{String, Chain}()
    end
end

"""
    structure(el)

Return the `MolecularStructure` that an `AbstractAtom`, `AbstractResidue`, `Chain`
or `Model` belongs to.
"""
structure(at::Atom) = structure(model(at))
structure(dis_at::DisorderedAtom) = structure(defaultatom(dis_at))
structure(res::Residue) = structure(model(res))
structure(dis_res::DisorderedResidue) = structure(defaultresidue(dis_res))
structure(ch::Chain) = structure(model(ch))
structure(mo::Model) = mo.structure
structure(struc::MolecularStructure) = struc

"""
    structurename(el)

Get the name of the `MolecularStructure` that a `StructuralElement` belongs to as
a `String`.
"""
structurename(el::Union{Model, Chain, AbstractResidue, AbstractAtom}) = structurename(structure(el))
structurename(struc::MolecularStructure) = struc.name

"""
    modelnumbers(struc)

Get the sorted model numbers from a `MolecularStructure` as a `Vector{Int}`.
"""
function modelnumbers(struc::MolecularStructure)
    return modelnumber.(sort(collect(values(models(struc)))))
end

"""
    models(struc)

Return the dictionary of `Model`s in a `MolecularStructure`.
"""
models(struc::MolecularStructure) = struc.models

"""
    defaultmodel(struc)

Get the default `Model` in a `MolecularStructure`.

This is the `Model` with the lowest model number.
"""
defaultmodel(struc::MolecularStructure) = first(sort(collect(values(models(struc)))))

# Sort lists of elements

# Sort atoms by serial
function Base.isless(at_one::AbstractAtom, at_two::AbstractAtom)
    return isless(serial(at_one), serial(at_two))
end

# Sort residues by chain, then resnumber, then ins code, then hetero
function Base.isless(res_one::AbstractResidue, res_two::AbstractResidue)
    if isless(chain(res_one), chain(res_two))
        return true
    elseif chainid(res_one) == chainid(res_two)
        if isless(resnumber(res_one), resnumber(res_two))
            return true
        elseif resnumber(res_one) == resnumber(res_two)
            if isless(inscode(res_one), inscode(res_two))
                return true
            elseif inscode(res_one) == inscode(res_two)
                if !ishetero(res_one) && ishetero(res_two)
                    return true
                end
            end
        end
    end
    return false
end

# Sort chains by chain ID
# Ordering is character sorting of subsequent letters except the empty chain ID
#   comes last
function Base.isless(ch_one::Chain, ch_two::Chain)
    # Deal with usual case of single letter comparison quickly
    if length(chainid(ch_one)) == 1 && length(chainid(ch_two)) == 1 &&
            chainid(ch_one) != " " && chainid(ch_two) != " "
        return Int(chainid(ch_one)[1]) < Int(chainid(ch_two)[1])
    end
    chid_one = strip(chainid(ch_one))
    chid_two = strip(chainid(ch_two))
    if length(chid_two) == 0
        return length(chid_one) > 0
    elseif length(chid_one) == 0
        return false
    elseif length(chid_one) != length(chid_two)
        return length(chid_one) < length(chid_two)
    end
    for (i, c) in enumerate(chid_one)
        if c == chid_two[i]
            continue
        else
            return Int(c) < Int(chid_two[i])
        end
    end
    return false
end

# Sort models by model number
function Base.isless(mo_one::Model, mo_two::Model)
    return isless(modelnumber(mo_one), modelnumber(mo_two))
end

"""
    sequentialresidues(res_first, res_second)

Determine if the second residue follows the first in sequence.

For this to be `true` the residues need to have the same chain ID, both need to
be standard/hetero residues and the residue number of the second needs to be one
greater than that of the first (or the residue numbers the same and the
insertion code of the second greater than the first).
"""
function sequentialresidues(res_first::AbstractResidue, res_second::AbstractResidue)
    if chainid(res_second) == chainid(res_first) &&
            ishetero(res_second) == ishetero(res_first)
        if resnumber(res_second) == resnumber(res_first) + 1
            return true
        elseif resnumber(res_second) == resnumber(res_first) &&
                inscode(res_second) > inscode(res_first)
            return true
        end
    end
    return false
end

# Iterators to yield sub elements when looping over an element

# Iterating over a MolecularStructure yields Models
Base.length(struc::MolecularStructure) = length(modelnumbers(struc))
Base.eltype(::Type{MolecularStructure}) = Model
function Base.iterate(struc::MolecularStructure, state=1)
    state <= length(struc) ? (struc[modelnumbers(struc)[state]], state + 1) : nothing
end

# Iterating over a Model yields Chains
Base.length(mo::Model) = length(chainids(mo))
Base.eltype(::Type{Model}) = Chain
function Base.iterate(mo::Model, state=1)
    state <= length(mo) ? (mo[chainids(mo)[state]], state + 1) : nothing
end

# Iterating over a Chain yields AbstractResidues
Base.length(ch::Chain) = length(resids(ch))
Base.eltype(::Type{Chain}) = AbstractResidue
function Base.iterate(ch::Chain, state=1)
    state <= length(ch) ? (ch[resids(ch)[state]], state + 1) : nothing
end

# Iterating over a Residue yields AbstractAtoms
Base.length(res::Residue) = length(atomnames(res, strip=false))
Base.eltype(::Type{Residue}) = AbstractAtom
function Base.iterate(res::Residue, state=1)
    state <= length(res) ? (res[atomnames(res, strip=false)[state]], state + 1) : nothing
end

# Iterating over a DisorderedResidue yields AbstractAtoms
# This is not necessarily intuitive, it may be expected to yield Residues
# However this way iterating over an AbstractResidue always yields AbstractAtoms
# Use collectresidues with the expand_disordered flag to unroll disordered
#   residues
Base.length(dis_res::DisorderedResidue) = length(atomnames(dis_res, strip=false))
Base.eltype(::Type{DisorderedResidue}) = AbstractAtom
function Base.iterate(dis_res::DisorderedResidue, state=1)
    if state <= length(dis_res)
        (defaultresidue(dis_res)[atomnames(dis_res, strip=false)[state]], state + 1)
    else
        nothing
    end
end

# Iterating over an Atom returns itself
# This is not necessarily intuitive, it may be expected to not be an iterator
# However this way iterating over an AbstractAtom always yields Atoms
Base.length(::Atom) = 1
Base.eltype(::Type{Atom}) = Atom
function Base.iterate(at::Atom, state=1)
    state <= 1 ? (at, state + 1) : nothing
end

# Iterating over a DisorderedAtom yields Atoms
Base.length(dis_at::DisorderedAtom) = length(altlocids(dis_at))
Base.eltype(::Type{DisorderedAtom}) = Atom
function Base.iterate(dis_at::DisorderedAtom, state=1)
    state <= length(dis_at) ? (dis_at[altlocids(dis_at)[state]], state + 1) : nothing
end

# Collection methods defined separately to iteration for speed
"""
    collect(el)

Returns a `Vector` of the sub-elements in a `StructuralElementOrList`, e.g.
`AbstractAtom`s in a `Residue` or `AbstractResidue`s in a `Chain`.
"""
Base.collect(struc::MolecularStructure) = [struc[mn] for mn in modelnumbers(struc)]
Base.collect(mo::Model) = [mo[cn] for cn in chainids(mo)]
Base.collect(ch::Chain) = AbstractResidue[ch[rn] for rn in resids(ch)]
Base.collect(res::Residue) = AbstractAtom[res.atoms[an] for an in atomnames(res, strip=false)]
function Base.collect(dis_res::DisorderedResidue)
    res = defaultresidue(dis_res)
    return AbstractAtom[res.atoms[an] for an in atomnames(res, strip=false)]
end
Base.collect(at::Atom) = [at]
Base.collect(dis_at::DisorderedAtom) = [dis_at[al] for al in altlocids(dis_at)]

"""
    collectmodels(el)

Returns a `Vector` of the models in a `StructuralElementOrList`.

Additional arguments are model selector functions - only models that return
`true` from all the functions are retained.
"""
collectmodels(struc::MolecularStructure) = collect(struc)

collectmodels(mo::Model) = [mo]

collectmodels(el::Union{Chain, AbstractResidue, AbstractAtom}) = [model(el)]

collectmodels(mos::AbstractVector{Model}) = mos

function collectmodels(els::AbstractVector{<:Union{Chain, AbstractResidue, AbstractAtom}})
    mo_list = Model[]
    for el in els
        if !(model(el) in mo_list)
            push!(mo_list, model(el))
        end
    end
    return mo_list
end

# One selector explicitly defined to prevent this being called without selectors
function collectmodels(el::StructuralElementOrList,
                    model_selector::Function,
                    model_selectors::Function...)
    return applyselectors(collectmodels(el), model_selector, model_selectors...)
end

"""
    countmodels(el)

Return the number of `Model`s in a `StructuralElementOrList` as an `Int`.

Additional arguments are model selector functions - only models that return
`true` from all the functions are counted.
"""
function countmodels(el::StructuralElementOrList, model_selectors::Function...)
    return length(collectmodels(el, model_selectors...))
end

countmodels(struc::MolecularStructure) = length(struc)

"""
    collectchains(el)

Returns a `Vector` of the chains in a `StructuralElementOrList`.

Additional arguments are chain selector functions - only chains that return
`true` from all the functions are retained.
"""
function collectchains(struc::MolecularStructure)
    if countmodels(struc) > 0
        return collectchains(defaultmodel(struc))
    else
        return Chain[]
    end
end

collectchains(mo::Model) = collect(mo)

collectchains(ch::Chain) = [ch]

collectchains(el::Union{AbstractResidue, AbstractAtom}) = [chain(el)]

function collectchains(mos::AbstractVector{Model})
    ch_list = Chain[]
    for mo in mos
        append!(ch_list, collectchains(mo))
    end
    return ch_list
end

collectchains(chs::AbstractVector{Chain}) = chs

function collectchains(els::AbstractVector{<:Union{AbstractResidue, AbstractAtom}})
    ch_list = Chain[]
    for el in els
        if !(chain(el) in ch_list)
            push!(ch_list, chain(el))
        end
    end
    return ch_list
end

function collectchains(el::StructuralElementOrList,
                    chain_selector::Function,
                    chain_selectors::Function...)
    return applyselectors(collectchains(el), chain_selector, chain_selectors...)
end

"""
    countchains(el)

Return the number of `Chain`s in a `StructuralElementOrList` as an `Int`.

Additional arguments are chain selector functions - only chains that return
`true` from all the functions are counted.
"""
function countchains(el::StructuralElementOrList, chain_selectors::Function...)
    return length(collectchains(el, chain_selectors...))
end

countchains(mo::Model) = length(mo)

"""
    collectresidues(el)

Returns a `Vector` of the residues in a `StructuralElementOrList`.

Additional arguments are residue selector functions - only residues that
return `true` from all the functions are retained.
The keyword argument `expand_disordered` (default `false`) determines whether to
return all copies of disordered residues separately.
"""
function collectresidues(struc::MolecularStructure; expand_disordered::Bool=false)
    if countmodels(struc) > 0
        return collectresidues(defaultmodel(struc); expand_disordered=expand_disordered)
    else
        return AbstractResidue[]
    end
end

function collectresidues(el::Union{Model, Vector{Model}, Vector{Chain}};
                            expand_disordered::Bool=false)
    res_list = AbstractResidue[]
    for sub_el in el
        append!(res_list, collectresidues(sub_el; expand_disordered=expand_disordered))
    end
    return res_list
end

# Note output is always Vector{AbstractResidue} unless input was Vector{Residue}
#   or Vector{DisorderedResidue}, in which case output is same type as input
#   type
function collectresidues(el::Union{Chain, Vector{<:AbstractResidue}};
                            expand_disordered::Bool=false)
    if expand_disordered
        res_list = AbstractResidue[]
        for res in el
            if isa(res, Residue)
                push!(res_list, res)
            else
                append!(res_list, collectresidues(res; expand_disordered=true))
            end
        end
        return res_list
    else
        return isa(el, Chain) ? collect(el) : el
    end
end

collectresidues(res::Residue; expand_disordered::Bool=false) = AbstractResidue[res]

function collectresidues(dis_res::DisorderedResidue; expand_disordered::Bool=false)
    if expand_disordered
        return AbstractResidue[disorderedres(dis_res, res_name) for res_name in resnames(dis_res)]
    else
        return AbstractResidue[dis_res]
    end
end

collectresidues(at::AbstractAtom; expand_disordered::Bool=false) = AbstractResidue[residue(at)]

function collectresidues(at_list::AbstractVector{<:AbstractAtom}; expand_disordered::Bool=false)
    res_list = AbstractResidue[]
    for at in at_list
        if !(residue(at) in res_list)
            push!(res_list, residue(at))
        end
    end
    return res_list
end

function collectresidues(el::StructuralElementOrList,
                    residue_selector::Function,
                    residue_selectors::Function...;
                    expand_disordered::Bool=false)
    return collectresidues(applyselectors(collectresidues(el), residue_selector,
                            residue_selectors...); expand_disordered=expand_disordered)
end

"""
    countresidues(el)

Return the number of residues in a `StructuralElementOrList` as an `Int`.

Additional arguments are residue selector functions - only residues that
return `true` from all the functions are counted.
The keyword argument `expand_disordered` (default `false`) determines whether to
return all copies of disordered residues separately.
"""
function countresidues(el::StructuralElementOrList,
                        residue_selectors::Function...;
                        expand_disordered::Bool=false)
    return length(collectresidues(el, residue_selectors...;
                                    expand_disordered=expand_disordered))
end

"""
    collectatoms(el)

Returns a `Vector` of the atoms in a `StructuralElementOrList`.

Additional arguments are atom selector functions - only atoms that return
`true` from all the functions are retained.
The keyword argument `expand_disordered` (default `false`) determines whether to
return all copies of disordered atoms separately.
"""
function collectatoms(struc::MolecularStructure; expand_disordered::Bool=false)
    if countmodels(struc) > 0
        return collectatoms(defaultmodel(struc); expand_disordered=expand_disordered)
    else
        return AbstractAtom[]
    end
end

function collectatoms(el::Union{Model, Chain, Vector{Model}, Vector{Chain},
                                Vector{<:AbstractResidue}};
                        expand_disordered::Bool=false)
    at_list = AbstractAtom[]
    for sub_el in el
        append!(at_list, collectatoms(sub_el; expand_disordered=expand_disordered))
    end
    return at_list
end

# Note output is always Vector{AbstractAtom} unless input was Vector{Atom} or
#   Vector{DisorderedAtom}, in which case output is same type as input type
function collectatoms(el::Union{Residue, Vector{<:AbstractAtom}};
                            expand_disordered::Bool=false)
    if expand_disordered
        at_list = AbstractAtom[]
        for at in el
            if isa(at, Atom)
                push!(at_list, at)
            else
                append!(at_list, collectatoms(at; expand_disordered=true))
            end
        end
        return at_list
    else
        return isa(el, Residue) ? collect(el) : el
    end
end

function collectatoms(dis_res::DisorderedResidue; expand_disordered::Bool=false)
    if expand_disordered
        return collectatoms(collectresidues(dis_res; expand_disordered=true);
                            expand_disordered=true)
    else
        return collectatoms(defaultresidue(dis_res))
    end
end

collectatoms(at::Atom; expand_disordered::Bool=false) = AbstractAtom[at]

function collectatoms(dis_at::DisorderedAtom; expand_disordered::Bool=false)
    if expand_disordered
        return AbstractAtom[at for at in dis_at]
    else
        return AbstractAtom[dis_at]
    end
end

function collectatoms(el::StructuralElementOrList,
                        atom_selector::Function,
                        atom_selectors::Function...;
                        expand_disordered::Bool=false)
    return collectatoms(applyselectors(collectatoms(el), atom_selector, atom_selectors...);
                        expand_disordered=expand_disordered)
end

"""
    countatoms(el)

Return the number of atoms in a `StructuralElementOrList` as an `Int`.

Additional arguments are atom selector functions - only atoms that return
`true` from all the functions are counted.
The keyword argument `expand_disordered` (default `false`) determines whether to
return all copies of disordered atoms separately.
"""
function countatoms(el::StructuralElementOrList,
                    atom_selectors::Function...;
                    expand_disordered::Bool=false)
    return length(collectatoms(el, atom_selectors...;
                                expand_disordered=expand_disordered))
end

# Add an atom represented in an AtomRecord to a Model
# Unsafe as sub-element lists are not updated (for speed)
# fixlists! should be run after all additions to update the sub-element lists
function unsafe_addatomtomodel!(mo::Model,
                    atom_rec::AtomRecord;
                    remove_disorder::Bool=false)
    # Add chain to model if necessary
    if !haskey(chains(mo), atom_rec.chain_id)
        mo[atom_rec.chain_id] = Chain(atom_rec.chain_id, mo)
    end
    ch = mo[atom_rec.chain_id]
    res_id = resid(atom_rec.het_atom, atom_rec.res_number, atom_rec.ins_code)
    # If residue does not exist in the chain, create a Residue
    if !haskey(residues(ch), res_id)
        ch[res_id] = Residue(
                    atom_rec.res_name,
                    atom_rec.res_number,
                    atom_rec.ins_code,
                    atom_rec.het_atom,
                    ch)
        res = ch[res_id]
    elseif isa(ch[res_id], Residue)
        # Residue exists in the chain and the residue names match
        # Add to that Residue
        if fullresname(ch[res_id]) == atom_rec.res_name
            res = ch[res_id]
        # Residue exists in the chain but the residue names do not match
        # Create a DisorderedResidue
        else
            ch[res_id] = DisorderedResidue(Dict(
                fullresname(ch[res_id]) => ch[res_id],
                atom_rec.res_name => Residue(
                    atom_rec.res_name,
                    atom_rec.res_number,
                    atom_rec.ins_code,
                    atom_rec.het_atom,
                    ch)
            ), fullresname(ch[res_id]))
            res = disorderedres(ch[res_id], atom_rec.res_name)
        end
    else
        # DisorderedResidue exists in the chain and the residue names match
        # Add to that DisorderedResidue
        if atom_rec.res_name in resnames(ch[res_id])
            res = disorderedres(ch[res_id], atom_rec.res_name)
        # DisorderedResidue exists in the chain and the residue names do not match
        # Create a new Residue in the DisorderedResidue
        else
            ch[res_id].names[atom_rec.res_name] = Residue(
                    atom_rec.res_name,
                    atom_rec.res_number,
                    atom_rec.ins_code,
                    atom_rec.het_atom,
                    ch)
            res = disorderedres(ch[res_id], atom_rec.res_name)
        end
    end
    at = Atom(
        atom_rec.serial,
        atom_rec.atom_name,
        atom_rec.alt_loc_id,
        atom_rec.coords,
        atom_rec.occupancy,
        atom_rec.temp_factor,
        atom_rec.element,
        atom_rec.charge,
        res)
    # If atom does not exist in the residue, create an Atom
    if !haskey(atoms(res), atom_rec.atom_name)
        res[atom_rec.atom_name] = at
    # Atom exists in the residue, atom names match and alt loc IDs are different
    elseif isa(res[atom_rec.atom_name], Atom) &&
            atom_rec.alt_loc_id != altlocid(res[atom_rec.atom_name])
        # If we are removing disorder and the new atom is preferred to the old one, replace the old one
        if remove_disorder &&
                choosedefaultaltlocid(at, res[atom_rec.atom_name]) == atom_rec.alt_loc_id
            res[atom_rec.atom_name] = at
        # If we are not removing disorder, create a new disordered atom container and add both atoms
        elseif !remove_disorder
            res[atom_rec.atom_name] = DisorderedAtom(Dict(
                atom_rec.alt_loc_id => at,
                altlocid(res[atom_rec.atom_name]) => res[atom_rec.atom_name]
            ), choosedefaultaltlocid(at, res[atom_rec.atom_name]))
        end
    # A disordered atom container already exists and the alt loc ID is not taken
    elseif isa(res[atom_rec.atom_name], DisorderedAtom) &&
            !(atom_rec.alt_loc_id in altlocids(res[atom_rec.atom_name]))
        # Add the new atom to the disordered atom container
        res[atom_rec.atom_name][atom_rec.alt_loc_id] = at
        # If the default alt loc requires changing, change it
        if choosedefaultaltlocid(defaultatom(res[atom_rec.atom_name]), at) != defaultaltlocid(res[atom_rec.atom_name])
            res[atom_rec.atom_name] = DisorderedAtom(
                        res[atom_rec.atom_name],
                        atom_rec.alt_loc_id)
        end
    else
        error("Two copies of the same atom have the same alternative location ID. Existing atom:\n" *
              "$(res[atom_rec.atom_name])\nNew atom record to add:\n$(atom_rec)")
    end
end

fullresname(res::Residue) = res.name

# Internal function to form ordered sub-element lists after parsing
function fixlists!(struc::MolecularStructure)
    for mo in struc
        for ch in mo
            ch.res_list = resid.(sort(collect(values(residues(ch)))))
            for res in ch
                if isa(res, Residue)
                    fixlists!(res)
                else
                    for res_name in resnames(res)
                        fixlists!(disorderedres(res, res_name))
                    end
                end
            end
        end
    end
end

function fixlists!(res::Residue)
    res.atom_list = fullatomname.(sort(collect(values(atoms(res)))))
end

fullatomname(at::Atom) = at.name
fullatomname(dis_at::DisorderedAtom) = defaultatom(dis_at).name

"""
    choosedefaultaltlocid(at_one, at_two)

Determine which of two `Atom`s representing a disorered atom better
qualifies as the default location.

The `Atom` with the highest occupancy is chosen; in the case of ties the
`Atom` with the lowest alternative location ID in alphabetical order is
chosen.
"""
function choosedefaultaltlocid(at_one::Atom, at_two::Atom)
    if occupancy(at_one) > occupancy(at_two) ||
            (occupancy(at_one) == occupancy(at_two) &&
            Int(altlocid(at_one)) < Int(altlocid(at_two)))
        return altlocid(at_one)
    else
        return altlocid(at_two)
    end
end

"Lookup table of amino acids, re-exported from BioSymbols."
const threeletter_to_aa = BioSymbols.threeletter_to_aa

# PDB file formats

"Protein Data Bank (PDB) file format."
struct PDB end

"Protein Data Bank (PDB) XML file format."
struct PDBXML end

"Protein Data Bank (PDB) mmCIF file format."
struct MMCIF end

"Protein Data Bank (PDB) MMTF file format."
struct MMTF end

"Mapping of Protein Data Bank (PDB) formats to their file extensions."
const pdbextension = Dict{Type, String}(
    PDB    => "pdb",
    PDBXML => "xml",
    MMCIF  => "cif",
    MMTF   => "mmtf",
)

"""
    generatechainid(entity_id)

Convert a positive `Integer` into a chain ID.

Goes A to Z, then AA to ZA, AB to ZB etc.
This is in line with Protein Data Bank (PDB) conventions.
"""
function generatechainid(entity_id::Integer)
    entity_id > 0 || throw(ArgumentError("Entity ID $entity_id is not positive"))
    divisor = entity_id
    out_string = ""
    while divisor > 0
        modulo = (divisor - 1) % 26
        out_string *= Char(65 + modulo)
        divisor = Int(floor((divisor - modulo) / 26))
    end
    return out_string
end

"""
    MMTFDict(filepath; gzip=false)
    MMTFDict(io; gzip=false)
    MMTFDict()

A Macromolecular Transmission Format (MMTF) dictionary.
Use of the dictionary requires the MMTF.jl package to be imported with
`import MMTF as MMTFPkg`.
Can be accessed using similar functions to a standard `Dict`.
Keys are field names as a `String` and values are various types.
To directly access the underlying dictionary of `MMTFDict` `d`, use
`d.dict`.
Call `MMTFDict` with a filepath or stream to read the dictionary from that
source.
The keyword argument `gzip` (default `false`) determines if the file is gzipped.
"""
struct MMTFDict <: AbstractDict{String, Any}
    dict::Dict{String, Any}
end

"""
    writemmtf(output, element, atom_selectors...; gzip=false)
    writemmtf(output, mmtf_dict; gzip=false)

Write a `StructuralElementOrList` or a `MMTFDict` to a MMTF file or output
stream.

Requires the MMTF.jl package to be imported with `import MMTF as MMTFPkg`.
Atom selector functions can be given as additional arguments - only atoms
that return `true` from all the functions are retained.
The keyword argument `expand_disordered` (default `true`) determines whether to
return all copies of disordered residues and atoms.
The keyword argument `gzip` (default `false`) determines if the file should be
gzipped.
"""
function writemmtf end

# Descriptive showing of elements on a single line

Base.show(io::IO, struc::MolecularStructure) = print(io,
    "MolecularStructure ",
    structurename(struc) != "" ? "$(structurename(struc)) " : "",
    "with $(countmodels(struc)) models, ",
    countchains(struc) != 0 ? "$(countchains(struc)) chains ($(join(chainids(struc), ","))), " : "0 chains, ",
    "$(countresidues(struc, standardselector)) residues, ",
    "$(countatoms(struc)) atoms"
)

Base.show(io::IO, mo::Model) = print(io,
    "Model $(modelnumber(mo)) with ",
    countchains(mo) != 0 ? "$(countchains(mo)) chains ($(join(chainids(mo), ","))), " : "0 chains, ",
    "$(countresidues(mo, standardselector)) residues, ",
    "$(countatoms(mo)) atoms"
)

Base.show(io::IO, ch::Chain) = print(io,
    "Chain $(chainid(ch)) with ",
    "$(countresidues(ch, standardselector)) residues, ",
    "$(countresidues(ch, heteroselector)) other molecules, ",
    "$(countatoms(ch)) atoms"
)

Base.show(io::IO, res::Residue) = print(io,
    "Residue $(resid(res, full=true)) with ",
    "name $(resname(res)), ",
    "$(countatoms(res)) atoms"
)

Base.show(io::IO, dis_res::DisorderedResidue) = print(io,
    "DisorderedResidue $(resid(dis_res, full=true)) with ",
    "names $(join(resnames(dis_res), ","))"
)

Base.show(io::IO, at::Atom) = print(io,
    "Atom $(atomname(at)) with ",
    "serial $(serial(at)), ",
    "coordinates $(coords(at))",
    altlocid(at) != ' ' ? ", alt loc ID $(altlocid(at))" : ""
)

Base.show(io::IO, dis_at::DisorderedAtom) = print(io,
    "DisorderedAtom $(atomname(dis_at)) with ",
    "alt loc IDs $(join(altlocids(dis_at), ","))"
)

Base.show(io::IO, at_rec::AtomRecord) = print(io,
    "AtomRecord $(strip(at_rec.atom_name)) with ",
    "serial $(at_rec.serial), ",
    "coordinates $(at_rec.coords)",
    at_rec.alt_loc_id != ' ' ? ", alt loc ID $(at_rec.alt_loc_id)" : ""
)
