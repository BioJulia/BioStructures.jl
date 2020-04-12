export
    StructuralElement,
    AbstractAtom,
    Atom,
    DisorderedAtom,
    AbstractResidue,
    Residue,
    DisorderedResidue,
    Chain,
    Model,
    ProteinStructure,
    AtomRecord,
    StructuralElementOrList,
    serial,
    atomname,
    altlocid,
    x,
    x!,
    y,
    y!,
    z,
    z!,
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
    applyselectors,
    applyselectors!,
    collectmodels,
    countmodels,
    collectchains,
    countchains,
    collectresidues,
    countresidues,
    collectatoms,
    countatoms,
    choosedefaultaltlocid,
    standardselector,
    heteroselector,
    atomnameselector,
    calphaatomnames,
    calphaselector,
    cbetaatomnames,
    cbetaselector,
    backboneatomnames,
    backboneselector,
    heavyatomselector,
    resnameselector,
    waterresnames,
    waterselector,
    notwaterselector,
    disorderselector,
    hydrogenselector,
    allselector,
    AminoAcidSequence,
    pairalign,
    DataFrame

"A macromolecular structural element."
abstract type StructuralElement end

"""
An atom that is part of a macromolecule - either an `Atom` or a
`DisorderedAtom`.
"""
abstract type AbstractAtom <: StructuralElement end

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
struct Residue <: AbstractResidue
    name::String
    number::Int
    ins_code::Char
    het_res::Bool # Does the residue consist of hetatoms?
    atom_list::Vector{String}
    atoms::Dict{String, AbstractAtom}
    chain::StructuralElement
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
struct Chain <: StructuralElement
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
struct ProteinStructure <: StructuralElement
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

ProteinStructure(name::AbstractString) = ProteinStructure(name, Dict())

ProteinStructure() = ProteinStructure("")

Model(number::Integer, struc::ProteinStructure) = Model(number, Dict(), struc)

Model(number::Integer) = Model(number, ProteinStructure())

Model() = Model(1)

Chain(id::AbstractString, mod::Model) = Chain(id, [], Dict(), mod)
Chain(id::Char, mod::Model) = Chain(string(id), [], Dict(), mod)

Chain(id::Union{AbstractString, Char}) = Chain(id, Model())

function Residue(name::AbstractString,
                number::Integer,
                ins_code::Char,
                het_res::Bool,
                ch::Chain)
    return Residue(name, number, ins_code, het_res, [], Dict(), ch)
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

# Accessing a Residue with an AbstractString returns the AbstractAtom with that
#   atom name
Base.getindex(res::Residue, atom_name::AbstractString) = findatombyname(res, atom_name)

function Base.setindex!(res::Residue, at::AbstractAtom, atom_name::AbstractString)
    res.atoms[atom_name] = at
    return res
end

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

# Accessing a Model with a Char or AbstractString returns the Chain with that
#   chain ID
Base.getindex(mod::Model, ch_id::AbstractString) = mod.chains[ch_id]
Base.getindex(mod::Model, ch_id::Char) = mod.chains[string(ch_id)]

function Base.setindex!(mod::Model, ch::Chain, ch_id::AbstractString)
    mod.chains[ch_id] = ch
    return mod
end

function Base.setindex!(mod::Model, ch::Chain, ch_id::Char)
    return Base.setindex!(mod, ch, string(ch_id))
end

# Accessing a ProteinStructure with an Integer returns the Model with that model
#   number
Base.getindex(struc::ProteinStructure, mod_n::Integer) = struc.models[mod_n]

function Base.setindex!(struc::ProteinStructure, mod::Model, mod_n::Integer)
    struc.models[mod_n] = mod
    return struc
end

# Accessing a ProteinStructure with a Char returns the Chain with that chain ID
#   on the default model
Base.getindex(struc::ProteinStructure, ch_id::AbstractString) = defaultmodel(struc)[ch_id]
Base.getindex(struc::ProteinStructure, ch_id::Char) = defaultmodel(struc)[string(ch_id)]

function Base.setindex!(struc::ProteinStructure, ch::Chain, ch_id::AbstractString)
    defaultmodel(struc)[ch_id] = ch
    return struc
end

function Base.setindex!(struc::ProteinStructure, ch::Chain, ch_id::Char)
    return Base.setindex!(struc, ch, string(ch_id))
end

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

Get the x coordinate of an `AbstractAtom` as a `Float64`.
"""
x(at::Atom) = at.coords[1]
x(dis_at::DisorderedAtom) = x(defaultatom(dis_at))

"""
    x!(at, val)

Set the x coordinate of an `AbstractAtom` to `val`.

For `DisorderedAtom`s only the default atom is updated.
"""
x!(at::Atom, x::Real) = (at.coords[1] = x; at)
x!(dis_at::DisorderedAtom, x::Real) = x!(defaultatom(dis_at), x)

"""
    y(at)

Get the y coordinate of an `AbstractAtom` as a `Float64`.
"""
y(at::Atom) = at.coords[2]
y(dis_at::DisorderedAtom) = y(defaultatom(dis_at))

"""
    y!(at, val)

Set the y coordinate of an `AbstractAtom` to `val`.

For `DisorderedAtom`s only the default atom is updated.
"""
y!(at::Atom, y::Real) = (at.coords[2] = y; at)
y!(dis_at::DisorderedAtom, y::Real) = y!(defaultatom(dis_at), y)

"""
    z(at)

Get the z coordinate of an `AbstractAtom` as a `Float64`.
"""
z(at::Atom) = at.coords[3]
z(dis_at::DisorderedAtom) = z(defaultatom(dis_at))

"""
    z!(at, val)

Set the z coordinate of an `AbstractAtom` to `val`.

For `DisorderedAtom`s only the default atom is updated.
"""
z!(at::Atom, z::Real) = (at.coords[3] = z; at)
z!(dis_at::DisorderedAtom, z::Real) = z!(defaultatom(dis_at), z)

"""
    coords(at)

Get the atomic coordinates of an `AbstractAtom` as a `Vector{Float64}`.
"""
coords(at::Atom) = at.coords
coords(dis_at::DisorderedAtom) = coords(defaultatom(dis_at))

"""
    coords!(at, new_coords)

Set the coordinates of an `AbstractAtom` to a `Vector` of 3 numbers.

For `DisorderedAtom`s only the default atom is updated.
"""
function coords!(at::Atom, new_coords::Vector{<:Real})
    if length(new_coords) != 3
        throw(ArgumentError("3 coordinates must be given"))
    end
    x!(at, new_coords[1])
    y!(at, new_coords[2])
    z!(at, new_coords[3])
    return at
end

function coords!(dis_at::DisorderedAtom, coords::Vector{<:Real})
    coords!(defaultatom(dis_at), coords)
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
model(mod::Model) = mod

"""
    modelnumber(el)

Get the model number of a `Model`, `Chain`, `AbstractResidue` or `AbstractAtom`
as an `Int`.
"""
modelnumber(mod::Model) = mod.number
modelnumber(el::Union{Chain, AbstractResidue, AbstractAtom}) = modelnumber(model(el))

"""
    chainids(mod)
    chainids(struc)

Get the sorted chain IDs of the chains in a `Model`, or the default `Model` of a
`ProteinStructure`, as a `Vector{String}`.
"""
chainids(mod::Model) = chainid.(sort(collect(values(chains(mod)))))

function chainids(struc::ProteinStructure)
    if countmodels(struc) > 0
        return chainids(defaultmodel(struc))
    else
        return String[]
    end
end

"""
    chains(mod)
    chains(struc)

Return the dictionary of `Chain`s in a `Model`, or the default `Model` of a
`ProteinStructure`.
"""
chains(mod::Model) = mod.chains

function chains(struc::ProteinStructure)
    if countmodels(struc) > 0
        return chains(defaultmodel(struc))
    else
        return Dict{String, Chain}()
    end
end

"""
    structure(el)

Return the `ProteinStructure` that an `AbstractAtom`, `AbstractResidue`, `Chain`
or `Model` belongs to.
"""
structure(at::Atom) = structure(model(at))
structure(dis_at::DisorderedAtom) = structure(defaultatom(dis_at))
structure(res::Residue) = structure(model(res))
structure(dis_res::DisorderedResidue) = structure(defaultresidue(dis_res))
structure(ch::Chain) = structure(model(ch))
structure(mod::Model) = mod.structure
structure(struc::ProteinStructure) = struc

"""
    structurename(el)

Get the name of the `ProteinStructure` that a `StructuralElement` belongs to as
a `String`.
"""
structurename(el::Union{Model, Chain, AbstractResidue, AbstractAtom}) = structurename(structure(el))
structurename(struc::ProteinStructure) = struc.name

"""
    modelnumbers(struc)

Get the sorted model numbers from a `ProteinStructure` as a `Vector{Int}`.
"""
function modelnumbers(struc::ProteinStructure)
    return modelnumber.(sort(collect(values(models(struc)))))
end

"""
    models(struc)

Return the dictionary of `Model`s in a `ProteinStructure`.
"""
models(struc::ProteinStructure) = struc.models

"""
    defaultmodel(struc)

Get the default `Model` in a `ProteinStructure`.

This is the `Model` with the lowest model number.
"""
defaultmodel(struc::ProteinStructure) = first(sort(collect(values(models(struc)))))

# Sort lists of elements

# Sort atoms by serial
function Base.isless(at_one::AbstractAtom, at_two::AbstractAtom)
    return isless(serial(at_one), serial(at_two))
end

# Sort residues by chain, then hetero, then resnumber, then ins code
function Base.isless(res_one::AbstractResidue, res_two::AbstractResidue)
    if isless(chain(res_one), chain(res_two))
        return true
    elseif chainid(res_one) == chainid(res_two)
        if !ishetero(res_one) && ishetero(res_two)
            return true
        elseif ishetero(res_one) == ishetero(res_two)
            if isless(resnumber(res_one), resnumber(res_two))
                return true
            elseif resnumber(res_one) == resnumber(res_two)
                if isless(inscode(res_one), inscode(res_two))
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
function Base.isless(mod_one::Model, mod_two::Model)
    return isless(modelnumber(mod_one), modelnumber(mod_two))
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

# Iterating over a ProteinStructure yields Models
Base.length(struc::ProteinStructure) = length(modelnumbers(struc))
Base.eltype(::Type{ProteinStructure}) = Model
function Base.iterate(struc::ProteinStructure, state=1)
    state <= length(struc) ? (struc[modelnumbers(struc)[state]], state + 1) : nothing
end

# Iterating over a Model yields Chains
Base.length(mod::Model) = length(chainids(mod))
Base.eltype(::Type{Model}) = Chain
function Base.iterate(mod::Model, state=1)
    state <= length(mod) ? (mod[chainids(mod)[state]], state + 1) : nothing
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
Base.collect(struc::ProteinStructure) = [struc[mn] for mn in modelnumbers(struc)]
Base.collect(mod::Model) = [mod[cn] for cn in chainids(mod)]
Base.collect(ch::Chain) = AbstractResidue[ch[rn] for rn in resids(ch)]
Base.collect(res::Residue) = AbstractAtom[res.atoms[an] for an in atomnames(res, strip=false)]
function Base.collect(dis_res::DisorderedResidue)
    res = defaultresidue(dis_res)
    return AbstractAtom[res.atoms[an] for an in atomnames(res, strip=false)]
end
Base.collect(at::Atom) = [at]
Base.collect(dis_at::DisorderedAtom) = [dis_at[al] for al in altlocids(dis_at)]

"""
    applyselectors(els, selectors...)

Returns a copy of a `Vector` of `StructuralElement`s with all elements that do
not return `true` from all the selector functions removed.
"""
function applyselectors(els::Vector{<:StructuralElement}, selectors::Function...)
    new_list = copy(els)
    applyselectors!(new_list, selectors...)
    return new_list
end

"""
    applyselectors!(els, selectors...)

Removes from a `Vector` of `StructuralElement`s all elements that do not return
`true` from all the selector functions.
"""
function applyselectors!(els::Vector{<:StructuralElement}, selectors::Function...)
    for selector in selectors
        filter!(selector, els)
    end
    return els
end

"""
    collectmodels(el)

Returns a `Vector` of the models in a `StructuralElementOrList`.

Additional arguments are model selector functions - only models that return
`true` from all the functions are retained.
"""
collectmodels(struc::ProteinStructure) = collect(struc)

collectmodels(mod::Model) = [mod]

collectmodels(el::Union{Chain, AbstractResidue, AbstractAtom}) = [model(el)]

collectmodels(mods::Vector{Model}) = mods

function collectmodels(els::Vector{<:Union{Chain, AbstractResidue, AbstractAtom}})
    mod_list = Model[]
    for el in els
        if !(model(el) in mod_list)
            push!(mod_list, model(el))
        end
    end
    return mod_list
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

countmodels(struc::ProteinStructure) = length(struc)

"""
    collectchains(el)

Returns a `Vector` of the chains in a `StructuralElementOrList`.

Additional arguments are chain selector functions - only chains that return
`true` from all the functions are retained.
"""
function collectchains(struc::ProteinStructure)
    if countmodels(struc) > 0
        return collectchains(defaultmodel(struc))
    else
        return Chain[]
    end
end

collectchains(mod::Model) = collect(mod)

collectchains(ch::Chain) = [ch]

collectchains(el::Union{AbstractResidue, AbstractAtom}) = [chain(el)]

function collectchains(mods::Vector{Model})
    ch_list = Chain[]
    for mod in mods
        append!(ch_list, collectchains(mod))
    end
    return ch_list
end

collectchains(chs::Vector{Chain}) = chs

function collectchains(els::Vector{<:Union{AbstractResidue, AbstractAtom}})
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

countchains(mod::Model) = length(mod)

"""
    collectresidues(el)

Returns a `Vector` of the residues in a `StructuralElementOrList`.

Additional arguments are residue selector functions - only residues that
return `true` from all the functions are retained.
The keyword argument `expand_disordered` (default `false`) determines whether to
return all copies of disordered residues separately.
"""
function collectresidues(struc::ProteinStructure; expand_disordered::Bool=false)
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

function collectresidues(at_list::Vector{<:AbstractAtom}; expand_disordered::Bool=false)
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
    return applyselectors(collectresidues(el; expand_disordered=expand_disordered),
                            residue_selector, residue_selectors...)
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
function collectatoms(struc::ProteinStructure; expand_disordered::Bool=false)
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
    return applyselectors(collectatoms(el; expand_disordered=expand_disordered),
                            atom_selector, atom_selectors...)
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
function unsafe_addatomtomodel!(mod::Model,
                    atom_rec::AtomRecord;
                    remove_disorder::Bool=false)
    # Add chain to model if necessary
    if !haskey(chains(mod), atom_rec.chain_id)
        mod[atom_rec.chain_id] = Chain(atom_rec.chain_id, mod)
    end
    ch = mod[atom_rec.chain_id]
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
function fixlists!(struc::ProteinStructure)
    for mod in struc
        for ch in mod
            append!(ch.res_list, resid.(sort(collect(values(residues(ch))))))
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
    append!(res.atom_list, fullatomname.(sort(collect(values(atoms(res))))))
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

"""
    standardselector(at)
    standardselector(res)

Determines if an `AbstractAtom` represents a standard atom, e.g. came from
a ATOM record in a Protein Data Bank (PDB) file, or if an `AbstractResidue`
represents a standard molecule, e.g. consists of ATOM records from a PDB
file.
"""
standardselector(at::AbstractAtom) = !ishetero(at)
standardselector(res::AbstractResidue) = !ishetero(res)

"""
    heteroselector(at)
    heteroselector(res)

Determines if an `AbstractAtom` represents a hetero atom, e.g. came from a
HETATM record in a Protein Data Bank (PDB) file, or if an `AbstractResidue`
represents a hetero molecule, e.g. consists of HETATM records from a PDB
file.
"""
heteroselector(at::AbstractAtom) = ishetero(at)
heteroselector(res::AbstractResidue) = ishetero(res)

"""
    atomnameselector(at, atom_names; strip=true)

Determines if an `AbstractAtom` has its atom name in the given `Set` or
`Vector`.
`strip` determines whether surrounding whitespace is stripped from the atom
name before it is checked in the list.
"""
function atomnameselector(at::AbstractAtom,
                    atom_names::Union{Set{String}, Vector{String}};
                    strip::Bool=true)
    return atomname(at, strip=strip) in atom_names
end

"`Set` of Cα atom names."
const calphaatomnames = Set(["CA"])

"""
    calphaselector(at)

Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
Cα atom.
"""
function calphaselector(at::AbstractAtom)
    return standardselector(at) && atomnameselector(at, calphaatomnames)
end

"`Set` of Cβ atom names."
const cbetaatomnames = Set(["CB"])

"""
    cbetaselector(at)

Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
Cβ atom, or a Cα atom in glycine.
"""
function cbetaselector(at::AbstractAtom)
    return standardselector(at) &&
        (atomnameselector(at, cbetaatomnames) ||
        (resname(at, strip=false) == "GLY" && atomnameselector(at, calphaatomnames)))
end

"`Set` of protein backbone atom names."
const backboneatomnames = Set(["CA", "N", "C", "O"])

"""
    backboneselector(at)

Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
protein backbone atom.
"""
function backboneselector(at::AbstractAtom)
    return standardselector(at) && atomnameselector(at, backboneatomnames)
end

"""
    heavyatomselector(at)

Determines if an `AbstractAtom` corresponds to a heavy (non-hydrogen) atom
and is not a hetero-atom.
"""
heavyatomselector(at::AbstractAtom) = standardselector(at) && !hydrogenselector(at)

"""
    resnameselector(res, res_names)
    resnameselector(at, res_names)

Determines if an `AbstractResidue` or `AbstractAtom` has its residue name
in the given `Set` or `Vector`.
"""
function resnameselector(el::Union{AbstractResidue, AbstractAtom},
                    res_names::Union{Set{String}, Vector{String}})
    return resname(el, strip=false) in res_names
end

"`Set` of residue names corresponding to water."
const waterresnames = Set(["HOH"])

"""
    waterselector(res)
    waterselector(at)

Determines if an `AbstractResidue` or `AbstractAtom` represents a water
molecule.
"""
function waterselector(el::Union{AbstractResidue, AbstractAtom})
    return resnameselector(el, waterresnames)
end

"""
    notwaterselector(res)
    notwaterselector(at)

Determines if an `AbstractResidue` or `AbstractAtom` does not represent a
water molecule.
"""
function notwaterselector(el::Union{AbstractResidue, AbstractAtom})
    return !waterselector(el)
end

"""
    disorderselector(at)
    disorderselector(res)

Determines whether an `AbstractAtom` or `AbstractResidue` is disordered,
i.e. has multiple locations in the case of atoms or multiple residue names
(point mutants) in the case of residues.
"""
disorderselector(at::AbstractAtom) = isdisorderedatom(at)
disorderselector(res::AbstractResidue) = isdisorderedres(res)

# Either the element is H or the element field is empty, the atom name contains
#   an H and there are no letters in the atom name before H
# For example atom names "1H" and "H1" would be hydrogens but "NH1" would not
"""
    hydrogenselector(at)

Determines if an `AbstractAtom` represents hydrogen.

Uses the element field where possible, otherwise uses the atom name.
"""
function hydrogenselector(at::AbstractAtom)
    at_name = atomname(at)
    return element(at) == "H" ||
        (element(at) == "" &&
        'H' in at_name &&
        !occursin(r"[a-zA-Z]", at_name[1:findfirst(isequal('H'), at_name) - 1]))
end

"""
    allselector(at)
    allselector(res)

Trivial selector that returns `true` for any `AbstractAtom` or
`AbstractResidue`.
Use it to select all atoms or residues.
"""
allselector(at::AbstractAtom) = true
allselector(res::AbstractResidue) = true

"""
    AminoAcidSequence(el)

Return the amino acid sequence of a protein.

Additional arguments are residue selector functions - only residues that
return `true` from all the functions are retained.
The `gaps` keyword argument determines whether to add gaps to the sequence
based on missing residue numbers (default `true`).
See BioSequences.jl for more on how to use sequences.
"""
function BioSequences.AminoAcidSequence(el::Union{StructuralElement, Vector{Model},
                                    Vector{Chain}, Vector{<:AbstractAtom}},
                        residue_selectors::Function...;
                        gaps::Bool=true)
    return AminoAcidSequence(collectresidues(el, residue_selectors...); gaps=gaps)
end

function BioSequences.AminoAcidSequence(res::Vector{<:AbstractResidue}; gaps::Bool=true)
    seq = AminoAcid[]
    for i in 1:length(res)
        if haskey(BioSymbols.threeletter_to_aa, resname(res[i], strip=false))
            push!(seq, BioSymbols.threeletter_to_aa[resname(res[i], strip=false)])
        else
            push!(seq, AA_X)
        end
        # Add gaps based on missing residue numbers
        if gaps && i + 1 <= length(res) && resnumber(res[i + 1]) - resnumber(res[i]) > 1 && chainid(res[i]) == chainid(res[i + 1])
            append!(seq, [AA_Gap for _ in 1:(resnumber(res[i + 1]) - resnumber(res[i]) - 1)])
        end
    end
    return AminoAcidSequence(seq)
end

"""
    pairalign(el1, el2, residue_selectors...)

Carries out a pairwise sequence alignment between the sequences of two
structural elements.

Additional arguments are residue selector functions - only residues that return
`true` from all the functions are retained.
The keyword arguments `scoremodel` (default
`AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1)`) and `aligntype`
(default `GlobalAlignment()`) determine the properties of the alignment.
"""
function BioAlignments.pairalign(el1::StructuralElementOrList,
                            el2::StructuralElementOrList,
                            residue_selectors::Function...;
                            scoremodel::AbstractScoreModel=AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1),
                            aligntype::BioAlignments.AbstractAlignment=GlobalAlignment())
    seq1 = AminoAcidSequence(el1, residue_selectors...; gaps=false)
    seq2 = AminoAcidSequence(el2, residue_selectors...; gaps=false)
    return pairalign(aligntype, seq1, seq2, scoremodel)
end

# DataFrame constructors for interoperability
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
function DataFrames.DataFrame(ats::Vector{<:AbstractAtom},
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
                    chainid(a), resnumber(a), inscode(a), x(a), y(a), z(a),
                    occupancy(a), tempfactor(a), element(a), charge(a),
                    modelnumber(a), isdisorderedatom(a)))
    end
    return df
end

function DataFrames.DataFrame(res::Vector{<:AbstractResidue},
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

# Descriptive showing of elements on a single line

Base.show(io::IO, struc::ProteinStructure) = print(io,
    "ProteinStructure ",
    structurename(struc) != "" ? "$(structurename(struc)) " : "",
    "with $(countmodels(struc)) models, ",
    countchains(struc) != 0 ? "$(countchains(struc)) chains ($(join(chainids(struc), ","))), " : "0 chains, ",
    "$(countresidues(struc, standardselector)) residues, ",
    "$(countatoms(struc)) atoms"
)

Base.show(io::IO, mod::Model) = print(io,
    "Model $(modelnumber(mod)) with ",
    countchains(mod) != 0 ? "$(countchains(mod)) chains ($(join(chainids(mod), ","))), " : "0 chains, ",
    "$(countresidues(mod, standardselector)) residues, ",
    "$(countatoms(mod)) atoms"
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
