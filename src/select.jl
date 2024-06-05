export
    applyselectors,
    applyselectors!,
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
    @sel_str

"""
    applyselectors(els, selectors...)

Returns a copy of a `Vector` of `StructuralElement`s with all elements that do
not return `true` from all the selector functions removed.
"""
function applyselectors(els::AbstractVector{<:StructuralElement}, selectors::Function...)
    new_list = copy(els)
    applyselectors!(new_list, selectors...)
    return new_list
end

"""
    applyselectors!(els, selectors...)

Removes from a `Vector` of `StructuralElement`s all elements that do not return
`true` from all the selector functions.
"""
function applyselectors!(els::AbstractVector{<:StructuralElement}, selectors::Function...)
    for selector in selectors
        filter!(selector, els)
    end
    return els
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

Determines if an `AbstractAtom` has its atom name in a list of names.
`strip` determines whether surrounding whitespace is stripped from the atom
name before it is checked in the list.
"""
function atomnameselector(at::AbstractAtom, atom_names; strip::Bool=true)
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
in a list of names.
"""
function resnameselector(el::Union{AbstractResidue, AbstractAtom}, res_names)
    return resname(el, strip=false) in res_names
end

"`Set` of residue names corresponding to water."
const waterresnames = Set(["HOH", "WAT"])

"""
    waterselector(res)
    waterselector(at)

Determines if an `AbstractResidue` or `AbstractAtom` represents a water
molecule, i.e. whether the residue name is in `waterresnames`.
"""
function waterselector(el::Union{AbstractResidue, AbstractAtom})
    return resnameselector(el, waterresnames)
end

"""
    notwaterselector(res)
    notwaterselector(at)

Determines if an `AbstractResidue` or `AbstractAtom` does not represent a
water molecule, i.e. whether the residue name is not in `waterresnames`.
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

#
# Dictionaries of residue, atom, and element types, as originally written
# in PDBTools.jl (at first sight, seems more general that the current)
#
#
# Data for natural protein residues
#
Base.@kwdef struct AminoAcidResidue
    name::String
    three_letter_code::String
    one_letter_code::String
    type::String
    polar::Bool
    hydrophobic::Bool
    mono_isotopic_mass::Float64
    mass::Float64
    charge::Int
    custom::Bool = true
end
#! format: off
const protein_residues = Dict{String,AminoAcidResidue}(
    "ALA" => AminoAcidResidue("Alanine",       "ALA", "A", "Aliphatic",  false, false,  71.037114,  71.0779,  0, false),
    "ARG" => AminoAcidResidue("Arginine",      "ARG", "R", "Basic",      true,  false, 156.101111, 156.1857,  1, false),
    "ASN" => AminoAcidResidue("Asparagine",    "ASN", "N", "Amide",      true,  false, 114.042927, 114.1026,  0, false),
    "ASP" => AminoAcidResidue("Aspartic acid", "ASP", "D", "Acidic",     true,  false, 115.026943, 115.0874, -1, false),
    "CYS" => AminoAcidResidue("Cysteine",      "CYS", "C", "Sulfuric",   false, false, 103.009185, 103.1429,  0, false),
    "GLN" => AminoAcidResidue("Glutamine",     "GLN", "Q", "Amide",      true,  false, 128.058578, 128.1292,  0, false),
    "GLU" => AminoAcidResidue("Glutamic acid", "GLU", "E", "Acidic",     true,  false, 129.042593, 129.1140, -1, false),
    "GLY" => AminoAcidResidue("Glycine",       "GLY", "G", "Aliphatic",  false, false,  57.021464,  57.0513,  0, false),
    "HIS" => AminoAcidResidue("Histidine",     "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "ILE" => AminoAcidResidue("Isoleucine",    "ILE", "I", "Aliphatic",  false, true,  113.084064, 113.1576,  0, false),
    "LEU" => AminoAcidResidue("Leucine",       "LEU", "L", "Aliphatic",  false, true,  113.084064, 113.1576,  0, false),
    "LYS" => AminoAcidResidue("Lysine",        "LYS", "K", "Basic",      true,  false, 128.094963, 128.1723,  1, false),
    "MET" => AminoAcidResidue("Methionine",    "MET", "M", "Sulfuric",   false, false, 131.040485, 131.1961,  0, false),
    "PHE" => AminoAcidResidue("Phenylalanine", "PHE", "F", "Aromatic",   false, true,  147.068414, 147.1739,  0, false),
    "PRO" => AminoAcidResidue("Proline",       "PRO", "P", "Cyclic",     false, false,  97.052764,  97.1152,  0, false),
    "SER" => AminoAcidResidue("Serine",        "SER", "S", "Hydroxylic", true,  false,  87.032028, 87.07730,  0, false),
    "THR" => AminoAcidResidue("Threonine",     "THR", "T", "Hydroxylic", true,  false, 101.047679, 101.1039,  0, false),
    "TRP" => AminoAcidResidue("Tryptophan",    "TRP", "W", "Aromatic",   false, true,  186.079313, 186.2099,  0, false),
    "TYR" => AminoAcidResidue("Tyrosine",      "TYR", "Y", "Aromatic",   true,  false, 163.063320, 163.1733,  0, false),
    "VAL" => AminoAcidResidue("Valine",        "VAL", "V", "Aliphatic",  false, true,   99.068414,  99.1311,  0, false),
    # Alternate protonation states for CHARMM and AMBER
    "ASPP" => AminoAcidResidue("Aspartic acid (protonated)", "ASP", "D", "Acidic", true,  false, 115.026943, 115.0874, 0, false),
    "GLUP" => AminoAcidResidue("Glutamic acid (protonated)", "GLU", "E", "Acidic", true,  false, 129.042593, 129.1140, 0, false),
    "HSD"  => AminoAcidResidue("Histidine (D)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "HSE"  => AminoAcidResidue("Histidine (E)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "HSP"  => AminoAcidResidue("Histidine (doubly protonated)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  1, false), 
    "HID"  => AminoAcidResidue("Histidine (D)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "HIE"  => AminoAcidResidue("Histidine (E)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "HIP"  => AminoAcidResidue("Histidine (doubly protonated)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  1, false), 
)
#! format: on
isprotein(atom::AbstractAtom) = haskey(protein_residues, resname(atom))
isacidic(atom::AbstractAtom) = isprotein(atom) && protein_residues[resname(atom)].type == "Acidic"
isaliphatic(atom::AbstractAtom) = isprotein(atom) && protein_residues[resname(atom)].type == "Aliphatic"
isaromatic(atom::AbstractAtom) = isprotein(atom) && protein_residues[resname(atom)].type == "Aromatic"
isbasic(atom::AbstractAtom) = isprotein(atom) && protein_residues[resname(atom)].type == "Basic"
ischarged(atom::AbstractAtom) = isprotein(atom) && protein_residues[resname(atom)].charge != 0
isneutral(atom::AbstractAtom) = isprotein(atom) && protein_residues[resname(atom)].charge == 0
ishydrophobic(atom::AbstractAtom) = isprotein(atom) && protein_residues[resname(atom)].hydrophobic
ispolar(atom::AbstractAtom) = isprotein(atom) && protein_residues[resname(atom)].polar
isnonpolar(atom::AbstractAtom) = isprotein(atom) && !ispolar(atom)
const backbone_atoms = ["N", "CA", "C", "O"]
isbackbone(atom::AbstractAtom; backbone_atoms=backbone_atoms) = isprotein(atom) && atomname(atom) in backbone_atoms
const not_side_chain_atoms = ["N", "CA", "C", "O", "HN", "H", "HA", "HT1", "HT2", "HT3"]
issidechain(atom::AbstractAtom; not_side_chain_atoms=not_side_chain_atoms) = isprotein(atom) && !(atomname(atom) in not_side_chain_atoms)
const water_residues = ["HOH", "OH2", "TIP", "TIP3", "TIP3P", "TIP4P", "TIP5P", "TIP7P", "SPC", "SPCE"]
iswater(atom::AbstractAtom) = strip(resname(atom)) in water_residues

#=
    Select

This structure acts a function when used within typical julia filtering functions, 
by converting a string selection into a call to query call. 

=#
struct Select <: Function
    sel::String
end
(s::Select)(at) = apply_query(parse_query(s.sel), at)
macro sel_str(str)
    Select(str)
end

# collect atoms function
function collectatoms(struc::StructuralElementOrList, sel::Select)
    atoms = collectatoms(struc)
    return filter!(sel, atoms)
end

# Comparison operators
const operators = (
    "=" => (x, y) -> isequal(x, y),
    "<" => (x, y) -> isless(x, y),
    ">" => (x, y) -> isless(y, x),
    "<=" => (x, y) -> (!isless(y, x)),
    ">=" => (x, y) -> (!isless(x, y)),
)

#
# Keywords
#
struct Keyword{F<:Function}
    ValueType::Type
    name::String
    getter::F
    operators::Tuple
end

#
# Macro keywords (functions without parameters)
#
struct MacroKeyword{F<:Function}
    name::String
    getter::F
end
(key::MacroKeyword)(::AbstractVector{<:AbstractString}) = key.getter

function (key::Keyword)(s::AbstractVector{<:AbstractString})
    # if 1.6 compatibility is not needed, this can be replaced by
    # (; getter, operators) = key
    getter, operators = key.getter, key.operators
    for op in operators
        if (i = findfirst(==(op.first), s)) !== nothing
            return el -> op.second(getter(el), parse_to_type(key, s[i+1]))
        end
    end
    # If no operator was found, assume that `=` was intended
    return el -> isequal(getter(el), parse_to_type(key, s[1]))
end

#=
    parse_to_type(key::Keyword, val::String)

Tries to parse `val` into the type of value expected by `key.ValueType`. 

=#
function parse_to_type(key::Keyword, val)
    if key.ValueType == String
        return val
    end
    try
        val = parse(key.ValueType, val)
        return val
    catch
        throw(ArgumentError("Could not parse $val for keyword $(key.name), expected $(key.ValueType)"))
    end
end

#
# Keywords
#
keywords = [
    Keyword(Int, "index", serial, operators),
    Keyword(Int, "serial", serial, operators),
    Keyword(Int, "resnumber", resnumber, operators),
    Keyword(Int, "resnum", resnumber, operators),
    Keyword(Int, "resid", resid, operators),
    Keyword(Float64, "occupancy", occupancy, operators),
    Keyword(Float64, "beta", tempfactor, operators),
    Keyword(Float64, "tempfactor", tempfactor, operators),
    Keyword(Int, "model", modelnumber, operators),
    Keyword(Int, "modelnumber", modelnumber, operators),
    Keyword(String, "name", atomname, operators),
    Keyword(String, "atomname", atomname, operators),
    #Keyword(String, "segname", segname, operators),
    Keyword(String, "resname", resname, operators),
    Keyword(String, "chain", chainid, operators),
    Keyword(String, "chainid", chainid, operators),
    Keyword(String, "element", element, operators),
    MacroKeyword("water", iswater),
    MacroKeyword("protein", isprotein),
    MacroKeyword("polar", ispolar),
    MacroKeyword("nonpolar", isnonpolar),
    MacroKeyword("basic", isbasic),
    MacroKeyword("acidic", isacidic),
    MacroKeyword("charged", ischarged),
    MacroKeyword("aliphatic", isaliphatic),
    MacroKeyword("aromatic", isaromatic),
    MacroKeyword("hydrophobic", ishydrophobic),
    MacroKeyword("neutral", isneutral),
    MacroKeyword("backbone", isbackbone),
    MacroKeyword("heavyatom", heavyatomselector),
    MacroKeyword("disordered", isdisorderedatom),
    MacroKeyword("sidechain", issidechain),
    MacroKeyword("all", el -> true),
]

#
# parse_query and apply_query are a very gentle contribution given by 
# CameronBieganek in https://discourse.julialang.org/t/parsing-selection-syntax/43632/9
# while explaining to me how to creat a syntex interpreter
#

#=
    parse_query(selection:String)

Calls `parse_query_vector` after splitting the selection string.

=#
parse_query(selection::String) = parse_query_vector(split(selection))

#=

    parse_query_vector(s::AbstractVector{<:AbstractString})

=#
function parse_query_vector(s)
    # or, and, not
    if (i = findfirst(==("or"), s)) !== nothing
        deleteat!(s, i)
        (|, parse_query_vector.((s[1:i-1], s[i:end]))...)
    elseif (i = findfirst(==("and"), s)) !== nothing
        deleteat!(s, i)
        (&, parse_query_vector.((s[1:i-1], s[i:end]))...)
    elseif (i = findfirst(==("not"), s)) !== nothing
        deleteat!(s, i)
        (!, parse_query_vector(s[i:end]))
    # keywords 
    else
        for key in keywords
            if (i = findfirst(==(key.name), s)) !== nothing
                deleteat!(s, i)
                return key(s)
            end
        end
        throw(ArgumentError(("Unable to parse selection string.")))
    end
end

function apply_query(q, a)
    if !(q isa Tuple)
        q(a)
    else
        f, args = Iterators.peel(q)
        f(apply_query.(args, Ref(a))...)
    end
end



