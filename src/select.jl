export
    applyselectors,
    applyselectors!,
    proteinresnames,
    standardselector,
    heteroselector,
    atomnameselector,
    calphaatomnames,
    calphaselector,
    cbetaatomnames,
    cbetaselector,
    backboneatomnames,
    backboneselector,
    sidechainselector,
    heavyatomselector,
    hydrogenselector,
    resnameselector,
    proteinselector,
    acidicresselector,
    aliphaticresselector,
    aromaticresselector,
    basicresselector,
    chargedresselector,
    neutralresselector,
    hydrophobicresselector,
    polarresselector,
    nonpolarresselector,
    waterresnames,
    waterselector,
    notwaterselector,
    disorderselector,
    sscodeselector,
    helixselector,
    sheetselector,
    coilselector,
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

# Data for natural protein amino acid residues, modified from PDBTools.jl
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
    custom::Bool=true
end

const amino_acid_data = Dict{String, AminoAcidResidue}(
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

"`Set` of residue names found in proteins and peptides."
const proteinresnames = Set(keys(amino_acid_data))

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

const notsidechainatomnames = Set(["N", "CA", "C", "O", "HN", "H", "HA", "HT1", "HT2", "HT3"])

"""
    sidechainselector(at)

Determines if an `AbstractAtom` is not a hetero-atom and corresponds to a
protein side chain atom.
"""
function sidechainselector(at::AbstractAtom)
    return standardselector(at) && !atomnameselector(at, notsidechainatomnames)
end

"""
    heavyatomselector(at)

Determines if an `AbstractAtom` corresponds to a heavy (non-hydrogen) atom
and is not a hetero-atom.
"""
heavyatomselector(at::AbstractAtom) = standardselector(at) && !hydrogenselector(at)

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
    resnameselector(res, res_names)
    resnameselector(at, res_names)

Determines if an `AbstractResidue` or `AbstractAtom` has its residue name
in a list of names.
"""
function resnameselector(el::Union{AbstractResidue, AbstractAtom}, res_names)
    return resname(el, strip=false) in res_names
end

"""
    proteinselector(res)
    proteinselector(at)

Determines if an `AbstractResidue` or `AbstractAtom` is part of a protein
or peptide based on the residue name.
"""
function proteinselector(el::Union{AbstractResidue, AbstractAtom})
    return resnameselector(el, proteinresnames)
end

"""
    acidicresselector(res)
    acidicresselector(at)

Determines if an `AbstractResidue` is, or an `AbstractAtom` is part of, an
acidic amino acid based on the residue name.
"""
function acidicresselector(el::Union{AbstractResidue, AbstractAtom})
    return proteinselector(el) && amino_acid_data[resname(el)].type == "Acidic"
end

"""
    aliphaticresselector(res)
    aliphaticresselector(at)

Determines if an `AbstractResidue` is, or an `AbstractAtom` is part of, an
aliphatic amino acid based on the residue name.
"""
function aliphaticresselector(el::Union{AbstractResidue, AbstractAtom})
    return proteinselector(el) && amino_acid_data[resname(el)].type == "Aliphatic"
end

"""
    aromaticresselector(res)
    aromaticresselector(at)

Determines if an `AbstractResidue` is, or an `AbstractAtom` is part of, an
aromatic amino acid based on the residue name.
"""
function aromaticresselector(el::Union{AbstractResidue, AbstractAtom})
    return proteinselector(el) && amino_acid_data[resname(el)].type == "Aromatic"
end

"""
    basicresselector(res)
    basicresselector(at)

Determines if an `AbstractResidue` is, or an `AbstractAtom` is part of, a
basic amino acid based on the residue name.
"""
function basicresselector(el::Union{AbstractResidue, AbstractAtom})
    return proteinselector(el) && amino_acid_data[resname(el)].type == "Basic"
end

"""
    chargedresselector(res)
    chargedresselector(at)

Determines if an `AbstractResidue` is, or an `AbstractAtom` is part of, a
charged amino acid based on the residue name.
"""
function chargedresselector(el::Union{AbstractResidue, AbstractAtom})
    return proteinselector(el) && amino_acid_data[resname(el)].charge != 0
end

"""
    neutralresselector(res)
    neutralresselector(at)

Determines if an `AbstractResidue` is, or an `AbstractAtom` is part of, a
neutral amino acid based on the residue name.
"""
function neutralresselector(el::Union{AbstractResidue, AbstractAtom})
    return proteinselector(el) && amino_acid_data[resname(el)].charge == 0
end

"""
    hydrophobicresselector(res)
    hydrophobicresselector(at)

Determines if an `AbstractResidue` is, or an `AbstractAtom` is part of, a
hydrophobic amino acid based on the residue name.
"""
function hydrophobicresselector(el::Union{AbstractResidue, AbstractAtom})
    return proteinselector(el) && amino_acid_data[resname(el)].hydrophobic
end

"""
    polarresselector(res)
    polarresselector(at)

Determines if an `AbstractResidue` is, or an `AbstractAtom` is part of, a
polar amino acid based on the residue name.
"""
function polarresselector(el::Union{AbstractResidue, AbstractAtom})
    return proteinselector(el) && amino_acid_data[resname(el)].polar
end

"""
    nonpolarresselector(res)
    nonpolarresselector(at)

Determines if an `AbstractResidue` is, or an `AbstractAtom` is part of, a
non-polar amino acid based on the residue name.
"""
function nonpolarresselector(el::Union{AbstractResidue, AbstractAtom})
    return proteinselector(el) && !amino_acid_data[resname(el)].polar
end

"`Set` of residue names corresponding to water."
const waterresnames = Set(["HOH", "OH2", "TIP", "TIP3", "TIP3P", "TIP4P",
                           "TIP5P", "TIP7P", "SPC", "SPCE", "WAT"])

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

"""
    sscodeselector(res, ss_codes)
    sscodeselector(at, ss_codes)

Determines if an `AbstractResidue` or `AbstractAtom` has its secondary structure code
in a list of secondary structure codes.
"""
function sscodeselector(el::Union{AbstractResidue, AbstractAtom}, ss_codes)
    return sscode(el) in ss_codes
end

"""
    helixselector(res)
    helixselector(at)

Determines if an `AbstractResidue` or `AbstractAtom` is part of an α-helix,
i.e. whether the secondary structure code is in `helixsscodes`.
"""
function helixselector(el::Union{AbstractResidue, AbstractAtom})
    return sscodeselector(el, helixsscodes)
end

"""
    sheetselector(res)
    sheetselector(at)

Determines if an `AbstractResidue` or `AbstractAtom` is part of a β-sheet,
i.e. whether the secondary structure code is in `sheetsscodes`.
"""
function sheetselector(el::Union{AbstractResidue, AbstractAtom})
    return sscodeselector(el, sheetsscodes)
end

"""
    coilselector(res)
    coilselector(at)

Determines if an `AbstractResidue` or `AbstractAtom` is part of a coil,
i.e. whether the secondary structure code is in `coilsscodes`.
"""
function coilselector(el::Union{AbstractResidue, AbstractAtom})
    return sscodeselector(el, coilsscodes)
end

"""
    allselector(at)
    allselector(res)

Trivial selector that returns `true` for any structural element.
"""
allselector(el) = true

# Acts as a function when used within typical julia filtering functions 
#   by converting a string selection into a query call
struct Select <: Function
    sel::String
end

(s::Select)(at) = apply_query(parse_query(s.sel), at)

"String selection syntax."
macro sel_str(str)
    Select(str)
end

const operators = (
    "="  => (x, y) -> isequal(x, y),
    "<"  => (x, y) -> isless(x, y),
    ">"  => (x, y) -> isless(y, x),
    "<=" => (x, y) -> (!isless(y, x)),
    ">=" => (x, y) -> (!isless(x, y)),
)

 struct Keyword{F <: Function}
    value_type::Type
    name::String
    getter::F
    operators::Tuple
end

# For functions without parameters
struct MacroKeyword{F <: Function}
    name::String
    getter::F
end

(key::MacroKeyword)(::AbstractVector{<:AbstractString}) = key.getter

function (key::Keyword)(s::AbstractVector{<:AbstractString})
    getter, operators = key.getter, key.operators
    for op in operators
        if (i = findfirst(==(op.first), s)) !== nothing
            return el -> op.second(getter(el), parse_to_type(key, s[i+1]))
        end
    end
    # If no operator was found, assume that `=` was intended
    return el -> isequal(getter(el), parse_to_type(key, s[1]))
end

function parse_to_type(key::Keyword, val)
    if key.value_type == String
        return val
    elseif key.value_type == Char && length(val) == 1
        return val[1]
    end
    try
        val = parse(key.value_type, val)
        return val
    catch
        throw(ArgumentError("Could not parse $val for keyword $(key.name), expected $(key.value_type)"))
    end
end

const keywords = [
    Keyword(Int    , "index"      , serial     , operators),
    Keyword(Int    , "serial"     , serial     , operators),
    Keyword(Int    , "resnumber"  , resnumber  , operators),
    Keyword(Int    , "resnum"     , resnumber  , operators),
    Keyword(Int    , "resid"      , resid      , operators),
    Keyword(Float64, "occupancy"  , occupancy  , operators),
    Keyword(Float64, "beta"       , tempfactor , operators),
    Keyword(Float64, "tempfactor" , tempfactor , operators),
    Keyword(Int    , "model"      , modelnumber, operators),
    Keyword(Int    , "modelnumber", modelnumber, operators),
    Keyword(String , "name"       , atomname   , operators),
    Keyword(String , "atomname"   , atomname   , operators),
    Keyword(String , "resname"    , resname    , operators),
    Keyword(String , "chain"      , chainid    , operators),
    Keyword(String , "chainid"    , chainid    , operators),
    Keyword(String , "element"    , element    , operators),
    Keyword(Char   , "inscode"    , inscode    , operators),
    Keyword(Char   , "sscode"     , sscode     , operators),
    Keyword(Float64, "x"          , x          , operators),
    Keyword(Float64, "y"          , y          , operators),
    Keyword(Float64, "z"          , z          , operators),
    MacroKeyword("standard"   , standardselector),
    MacroKeyword("hetero"     , heteroselector),
    MacroKeyword("backbone"   , backboneselector),
    MacroKeyword("sidechain"  , sidechainselector),
    MacroKeyword("heavyatom"  , heavyatomselector),
    MacroKeyword("hydrogen"   , hydrogenselector),
    MacroKeyword("protein"    , proteinselector),
    MacroKeyword("acidic"     , acidicresselector),
    MacroKeyword("aliphatic"  , aliphaticresselector),
    MacroKeyword("aromatic"   , aromaticresselector),
    MacroKeyword("basic"      , basicresselector),
    MacroKeyword("charged"    , chargedresselector),
    MacroKeyword("neutral"    , neutralresselector),
    MacroKeyword("hydrophobic", hydrophobicresselector),
    MacroKeyword("polar"      , polarresselector),
    MacroKeyword("nonpolar"   , nonpolarresselector),
    MacroKeyword("water"      , waterselector),
    MacroKeyword("disordered" , disorderselector),
    MacroKeyword("helix"      , helixselector),
    MacroKeyword("sheet"      , sheetselector),
    MacroKeyword("coil"       , coilselector),
    MacroKeyword("all"        , allselector),
]

# See https://discourse.julialang.org/t/parsing-selection-syntax/43632/9
parse_query(selection::AbstractString) = parse_query_vector(split(selection))

function parse_query_vector(s::AbstractVector{<:AbstractString})
    # or, and, not
    if (i = findfirst(==("or"), s)) !== nothing
        deleteat!(s, i)
        return (|, parse_query_vector.((s[1:i-1], s[i:end]))...)
    elseif (i = findfirst(==("and"), s)) !== nothing
        deleteat!(s, i)
        return (&, parse_query_vector.((s[1:i-1], s[i:end]))...)
    elseif (i = findfirst(==("not"), s)) !== nothing
        deleteat!(s, i)
        return (!, parse_query_vector(s[i:end]))
    # Keywords 
    else
        for key in keywords
            if (i = findfirst(==(key.name), s)) !== nothing
                deleteat!(s, i)
                return key(s)
            end
        end
        throw(ArgumentError(("Unable to parse selection string: $s")))
    end
end

function apply_query(q, a)
    if !(q isa Tuple)
        return q(a)
    else
        f, args = Iterators.peel(q)
        return f(apply_query.(args, Ref(a))...)
    end
end
