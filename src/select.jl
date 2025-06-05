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
# by converting a string selection into a query call
struct Select{Q} <: Function
    query_string::String
    query::Q
end

function Select(query_string::AbstractString) 
    query = parse_query(query_string)
    return Select(query_string, query)
end

(s::Select)(at) = apply_query(s.query, at)

Base.show(io::IO, ::MIME"text/plain", s::Select) = print(io, """Select("$(s.query_string)")""")

# Parse selection string allowing interpolation in sel macro:
# https://discourse.julialang.org/t/str-string-interpolation/125766/11?u=lmiq
_select(args...) = Select(string(args...))

"String selection syntax."
macro sel_str(s)
    ex = Expr(:call, GlobalRef(BioStructures, :_select))
    i = firstindex(s)
    buf = IOBuffer(maxsize=ncodeunits(s))
    while i <= ncodeunits(s)
        c = @inbounds s[i]
        i = nextind(s, i)
        if c === '$'
            position(buf) > 0 && push!(ex.args, String(take!(buf)))
            val, i = Meta.parse(s, i; greedy=false)
            Meta.isexpr(val, :incomplete) && error(val.args[1])
            val !== nothing && push!(ex.args, val)
        else
            print(buf, c)
        end
    end
    position(buf) > 0 && push!(ex.args, String(take!(buf)))
    return esc(ex)
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
]

const macro_keywords = [
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

#=
    parse_query(selection:String)

Calls `parse_query_vector` after splitting the selection string.

=#
function parse_query(selection::String) 
    s = replace(selection, "(" => " ( ", ")" => " ) ")
    return parse_query_vector(split(s))
end

function apply_query(q, a)
    if !(q isa Tuple)
        q(a)
    else
        f, args = Iterators.peel(q)
        f(apply_query.(args, Ref(a))...)
    end
end

parse_error(str) = throw(ArgumentError(str))

#
# Obs: the following code were generated by Gemini 2.5-Pro, with modifications, 
# and then tested. 
#

# New helper functions
function is_operator(token::AbstractString)
    return token == "and" || token == "or" || token == "not"
end

function is_fully_enclosed(tokens::AbstractVector{<:AbstractString})
    level = 0
    # Check if the first '(' matches the last ')' without level becoming zero in between
    # for any token except the last one.
    for i in firstindex(tokens):(lastindex(tokens)-1)
        if tokens[i] == "("
            level += 1
        elseif tokens[i] == ")"
            level -= 1
            if level == 0 # Closed too early, means not fully enclosed by the outermost pair
                 return false
            end
        end
    end
    # After iterating up to tokens[end-1], level should be 1 if tokens[begin] was '('
    # and it correctly matches tokens[end]. If level is not 1, it means mismatched parentheses within.
    return level == 1
end

function find_operator_at_level_zero(op_str::String, tokens::AbstractVector{<:AbstractString})
    level = 0
    # Find first occurrence from left to right (maintaining current style)
    for i in eachindex(tokens)
        if tokens[i] == "("
            level += 1
        elseif tokens[i] == ")"
            level -= 1
            if level < 0
                parse_error("Mismatched parentheses: too many closing parentheses.")
            end
        elseif tokens[i] == op_str && level == 0
            return i
        end
    end
    if level != 0
        parse_error("Mismatched parentheses: not enough closing parentheses.")
    end
    return 0 # Not found at level zero
end

# Modified parse_query_vector
function parse_query_vector(s_vec_const::AbstractVector{<:AbstractString})
    s_vec = s_vec_const # Operate on slices or copies, not modifying original array passed around

    if isempty(s_vec)
        parse_error("Empty query segment.")
    end

    # Handle expressions fully enclosed in matching parentheses
    # e.g. "(A and B)" should be parsed by parsing "A and B"
    temp_s_vec = s_vec # Use a temporary variable for iterative stripping
    while length(temp_s_vec) > 1 && temp_s_vec[begin] == "(" && temp_s_vec[end] == ")" && is_fully_enclosed(temp_s_vec)
        temp_s_vec = temp_s_vec[begin+1:end-1]
        if isempty(temp_s_vec)
            parse_error("Empty parentheses in query: '()'")
        end
    end
    s_vec = temp_s_vec # Assign the stripped version back

    # Operator precedence: OR, then AND, then NOT (as in original code for splitting)
    # Find 'or' not within parentheses
    if (i = find_operator_at_level_zero("or", s_vec)) > 0
        left_tokens = s_vec[begin:i-1]
        right_tokens = s_vec[i+1:end]
        if isempty(left_tokens) || isempty(right_tokens)
            parse_error("Syntax error near 'or'. Missing operand.")
        end
        return (|, parse_query_vector(left_tokens), parse_query_vector(right_tokens))

    elseif (i = find_operator_at_level_zero("and", s_vec)) > 0
        left_tokens = s_vec[begin:i-1]
        right_tokens = s_vec[i+1:end]
        if isempty(left_tokens) || isempty(right_tokens)
            parse_error("Syntax error near 'and'. Missing operand.")
        end
        return (&, parse_query_vector(left_tokens), parse_query_vector(right_tokens))

    elseif s_vec[begin] == "not"
        if length(s_vec) == 1
            parse_error("Syntax error near 'not'. Missing operand.")
        end
        remaining_tokens = s_vec[begin+1:end]
        if isempty(remaining_tokens) # Should be caught by length check, but defensive
            parse_error("Syntax error near 'not'. Missing operand.")
        end
        # Prevent "not and", "not or", "not not" if "not" is not a general prefix operator in this DSL
        if is_operator(remaining_tokens[begin]) && remaining_tokens[begin] != "not" # allow "not not" if desired, though unusual
             parse_error("Operator '$(remaining_tokens[begin])' cannot directly follow 'not'.")
        end
        return (!, parse_query_vector(remaining_tokens))

    # Base case: No top-level logical operators. Must be a keyword phrase.
    else
        #if isempty(s_vec) # Should not happen if initial checks are correct
        #    parse_error("Unexpected empty query segment.")
        #end
        token_keyword_name = s_vec[begin]

        # Standard Keywords (e.g., "name", "resnum", "index")
        for key_obj in keywords # key_obj is of type Keyword
            if token_keyword_name == key_obj.name
                if length(s_vec) == 1 # Keyword name token only, no arguments
                    parse_error("Keyword '$(key_obj.name)' requires at least one argument.")
                end
                
                keyword_args = s_vec[begin+1:end] # Arguments following the keyword name
        
                is_operator_syntax_match = false
                if !isempty(keyword_args)
                    first_arg = keyword_args[1]
                    for op_tuple in key_obj.operators # e.g., ("<", isless)
                        operator_string = op_tuple[1]
                        if first_arg == operator_string
                            # Expected form: "keyword operator value", so keyword_args should be ["operator", "value"] (length 2)
                            if length(keyword_args) == 2
                                is_operator_syntax_match = true
                            else
                                parse_error(
                                    "Malformed operator expression for keyword '$(key_obj.name)'. "*
                                    "Expected 'keyword $operator_string value'. Got: $(join(s_vec, " "))"
                                )
                            end
                            break # Operator string found and processed
                        end
                    end
                end
        
                if is_operator_syntax_match
                    # Case: "keyword operator value", e.g., "resnum < 13"
                    # keyword_args will be ["<", "13"]. The Keyword functor handles this structure.
                    return key_obj(keyword_args)
                else
                    # Case: Not a recognized "keyword operator value" structure.
                    # This implies implicit equality: "keyword value" or "keyword value1 value2 ..." (for OR expansion).
        
                    #if isempty(keyword_args) # Should have been caught by length(s_vec) == 1
                    #     parse_error("No arguments provided for keyword '$(key_obj.name)'.") # Should be unreachable
                    #end
        
                    # Sanity check for multi-value: ensure no operators are present in the value list.
                    # E.g. "resnum 10 < 20" is an error here because "10" is not an operator,
                    # but "<" appears later in a context expecting only values.
                    for arg_val in keyword_args
                        for op_tuple in key_obj.operators
                            if arg_val == op_tuple[1] # op_tuple[1] is the operator string
                                parse_error(
                                    "Syntax error for keyword '$(key_obj.name)'. Operator '$(op_tuple[1])' found in an unexpected position. "*
                                    "Arguments: $(join(keyword_args, " ")). Operator expressions must be 'keyword $(op_tuple[1]) value'."
                                )
                            end
                        end
                    end
        
                    # Proceed with implicit equality (single value or multi-value OR).
                    if length(keyword_args) == 1
                        # e.g., "name CA" -> keyword_args = ["CA"]
                        # The Keyword functor handles this as implicit equality.
                        return key_obj(keyword_args)
                    else
                        # Multi-value implicit OR case, e.g., "resname ARG GLU ASP"
                        # keyword_args = ["ARG", "GLU", "ASP"]
                        current_expr_tree = key_obj([keyword_args[end]]) # Process the last value
                        for k_idx in (length(keyword_args)-1):-1:firstindex(keyword_args) # Iterate remaining values
                            current_expr_tree = (|, key_obj([keyword_args[k_idx]]), current_expr_tree)
                        end
                        return current_expr_tree
                    end
                end
            end
        end
        
        # Macro Keywords (e.g., "protein", "water")
        for key_obj in macro_keywords
            if token_keyword_name == key_obj.name
                if length(s_vec) > 1
                    parse_error("Macro keyword '$(key_obj.name)' does not take arguments. Unexpected tokens: $(join(s_vec[begin+1:end], " "))")
                end
                # MacroKeyword functor expects an argument list (empty for macros)
                return key_obj(String[]) 
            end
        end
        
        parse_error("Unknown keyword or invalid syntax at: '$(join(s_vec, " "))'")
    end
end