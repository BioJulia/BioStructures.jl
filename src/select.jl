using TestItems # suggested to be used all around
#
# Function to perform the most important selections that are used in solute-solvent
# analysis
#
export @sel_str
const TESTPDB = joinpath(@__DIR__, "..", "test", "structure.pdb")

#
# Dictionaries of residue, atom, and element types, as originally written
# in PDBTools.jl (at first sight, seems more general that the current)
#
#
# Data for natural protein residues
#
Base.@kwdef struct ProteinResidue
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
const protein_residues = Dict{String,ProteinResidue}(
    "ALA" => ProteinResidue("Alanine",       "ALA", "A", "Aliphatic",  false, false,  71.037114,  71.0779,  0, false),
    "ARG" => ProteinResidue("Arginine",      "ARG", "R", "Basic",      true,  false, 156.101111, 156.1857,  1, false),
    "ASN" => ProteinResidue("Asparagine",    "ASN", "N", "Amide",      true,  false, 114.042927, 114.1026,  0, false),
    "ASP" => ProteinResidue("Aspartic acid", "ASP", "D", "Acidic",     true,  false, 115.026943, 115.0874, -1, false),
    "CYS" => ProteinResidue("Cysteine",      "CYS", "C", "Sulfuric",   false, false, 103.009185, 103.1429,  0, false),
    "GLN" => ProteinResidue("Glutamine",     "GLN", "Q", "Amide",      true,  false, 128.058578, 128.1292,  0, false),
    "GLU" => ProteinResidue("Glutamic acid", "GLU", "E", "Acidic",     true,  false, 129.042593, 129.1140, -1, false),
    "GLY" => ProteinResidue("Glycine",       "GLY", "G", "Aliphatic",  false, false,  57.021464,  57.0513,  0, false),
    "HIS" => ProteinResidue("Histidine",     "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "ILE" => ProteinResidue("Isoleucine",    "ILE", "I", "Aliphatic",  false, true,  113.084064, 113.1576,  0, false),
    "LEU" => ProteinResidue("Leucine",       "LEU", "L", "Aliphatic",  false, true,  113.084064, 113.1576,  0, false),
    "LYS" => ProteinResidue("Lysine",        "LYS", "K", "Basic",      true,  false, 128.094963, 128.1723,  1, false),
    "MET" => ProteinResidue("Methionine",    "MET", "M", "Sulfuric",   false, false, 131.040485, 131.1961,  0, false),
    "PHE" => ProteinResidue("Phenylalanine", "PHE", "F", "Aromatic",   false, true,  147.068414, 147.1739,  0, false),
    "PRO" => ProteinResidue("Proline",       "PRO", "P", "Cyclic",     false, false,  97.052764,  97.1152,  0, false),
    "SER" => ProteinResidue("Serine",        "SER", "S", "Hydroxylic", true,  false,  87.032028, 87.07730,  0, false),
    "THR" => ProteinResidue("Threonine",     "THR", "T", "Hydroxylic", true,  false, 101.047679, 101.1039,  0, false),
    "TRP" => ProteinResidue("Tryptophan",    "TRP", "W", "Aromatic",   false, true,  186.079313, 186.2099,  0, false),
    "TYR" => ProteinResidue("Tyrosine",      "TYR", "Y", "Aromatic",   true,  false, 163.063320, 163.1733,  0, false),
    "VAL" => ProteinResidue("Valine",        "VAL", "V", "Aliphatic",  false, true,   99.068414,  99.1311,  0, false),
    # Alternate protonation states for CHARMM and AMBER
    "ASPP" => ProteinResidue("Aspartic acid (protonated)", "ASP", "D", "Acidic", true,  false, 115.026943, 115.0874, 0, false),
    "GLUP" => ProteinResidue("Glutamic acid (protonated)", "GLU", "E", "Acidic", true,  false, 129.042593, 129.1140, 0, false),
    "HSD"  => ProteinResidue("Histidine (D)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "HSE"  => ProteinResidue("Histidine (E)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "HSP"  => ProteinResidue("Histidine (doubly protonated)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  1, false), 
    "HID"  => ProteinResidue("Histidine (D)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "HIE"  => ProteinResidue("Histidine (E)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  0, false), 
    "HIP"  => ProteinResidue("Histidine (doubly protonated)", "HIS", "H", "Aromatic",   true,  false, 137.058912, 137.1393,  1, false), 
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


# Function that returns true for all atoms: the default selection
all(atoms::AbstractVector{<:AbstractAtom}) = true

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
        parse_error(
            "Could not parse $val for keyword $(key.name), expected $(key.ValueType)",
        )
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
    MacroKeyword("all", all),
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

@testitem "Selections" begin

    using BioStructures
    struc = read(BioStructures.TESTPDB, PDB)

    @test length(collectatoms(struc, sel"name CA")) == 104
    sel = collectatoms(struc, sel"index = 13")
    @test length(sel) == 1
    @test serial(sel[1]) == 13

    @test length(collectatoms(struc, sel"index > 1 and index < 13")) == 11
    @test length(collectatoms(struc, sel"protein")) == 1463
    @test length(collectatoms(struc, sel"water")) == 135
    @test length(collectatoms(struc, sel"resname GLY")) == 84
    #@test length(filter(sel"segname PROT", struc)) == 1463
    # residue should be the incremental residue number
    #@test length(filter(sel"residue = 2", struc)) == 11
    @test length(collectatoms(struc, sel"protein and resnum = 2")) == 11
    @test length(collectatoms(struc, sel"neutral")) == 1233
    @test length(collectatoms(struc, sel"charged")) == 230
    @test length(collectatoms(struc, sel"sidechain")) == 854
    @test length(collectatoms(struc, sel"acidic")) == 162
    @test length(collectatoms(struc, sel"basic")) == 68
    @test length(collectatoms(struc, sel"hydrophobic")) == 327
    @test length(collectatoms(struc, sel"not hydrophobic")) == 1286
    @test length(collectatoms(struc, sel"aliphatic")) == 379
    @test length(collectatoms(struc, sel"aromatic")) == 344
    @test length(collectatoms(struc, sel"polar")) == 880
    @test length(collectatoms(struc, sel"nonpolar")) == 583
    @test length(collectatoms(struc, sel"backbone")) == 415
    @test length(collectatoms(struc, sel"element H")) == 538

    @test_throws ArgumentError collectatoms(struc, sel"abc")

end # testitem

