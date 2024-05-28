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

"""
    Select

This structure acts a function when used within typical julia filtering functions, 
by converting a string selection into a call to query call. 

# Example

```jldoctest
julia> using PDBTools

julia> atoms = readPDB(PDBTools.TESTPDB, "protein");

julia> findfirst(Select("name CA"), atoms)
5

julia> filter(Select("name CA and residue 1"), atoms)
   Array{Atoms,1} with 1 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       5   CA     ALA     A        1        1   -8.483  -14.912   -6.726  1.00  0.00     1    PROT         5

```
"""
struct Select <: Function
    sel::String
end
(s::Select)(at) = apply_query(parse_query(s.sel), at)

macro sel_str(str)
    Select(str)
end

function atom_filter(by::F, struc::ProteinStructure) where {F<:Function}
    sel = AbstractAtom[]
    for mod in struc
        for ch in mod
            for res in ch
                for at in res
                    if by(at)
                        push!(sel, at)
                    end
                end
            end
        end
    end
    return sel
end

function atom_filter(by::F, atoms::AbstractVector{<:AbstractAtom}) where {F<:Function}
    sel = AbstractAtom[]
    for at in atoms
        if by(at)
            push!(sel, at)
        end
    end
    return sel
end

# Overload Base functions to forward the selection syntax
Base.filter(by::Select, struc::ProteinStructure) = atom_filter(by, struc)
Base.findall(by::Select, struc::ProteinStructure) = serial.(atom_filter(by, struc))

# collect atoms function 
collectatoms(atoms::AbstractVector{<:AbstractAtom}, sel::Select) = atom_filter(sel, atoms)

# Function that returns true for all atoms: the default selection
all(atoms) = true

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

function (key::Keyword)(s::AbstractVector{<:AbstractString})
    # if 1.6 compatibility is not needed, this can be replaced by
    # (; name, getter, operators) = key
    name, getter, operators = key.name, key.getter, key.operators
    for op in operators
        if (i = has_key(op.first, s)) > 0
            return el -> op.second(getter(el), parse_to_type(key, s[i+1]))
        end
    end
    # If no operator was found, assume that `=` was intended
    i = has_key(name, s)
    return el -> isequal(getter(el), parse_to_type(key, s[i+1]))
end

#=
    FunctionalKeyword{T}

This is a structure that will store a keyword that depends on an external function
requiring an operator and an argument. 

## Example:

```
element_keyword = FunctionalKeyword(String,      
                                    "element",
                                    element,
                                    ("=",isequal))

```
will define a keyword "element" to be used as `element C`, which will return
`true` if there is an `element` function such that `element(atom) == C`.

=#
struct FunctionalKeyword{FunctionType}
    ValueType::Type
    name::String
    by::FunctionType
    operators::Tuple
end

function (key::FunctionalKeyword)(s::AbstractVector{<:AbstractString})
    # if 1.6 compatibility is not needed, this can be replaced by
    # (; name, by, operators) = key
    name, by, operators = key.name, key.by, key.operators
    for op in operators
        if (i = has_key(op.first, s)) > 0
            return el -> op.second(by(el), parse_to_type(key, s[i+1]))
        end
    end
    # If no operator was found, assume that `=` was intended
    i = has_key(name, s)
    return el -> isequal(by(el), parse_to_type(key, s[i+1]))
end

#
# Macro keywords (functions without parameters)
#
struct MacroKeyword{FunctionType}
    name::String
    by::FunctionType
end

function (key::MacroKeyword)(s::AbstractVector{<:AbstractString})
    return key.by
end

#=
    parse_to_type(key::Keyword, val::String)

Tries to parse `val` into the type of value expected by `key.ValueType`. 

=#
function parse_to_type(key::Union{Keyword,FunctionalKeyword}, val)
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
]

macro_keywords = [
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

functional_keywords = [FunctionalKeyword(String, "element", element, operators)]

#
# parse_query and apply_query are a very gentle contribution given by 
# CameronBieganek in https://discourse.julialang.org/t/parsing-selection-syntax/43632/9
# while explaining to me how to creat a syntex interpreter
#
#=
    has_key(key::String, s::AbstractVector{<:AbstractString})

Returns the first index of the vector `s` in which where `key` is found, or 0. 

## Example:

```julia-repl

julia> PDBTools.has_key("or",["name","CA","or","index","1"])
3

julia> PDBTools.has_key("and",["name","CA","or","index","1"])
0

```

=#
function has_key(key::String, s::AbstractVector{<:AbstractString})
    i = findfirst(isequal(key), s)
    if isnothing(i)
        0
    else
        i
    end
end

#=
    parse_query(selection:String)

Calls `parse_query_vector` after splitting the selection string.

=#
parse_query(selection::String) = parse_query_vector(split(selection))

#=

    parse_query_vector(s::AbstractVector{<:AbstractString})

=#
function parse_query_vector(s)
    if (i = has_key("or", s)) > 0
        deleteat!(s, i)
        (|, parse_query_vector.((s[1:i-1], s[i:end]))...)
    elseif (i = has_key("and", s)) > 0
        deleteat!(s, i)
        (&, parse_query_vector.((s[1:i-1], s[i:end]))...)
    elseif (i = has_key("not", s)) > 0
        deleteat!(s, i)
        (!, parse_query_vector(s[i:end]))

        # keywords 
    else
        for key in keywords
            if (i = has_key(key.name, s)) > 0
                deleteat!(s, i)
                return key(s)
            end
        end
        for key in macro_keywords
            if (i = has_key(key.name, s)) > 0
                deleteat!(s, i)
                return key(s)
            end
        end
        for key in functional_keywords
            if (i = has_key(key.name, s)) > 0
                deleteat!(s, i)
                return key(s)
            end
        end
        parse_error("Unable to parse selection string.")
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

#
# Simple error message (https://discourse.julialang.org/t/a-julia-equivalent-to-rs-stop/36568/13)
#
struct NoBackTraceException
    exc::Exception
end

function Base.showerror(io::IO, ex::NoBackTraceException, bt; backtrace=true)
    Base.with_output_color(get(io, :color, false) ? Base.error_color() : :nothing, io) do io
        showerror(io, ex.exc)
    end
end
parse_error(str) = throw(NoBackTraceException(ErrorException(str)))


@testitem "Selections" begin

    using BioStructures
    atoms = read(BioStructures.TESTPDB, PDB)

    @test length(filter(sel"name CA", atoms)) == 104
    sel = filter(sel"index = 13", atoms)
    @test length(sel) == 1
    @test serial(sel[1]) == 13

    @test length(filter(sel"index > 1 and index < 13", atoms)) == 11
    @test length(filter(sel"protein", atoms)) == 1463
    @test length(filter(sel"water", atoms)) == 58014
    @test length(filter(sel"resname GLY", atoms)) == 84
    #@test length(filter(sel"segname PROT", atoms)) == 1463
    # residue should be the incremental residue number
    #@test length(filter(sel"residue = 2", atoms)) == 11
    @test length(filter(sel"protein and resnum = 2", atoms)) == 11
    @test length(filter(sel"neutral", atoms)) == 1233
    @test length(filter(sel"charged", atoms)) == 230
    @test length(filter(sel"sidechain", atoms)) == 854
    @test length(filter(sel"acidic", atoms)) == 162
    @test length(filter(sel"basic", atoms)) == 68
    @test length(filter(sel"hydrophobic", atoms)) == 327
    @test length(filter(sel"hydrophobic", atoms)) == 327
    @test length(filter(sel"aliphatic", atoms)) == 379
    @test length(filter(sel"aromatic", atoms)) == 344
    @test length(filter(sel"polar", atoms)) == 880
    @test length(filter(sel"nonpolar", atoms)) == 583

    @test collectatoms(atoms, sel"resname LYS") == filter(sel"resname LYS", atoms)

end # testitem

