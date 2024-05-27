#
# Function to perform the most important selections that are used in solute-solvent
# analysis
#
export select, Select, @sel_str

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

# Function that returns true for all atoms: the default selection
all(atoms) = true

# Main function: receives the atoms vector and a julia function to select
select(selection::String, struc) = atom_filter(Select(selection), struc)

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
    (;name, getter, operators) = key
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
    (;name, by, operators) = key
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
# Keywords for PDBTools
#

keywords = [
    Keyword(Int, "serial", serial, operators),
    Keyword(Int, "resnumber", resnumber, operators),
    Keyword(Int, "resid", resid, operators),
    Keyword(Float64, "occupancy", occupancy, operators),
    Keyword(Float64, "beta",  tempfactor, operators),
    Keyword(Float64, "tempfactor", tempfactor, operators),
    Keyword(Int, "model", modelnumber, operators),
    Keyword(Int, "modelnumber", modelnumber, operators),
    Keyword(String, "name", atomname, operators),
    Keyword(String, "atomname", atomname, operators),
#    Keyword(String, "segname", segname, operators),
    Keyword(String, "resname", resname, operators),
    Keyword(String, "chain", chainid, operators),
    Keyword(String, "chainid", chainid, operators),
]

macro_keywords = [
    MacroKeyword("water", at -> resname(at) in waterresnames),
#    MacroKeyword("water", at -> resname(at) in waterresnames),
#    MacroKeyword("protein", at -> resname(at) in proteinresnames),
#    MacroKeyword("polar", ispolar),
#    MacroKeyword("nonpolar", isnonpolar),
#    MacroKeyword("basic", isbasic),
#    MacroKeyword("acidic", isacidic),
#    MacroKeyword("charged", ischarged),
#    MacroKeyword("aliphatic", isaliphatic),
#    MacroKeyword("aromatic", isaromatic),
#    MacroKeyword("hydrophobic", ishydrophobic),
#    MacroKeyword("neutral", isneutral),
#    MacroKeyword("backbone", isbackbone),
#    MacroKeyword("sidechain", issidechain),
#    MacroKeyword("all", all),
]

#functional_keywords = [FunctionalKeyword(String, "element", element, operators)]

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


#@testitem "Selections" begin
#
#    atoms = readPDB(PDBTools.TESTPDB)
#
#    @test length(select(atoms, "name CA")) == 104
#    sel = select(atoms, "index = 13")
#    @test length(sel) == 1
#    @test sel[1].index == 13
#    @test sel[1].index_pdb == 13
#
#    @test length(select(atoms, "index > 1 and index < 13")) == 11
#
#    @test length(select(atoms, "protein")) == 1463
#
#    @test length(select(atoms, "water")) == 58014
#
#    @test length(select(atoms, "resname GLY")) == 84
#
#    @test length(select(atoms, "segname PROT")) == 1463
#
#    @test length(select(atoms, "residue = 2")) == 11
#
#    @test length(select(atoms, "neutral")) == 1233
#
#    @test length(select(atoms, "charged")) == 230
#
#    @test length(select(atoms, "sidechain")) == 854
#
#    @test length(select(atoms, "acidic")) == 162
#
#    @test length(select(atoms, "basic")) == 68
#
#    @test length(select(atoms, "hydrophobic")) == 327
#
#    @test length(select(atoms, "hydrophobic")) == 327
#
#    @test length(select(atoms, "aliphatic")) == 379
#
#    @test length(select(atoms, "aromatic")) == 344
#
#    @test length(select(atoms, "polar")) == 880
#
#    @test length(select(atoms, "nonpolar")) == 583
#
#    @test maxmin(atoms, "chain A").xlength ≈ [83.083, 83.028, 82.7]
#
#    # Test editing a field
#    atoms[1].index = 0
#    @test atoms[1].index == 0
#
#    # Test residue iterator
#    someresidues = select(atoms, "residue < 15")
#    let
#        n = 0
#        m = 0.0
#        for res in eachresidue(someresidues)
#            if name(res) == "SER"
#                n += 1
#                for atom in res
#                    m += mass(atom)
#                end
#            end
#        end
#        @test n == 4
#        @test m ≈ 348.31340000000006
#    end
#
#    # Residue properties (discontinous set)
#    lessresidues = select(someresidues, "residue < 3 or residue > 12")
#    @test residue.(eachresidue(lessresidues)) == [1, 2, 13, 14]
#    @test resnum.(eachresidue(lessresidues)) == [1, 2, 13, 14]
#    @test name.(eachresidue(lessresidues)) == ["ALA", "CYS", "SER", "SER"]
#    @test resname.(eachresidue(lessresidues)) == ["ALA", "CYS", "SER", "SER"]
#    @test chain.(eachresidue(lessresidues)) == ["A", "A", "A", "A"]
#    @test model.(eachresidue(lessresidues)) == [1, 1, 1, 1]
#
#end # testitem