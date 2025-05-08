export
    PDBParseError,
    read,
    spaceatomname,
    pdbline,
    writepdb

"Error arising from parsing a Protein Data Bank (PDB) file."
struct PDBParseError <: Exception
    message::String
    line_number::Int
    line::String
end

function Base.showerror(io::IO, e::PDBParseError)
    return print(io,
            "PDBParseError: ",
            e.message,
            " at line ",
            e.line_number,
            " of file:\n",
            e.line)
end

"""
    read(filepath::AbstractString, format::Type; <keyword arguments>)
    read(input::IO, format::Type; <keyword arguments>)

Read a Protein Data Bank (PDB) file and return a `MolecularStructure`.

# Arguments
- `format::Type`: the format of the PDB file; options are PDBFormat, MMCIFFormat
    and MMTFFormat. MMTFFormat requires the MMTF.jl package to be imported.
- `structure_name::AbstractString`: the name given to the returned
    `MolecularStructure`; defaults to the file name.
- `remove_disorder::Bool=false`: whether to remove atoms with alt loc ID not ' '
    or 'A'.
- `read_std_atoms::Bool=true`: whether to read standard ATOM records.
- `read_het_atoms::Bool=true`: whether to read HETATOM records.
- `run_dssp::Bool=false`: whether to run DSSP to assign secondary structure.
    Requires the DSSP_jll.jl package to be imported if set to `true`.
- `run_stride::Bool=false`: whether to run STRIDE to assign secondary structure.
    Requires the STRIDE_jll.jl package to be imported if set to `true`.
- `gzip::Bool=false`: whether the input is gzipped, not available for PDB
    format.
"""
function Base.read(input::IO,
            ::Type{PDBFormat};
            structure_name::AbstractString="",
            remove_disorder::Bool=false,
            read_std_atoms::Bool=true,
            read_het_atoms::Bool=true,
            run_dssp::Bool=false,
            run_stride::Bool=false)
    # Define MolecularStructure and add to it incrementally
    struc = MolecularStructure(structure_name)
    struc[1] = Model(1, struc)
    # Entries outside of a MODEL/ENDMDL block are added to model 1
    curr_model = 1
    line_n = 0
    for line in eachline(input)
        line_n += 1
        # Read ATOM and HETATM records as required
        if (read_std_atoms && startswith(line, "ATOM  ")) ||
                (read_het_atoms && startswith(line, "HETATM"))
            unsafe_addatomtomodel!(
                    struc[curr_model],
                    AtomRecord(line, line_n), remove_disorder=remove_disorder)
        # Read MODEL record
        elseif startswith(line, "MODEL ")
            try
                curr_model = parse(Int, line[11:min(14, end)])
            catch
                throw(PDBParseError(
                    "could not read model serial number", line_n, line))
            end
            # Create model if required
            if !haskey(models(struc), curr_model)
                struc[curr_model] = Model(curr_model, struc)
            end
        # Read ENDMDL record
        elseif startswith(line, "ENDMDL")
            curr_model = 1
        end
    end
    # Remove any models that were not added to
    for mo in struc
        if countchains(mo) == 0
            delete!(models(struc), modelnumber(mo))
        end
    end
    # Generate lists for iteration
    fixlists!(struc)

    if run_dssp && run_stride
        throw(ArgumentError("run_dssp and run_stride cannot both be true"))
    end
    if run_dssp
        rundssp!(struc)
    end
    if run_stride
        runstride!(struc)
    end

    return struc
end

function Base.read(filepath::AbstractString,
            t::Type{<:Union{PDBFormat, MMCIFFormat, MMTFFormat}};
            structure_name::AbstractString=splitdir(filepath)[2],
            kwargs...)
    open(filepath) do input
        read(input, t; structure_name=structure_name, kwargs...)
    end
end

# Constructor from PDB ATOM/HETATM line
function AtomRecord(pdb_line::String, line_n::Integer=1)
    n = length(pdb_line)
    n >= 54 || throw(PDBParseError("line too short", line_n, pdb_line))
    AtomRecord(
        pdb_line[1] == 'H', # This assumes the line has already been checked as an ATOM/HETATM record
        parseserial(pdb_line, line_n),
        parseatomname(pdb_line, line_n),
        parsealtloc(pdb_line, line_n),
        parseresname(pdb_line, line_n),
        parsechainid(pdb_line, line_n),
        parseresnumber(pdb_line, line_n),
        parseinscode(pdb_line, line_n),
        SVector{3,Float64}((
            parsecoordx(pdb_line, line_n),
            parsecoordy(pdb_line, line_n),
            parsecoordz(pdb_line, line_n)
        )),
        n >= 60 ? parseoccupancy(pdb_line) : 1.0,
        n >= 66 ? parsetempfac(pdb_line) : 0.0,
        n >= 78 ? parseelement(pdb_line) : "  ",
        n >= 80 ? parsecharge(pdb_line) : "  ",
    )
end

function parseserial(line::String, line_n::Integer=1)
    ret = tryparse(Int, line[7:11])
    if ret === nothing
        throw(PDBParseError("could not read atom serial number", line_n, line))
    end
    return ret
end

function parseatomname(line::String, line_n::Integer=1)
    return line[13:16]
end

function parsealtloc(line::String, line_n::Integer=1)
    return line[17]
end

function parseresname(line::String, line_n::Integer=1)
    return line[18:20]
end

function parsechainid(line::String, line_n::Integer=1)
    return string(line[22])
end

function parseresnumber(line::String, line_n::Integer=1)
    ret = tryparse(Int, line[23:26])
    if ret === nothing
        throw(PDBParseError("could not read residue number", line_n, line))
    end
    return ret
end

function parseinscode(line::String, line_n::Integer=1)
    return line[27]
end

function parsecoordx(line::String, line_n::Integer=1)
    ret = tryparse(Float64, line[31:38])
    if ret === nothing
        throw(PDBParseError("could not read x coordinate", line_n, line))
    end
    return ret
end

function parsecoordy(line::String, line_n::Integer=1)
    ret = tryparse(Float64, line[39:46])
    if ret === nothing
        throw(PDBParseError("could not read y coordinate", line_n, line))
    end
    return ret
end

function parsecoordz(line::String, line_n::Integer=1)
    ret = tryparse(Float64, line[47:54])
    if ret === nothing
        throw(PDBParseError("could not read z coordinate", line_n, line))
    end
    return ret
end

function parseoccupancy(line::String)
    ret = tryparse(Float64, line[55:60])
    return (ret === nothing ? 1.0 : ret)
end

function parsetempfac(line::String)
    ret = tryparse(Float64, line[61:66])
    return (ret === nothing ? 0.0 : ret)
end

function parseelement(line::String)
    return line[77:78]
end

function parsecharge(line::String)
    return line[79:80]
end

# Form a string of a certain length from a value by adding spaces to the left
# Throws an error if the value is too long
function spacestring(val_in, new_length::Integer)
    string_out = string(val_in)
    if length(string_out) > new_length
        throw(ArgumentError("Cannot fit value \"$string_out\" into $new_length space(s)"))
    end
    return lpad(string_out, new_length)
end

"""
    spaceatomname(at::Atom)

Space an `Atom` name such that the last element letter (generally) appears in
the second column.

If the `element` property of the `Atom` is set it is used to get the element,
otherwise the name starts from the second column where possible.
This function is generally not required as spacing is recorded when atom names
are read in from a Protein Data Bank (PDB) file.
However this spacing can be important, for example distinguising between Cα and
calcium atoms.
"""
function spaceatomname(at::Atom)
    at_name = atomname(at, strip=false)
    chars = length(at_name)
    if chars == 4
        return at_name
    end
    strip_el = element(at, strip=true)
    if chars > 4
        throw(ArgumentError("Atom name is greater than four characters: \"$at_name\""))
    end
    # In the absence of the element, the first index goes in column two
    if strip_el == "" || findfirst(isequal(strip_el[1]), at_name) == nothing
        cent_ind = 1
    # The last letter of the element goes in column two where possible
    else
        cent_ind = something(findfirst(isequal(strip_el[1]), at_name), 0) + length(strip_el) - 1
    end
    if cent_ind > 2
        throw(ArgumentError("Atom name is too long to space correctly: \"$at_name\""))
    end
    if cent_ind == 1 && chars < 4
        out_string = " $at_name"
    else
        out_string = "$at_name"
    end
    return rpad(out_string, 4)
end

# Decimal places to format output to
const coordspec = FormatSpec(".3f")
const floatspec = FormatSpec(".2f")

"""
    pdbline(at::Atom)
    pdbline(at::DisorderedAtom)
    pdbline(at::AtomRecord)

Form a Protein Data Bank (PDB) format ATOM or HETATM record as a `String` from
an `Atom`, `DisorderedAtom` or `AtomRecord`.

This will throw an `ArgumentError` if a value cannot fit into the allocated
space, e.g. the chain ID is longer than one character or the atom serial is
greater than 99999.
In this case consider using `writemmcif` or `writemmtf` to write a mmCIF file or
a MMTF file.
"""
function pdbline(at::Atom)
    return (ishetero(at) ? "HETATM" : "ATOM  ") *
            # This will throw an error for serial numbers over 99999
            spacestring(serial(at), 5) *
            " " *
            spaceatomname(at) *
            string(altlocid(at)) *
            spacestring(resname(at, strip=false), 3) *
            " " *
            chainid(at) * # This is checked during writing for being one character
            spacestring(resnumber(at), 4) *
            string(inscode(at)) *
            "   " *
            # This will throw an error for large coordinate values, e.g. -1000.123
            spacestring(pyfmt(coordspec, round(x(at), digits=3)), 8) *
            spacestring(pyfmt(coordspec, round(y(at), digits=3)), 8) *
            spacestring(pyfmt(coordspec, round(z(at), digits=3)), 8) *
            spacestring(pyfmt(floatspec, round(occupancy(at), digits=2)), 6) *
            # This will throw an error for large temp facs, e.g. 1000.12
            spacestring(pyfmt(floatspec, round(tempfactor(at), digits=2)), 6) *
            "          " *
            spacestring(element(at, strip=false), 2) *
            spacestring(charge(at, strip=false), 2)
end

pdbline(dis_at::DisorderedAtom) = pdbline(defaultatom(dis_at))

function pdbline(at_rec::AtomRecord)
    return (at_rec.het_atom ? "HETATM" : "ATOM  ") *
            # This will throw an error for serial numbers over 99999
            spacestring(at_rec.serial, 5) *
            " " *
            spacestring(at_rec.atom_name, 4) *
            string(at_rec.alt_loc_id) *
            spacestring(at_rec.res_name, 3) *
            " " *
            at_rec.chain_id * # This is checked during writing for being one character
            spacestring(at_rec.res_number, 4) *
            string(at_rec.ins_code) *
            "   " *
            # This will throw an error for large coordinate values, e.g. -1000.123
            spacestring(pyfmt(coordspec, round(at_rec.coords[1], digits=3)), 8) *
            spacestring(pyfmt(coordspec, round(at_rec.coords[2], digits=3)), 8) *
            spacestring(pyfmt(coordspec, round(at_rec.coords[3], digits=3)), 8) *
            spacestring(pyfmt(floatspec, round(at_rec.occupancy, digits=2)), 6) *
            # This will throw an error for large temp facs, e.g. 1000.12
            spacestring(pyfmt(floatspec, round(at_rec.temp_factor, digits=2)), 6) *
            "          " *
            spacestring(at_rec.element, 2) *
            spacestring(at_rec.charge, 2)
end

function checkchainerror(el::Union{Chain, AbstractResidue, AbstractAtom})
    if length(chainid(el)) != 1
        throw(ArgumentError("Cannot write valid PDB file with non-single " *
            "character chain ID \"$(chainid(el))\", consider writing a mmCIF " *
            "or MMTF file"))
    end
end

"""
    writepdb(output, element, atom_selectors...)

Write a `StructuralElementOrList` to a Protein Data Bank (PDB) format file or
output stream.

Only ATOM, HETATM, MODEL and ENDMDL records are written - there is no header and
there are no TER records.
Atom selector functions can be given as additional arguments - only atoms that
return `true` from all the functions are retained.
The keyword argument `expand_disordered` (default `true`) determines whether to
return all copies of disordered residues and atoms.
"""
function writepdb(filepath::AbstractString,
                el::StructuralElementOrList,
                atom_selectors::Function...;
                expand_disordered::Bool=true)
    open(filepath, "w") do output
        writepdb(output, el, atom_selectors...;
                    expand_disordered=expand_disordered)
    end
end

function writepdb(output::IO,
                    el::Union{MolecularStructure, Vector{Model}},
                    atom_selectors::Function...;
                    expand_disordered::Bool=true)
    # If there are multiple models, write out MODEL/ENDMDL lines
    if length(el) > 1
        for mo in sort(collect(el))
            ats = collectatoms(mo, atom_selectors...; expand_disordered=expand_disordered)
            if length(ats) == 0
                continue
            end
            println(output, "MODEL     ", spacestring(modelnumber(mo), 4), repeat(" ", 66))
            writepdb(output, ats; expand_disordered=expand_disordered)
            println(output, "ENDMDL$(repeat(" ", 74))")
        end
    # If there is only one model, do not write out MODEL/ENDMDL lines
    elseif length(el) == 1
        writepdb(output, first(el), atom_selectors...;
                    expand_disordered=expand_disordered)
    end
end

function writepdb(output::IO,
                el::Union{Model, Chain, AbstractResidue, AbstractAtom,
                            Vector{Chain}, Vector{<:AbstractResidue},
                            Vector{<:AbstractAtom}},
                atom_selectors::Function...;
                expand_disordered::Bool=true)
    for at in collectatoms(el, atom_selectors...;
                            expand_disordered=expand_disordered)
        checkchainerror(at)
        println(output, pdbline(at))
    end
end
