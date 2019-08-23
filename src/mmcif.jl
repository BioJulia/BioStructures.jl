export
    MMCIFDict,
    writemmcif

# mmCIF special characters
const quotechars = Set(['\'', '\"'])
const whitespacechars = Set([' ', '\t'])
const missingvals = Set([".", "?"])
const specialchars = Set(['_', '#', '\$', '[', ']', ';'])
const specialwords = Set(["loop_", "stop_", "global_"])

# If certain entries should have a certain order of keys, that is specified here
const mmciforder = Dict(
    "_atom_site"=> [
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_iso_or_equiv",
        "pdbx_formal_charge",
        "auth_seq_id",
        "auth_comp_id",
        "auth_asym_id",
        "auth_atom_id",
        "pdbx_PDB_model_num",
    ]
)


"""
A macromolecular Crystallographic Information File (mmCIF) dictionary.

Can be accessed using similar functions to a standard `Dict`.
Keys are field names as a `String` and values are always `Vector{String}`, even
for multiple components or numerical data.
To directly access the underlying dictionary of `MMCIFDict` `d`, use
`d.dict`.
Call `MMCIFDict` with a filepath or stream to read the dictionary from that
source.
"""
struct MMCIFDict
    dict::Dict{String, Vector{String}}
end

MMCIFDict() = MMCIFDict(Dict())

Base.getindex(mmcif_dict::MMCIFDict, field::AbstractString) = mmcif_dict.dict[field]

function Base.setindex!(mmcif_dict::MMCIFDict,
                    val::Vector{String},
                    field::AbstractString)
    mmcif_dict.dict[field] = val
    return mmcif_dict
end

Base.keys(mmcif_dict::MMCIFDict) = keys(mmcif_dict.dict)
Base.values(mmcif_dict::MMCIFDict) = values(mmcif_dict.dict)
Base.haskey(mmcif_dict::MMCIFDict, key) = haskey(mmcif_dict.dict, key)

function Base.show(io::IO, mmcif_dict::MMCIFDict)
    print(io, "mmCIF dictionary with $(length(keys(mmcif_dict))) fields")
end


# Split a mmCIF line into tokens
# See https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax for syntax
function splitline(s::AbstractString)
    tokens = String[]
    in_token = false
    # Quote character of the currently open quote, or ' ' if no quote open
    quote_open_char = ' '
    start_i = 0
    for (i, c) in enumerate(s)
        if c in whitespacechars
            if in_token && quote_open_char == ' '
                in_token = false
                push!(tokens, s[start_i:i-1])
            end
        elseif c in quotechars
            if quote_open_char == ' '
                if in_token
                    throw(ArgumentError("Opening quote in middle of word: $s"))
                end
                quote_open_char = c
                in_token = true
                start_i = i + 1
            elseif c == quote_open_char && (i == length(s) || s[i+1] in whitespacechars)
                quote_open_char = ' '
                in_token = false
                push!(tokens, s[start_i:i-1])
            end
        elseif c == '#' && quote_open_char == ' '
            return tokens
        elseif !in_token
            in_token = true
            start_i = i
        end
    end
    if in_token
        push!(tokens, s[start_i:end])
    end
    if quote_open_char != ' '
        throw(ArgumentError("Line ended with quote open: $s"))
    end
    return tokens
end


# Get tokens from a mmCIF file
function tokenizecif(f::IO)
    tokens = String[]
    for line in eachline(f)
        if startswith(line, "#")
            continue
        elseif startswith(line, ";")
            token_buffer = [rstrip(line[2:end])]
            for inner_line in eachline(f)
                inner_strip = rstrip(inner_line)
                if inner_strip == ";"
                    break
                end
                push!(token_buffer, inner_strip)
            end
            push!(tokens, join(token_buffer, "\n"))
        else
            append!(tokens, splitline(line))
        end
    end
    return tokens
end

# Get tokens from a mmCIF file corresponding to atom site records only
# This will fail if there is only a single atom record in the file
#   and it is not in the loop format
function tokenizecifstructure(f::IO)
    tokens = String[]
    reading = false
    in_keys = true
    for line in eachline(f)
        if (!reading && !startswith(line, "_atom_site.")) || startswith(line, "#")
            continue
        elseif startswith(line, "_atom_site.")
            reading = true
            push!(tokens, rstrip(line))
        elseif startswith(line, ";")
            in_keys = false
            token_buffer = [rstrip(line[2:end])]
            for inner_line in eachline(f)
                inner_strip = rstrip(inner_line)
                if inner_strip == ";"
                    break
                end
                push!(token_buffer, inner_strip)
            end
            push!(tokens, join(token_buffer, "\n"))
        elseif (!in_keys && startswith(line, "_")) || startswith(line, "loop_")
            break
        else
            in_keys = false
            append!(tokens, splitline(line))
        end
    end
    length(tokens) > 0 ? pushfirst!(tokens, "loop_") : nothing
    return tokens
end


function MMCIFDict(mmcif_filepath::AbstractString)
    open(mmcif_filepath) do f
        MMCIFDict(f)
    end
end

# Read a mmCIF file into a MMCIFDict
function MMCIFDict(f::IO)
    mmcif_dict = MMCIFDict()
    tokens = tokenizecif(f)
    # Data label token is read first
    if length(tokens) == 0
        return mmcif_dict
    end
    data_token = first(tokens)
    mmcif_dict[data_token[1:5]] = [data_token[6:end]]
    return populatedict!(mmcif_dict, tokens[2:end])
end


# Add tokens to a mmCIF dictionary
function populatedict!(mmcif_dict::MMCIFDict, tokens::Vector{String})
    key = ""
    keys = String[]
    loop_flag = false
    i = 0 # Value counter
    n = 0 # Key counter
    for token in tokens
        if lowercase(token) == "loop_"
            loop_flag = true
            keys = String[]
            i = 0
            n = 0
            continue
        elseif loop_flag
            # The second condition checks we are in the first column
            # Some mmCIF files (e.g. 4q9r) have values in later columns
            #   starting with an underscore and we don't want to read
            #   these as keys
            if startswith(token, "_") && (n == 0 || i % n == 0)
                if i > 0
                    loop_flag = false
                else
                    mmcif_dict[token] = String[]
                    push!(keys, token)
                    n += 1
                    continue
                end
            else
                try
                    push!(mmcif_dict[keys[i % n + 1]], token)
                catch ex
                    # A zero division error means we have not found any keys
                    if isa(ex, DivideError)
                        throw(ArgumentError("Loop keys not found, token: \"$token\""))
                    else
                        rethrow()
                    end
                end
                i += 1
                continue
            end
        end
        if key == ""
            key = token
        else
            # Single values are read in as an array for consistency
            mmcif_dict[key] = [token]
            key = ""
        end
    end
    return mmcif_dict
end


function Base.read(input::IO,
            ::Type{MMCIF};
            structure_name::AbstractString="",
            remove_disorder::Bool=false,
            read_std_atoms::Bool=true,
            read_het_atoms::Bool=true)
    mmcif_dict = MMCIFDict()
    populatedict!(mmcif_dict, tokenizecifstructure(input))
    # Define ProteinStructure and add to it incrementally
    struc = ProteinStructure(structure_name)
    if haskey(mmcif_dict, "_atom_site.id")
        for i in 1:length(mmcif_dict["_atom_site.id"])
            if (read_std_atoms && mmcif_dict["_atom_site.group_PDB"][i] == "ATOM") ||
                    (read_het_atoms && mmcif_dict["_atom_site.group_PDB"][i] == "HETATM")
                model_n = parse(Int, mmcif_dict["_atom_site.pdbx_PDB_model_num"][i])
                # Create model if required
                if !haskey(models(struc), model_n)
                    struc[model_n] = Model(model_n, struc)
                end
                unsafe_addatomtomodel!(
                    struc[model_n],
                    AtomRecord(mmcif_dict, i),
                    remove_disorder=remove_disorder)
            end
        end
        # Generate lists for iteration
        fixlists!(struc)
    end
    return struc
end


# Constructor from mmCIF ATOM/HETATM line
AtomRecord(d::MMCIFDict, i::Integer) = AtomRecord(
    d["_atom_site.group_PDB"][i] == "HETATM",
    parse(Int, d["_atom_site.id"][i]),
    d["_atom_site.auth_atom_id"][i],
    d["_atom_site.label_alt_id"][i] in missingvals ? ' ' : d["_atom_site.label_alt_id"][i][1],
    d["_atom_site.auth_comp_id"][i],
    d["_atom_site.auth_asym_id"][i],
    parse(Int, d["_atom_site.auth_seq_id"][i]),
    d["_atom_site.pdbx_PDB_ins_code"][i] in missingvals ? ' ' : d["_atom_site.pdbx_PDB_ins_code"][i][1],
    [
        parse(Float64, d["_atom_site.Cartn_x"][i]),
        parse(Float64, d["_atom_site.Cartn_y"][i]),
        parse(Float64, d["_atom_site.Cartn_z"][i])
    ],
    d["_atom_site.occupancy"][i] in missingvals ? 1.0 : parse(Float64, d["_atom_site.occupancy"][i]),
    d["_atom_site.B_iso_or_equiv"][i] in missingvals ? 0.0 : parse(Float64, d["_atom_site.B_iso_or_equiv"][i]),
    d["_atom_site.type_symbol"][i] in missingvals ? "  " : d["_atom_site.type_symbol"][i],
    d["_atom_site.pdbx_formal_charge"][i] in missingvals ? "  " : d["_atom_site.pdbx_formal_charge"][i],
)


# Format a mmCIF data value by enclosing with quotes or semicolon lines where
#   appropriate. See
#   https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax for syntax.
function formatmmcifcol(val::AbstractString, col_width::Integer=length(val))
    # If there is a newline or quotes cannot be contained, use semicolon
    #   and newline construct
    if requiresnewline(val)
        return "\n;$val\n;\n"
    elseif requiresquote(val)
        # Choose quote character
        if occursin("' ", val)
            return rpad("\"$val\"", col_width)
        else
            return rpad("'$val'", col_width)
        end
    # Safe to not quote
    # Numbers must not be quoted
    else
        return rpad(val, col_width)
    end
end

function requiresnewline(val)
    return occursin("\n", val) || (occursin("' ", val) && occursin("\" ", val))
end

function requiresquote(val)
    return occursin(" ", val) || occursin("'", val) || occursin("\"", val) ||
        val[1] in specialchars || startswith(lowercase(val), "data_") ||
        startswith(lowercase(val), "save_") || val in specialwords
end


"""
    writemmcif(output, element, atom_selectors...)
    writemmcif(output, mmcif_dict)

Write a `StructuralElementOrList` or a `MMCIFDict` to a mmCIF format file
or output stream.

Atom selector functions can be given as additional arguments - only atoms
that return `true` from all the functions are retained.
"""
function writemmcif(filepath::AbstractString, mmcif_dict::MMCIFDict)
    open(filepath, "w") do output
        writemmcif(output, mmcif_dict)
    end
end

function writemmcif(output::IO, mmcif_dict::MMCIFDict)
    # Form dictionary where key is first part of mmCIF key and value is list
    #   of corresponding second parts
    key_lists = Dict{String, Vector{String}}()
    data_val = String[]
    for key in keys(mmcif_dict)
        if lowercase(key) == "data_"
            data_val = mmcif_dict[key]
        else
            s = split(key, ".")
            if length(s) == 2
                if s[1] in keys(key_lists)
                    push!(key_lists[s[1]], s[2])
                else
                    key_lists[s[1]] = [s[2]]
                end
            else
                throw(ArgumentError("Invalid key in mmCIF dictionary: $key"))
            end
        end
    end

    # Re-order lists if an order has been specified
    # Not all elements from the specified order are necessarily present
    for (key, key_list) in key_lists
        if key in keys(mmciforder)
            inds = Int[]
            for i in key_list
                f = something(findfirst(isequal(i), mmciforder[key]), 0)
                if f == 0
                    # Unrecognised key - add at end
                    f = length(mmciforder[key]) + 1
                end
                push!(inds, f)
            end
            key_lists[key] = [k for (_, k) in sort(collect(zip(inds, key_list)))]
        end
    end

    # Write out top data_ line
    if data_val != String[]
        println(output, "data_$(first(data_val))\n#")
    else
        println(output, "data_unknown\n#")
    end

    for (key, key_list) in key_lists
        # Pick a sample mmCIF value
        sample_val = mmcif_dict["$key.$(first(key_list))"]
        n_vals = length(sample_val)
        # Check the mmCIF dictionary has consistent list sizes
        for i in key_list
            if length(mmcif_dict["$key.$i"]) != n_vals
                throw(ArgumentError("Inconsistent list sizes in mmCIF dictionary: $key.$i"))
            end
        end
        # If the value has a single component, write as key-value pairs
        if length(sample_val) == 1
            # Find the maximum key length
            m = maximum(length.(key_list))
            for i in key_list
                println(output, "$(rpad("$key.$i", length(key)+m+4))$(formatmmcifcol(first(mmcif_dict["$key.$i"])))")
            end
        # If the value has more than one component, write as keys then a value table
        else
            println(output, "loop_")
            col_widths = Dict{String, Int}()
            # Write keys and find max widths for each set of values
            # Technically the max of the sum of the column widths is 2048
            for i in key_list
                println(output, "$key.$i")
                col_widths[i] = 0
                for val in mmcif_dict["$key.$i"]
                    len_val = length(val)
                    # If the value requires quoting it will add 2 characters
                    if requiresquote(val) && !requiresnewline(val)
                        len_val += 2
                    end
                    if len_val > col_widths[i]
                        col_widths[i] = len_val
                    end
                end
            end

            # Write the values as rows
            for i in 1:n_vals
                for col in key_list
                    print(output, formatmmcifcol(mmcif_dict["$key.$col"][i], col_widths[col]+1))
                end
                println(output)
            end
        end
        println(output, "#")
    end
end

function writemmcif(filepath::AbstractString,
                el::StructuralElementOrList,
                atom_selectors::Function...)
    open(filepath, "w") do output
        writemmcif(output, el, atom_selectors...)
    end
end

function writemmcif(output::IO,
                el::Union{ProteinStructure, Model, Chain, AbstractResidue,
                    Vector{Model}, Vector{Chain}, Vector{AbstractResidue},
                    Vector{Residue}, Vector{DisorderedResidue}},
                atom_selectors::Function...)
    # Create an empty dictionary and add atoms one at a time
    atom_dict = Dict{String, Vector{String}}(
            ["_atom_site.$i"=> String[] for i in mmciforder["_atom_site"]])
    atom_dict["data_"] = [structurename(first(el))]

    # Ensure multiple models get written out correctly
    loop_el = el
    if isa(el, ProteinStructure)
        loop_el = collectmodels(el)
    end

    # Collect residues then expand out disordered residues
    for res in collectresidues(loop_el)
        model_n = string(modelnumber(res))
        chain_id = strip(chainid(res)) == "" ? "." : chainid(res)
        res_n = string(resnumber(res))
        het = ishetero(res) ? "HETATM" : "ATOM"
        if isa(res, Residue)
            res_name = resname(res)
            for at in collectatoms(res, atom_selectors...)
                for atom_record in at
                    appendatom!(atom_dict, atom_record, model_n, chain_id,
                        res_n, res_name, het)
                end
            end
        else
            for res_name in resnames(res)
                for at in collectatoms(disorderedres(res, res_name), atom_selectors...)
                    for atom_record in at
                        appendatom!(atom_dict, atom_record, model_n, chain_id,
                            res_n, res_name, het)
                    end
                end
            end
        end
    end

    # Now the MMCIFDict has been generated, write it out to the file
    return writemmcif(output, MMCIFDict(atom_dict))
end

function writemmcif(output::IO, at::AbstractAtom, atom_selectors::Function...)
    atom_dict = Dict{String, Vector{String}}(
            ["_atom_site.$i"=> String[] for i in mmciforder["_atom_site"]])
    atom_dict["data_"] = [structurename(at)]
    for atom_record in at
        appendatom!(atom_dict, atom_record, string(modelnumber(at)),
            strip(chainid(at)) == "" ? "." : chainid(at),
            string(resnumber(at)), resname(at),
            ishetero(at) ? "HETATM" : "ATOM")
    end
    return writemmcif(output, MMCIFDict(atom_dict))
end

function writemmcif(output::IO,
                ats::Vector{<:AbstractAtom},
                atom_selectors::Function...)
    atom_dict = Dict{String, Vector{String}}(
            ["_atom_site.$i"=> String[] for i in mmciforder["_atom_site"]])
    atom_dict["data_"] = [structurename(ats[1])]
    for at in collectatoms(ats, atom_selectors...)
        for atom_record in at
            appendatom!(atom_dict, atom_record, string(modelnumber(at)),
                strip(chainid(at)) == "" ? "." : chainid(at),
                string(resnumber(at)), resname(at),
                ishetero(at) ? "HETATM" : "ATOM")
        end
    end
    return writemmcif(output, MMCIFDict(atom_dict))
end


# Add an atom record to a growing atom dictionary
function appendatom!(atom_dict, at, model_n, chain_id, res_n, res_name, het)
    push!(atom_dict["_atom_site.group_PDB"], het)
    push!(atom_dict["_atom_site.id"], string(serial(at)))
    push!(atom_dict["_atom_site.type_symbol"], length(element(at)) == 0 ? "?" : element(at))
    push!(atom_dict["_atom_site.label_atom_id"], atomname(at))
    push!(atom_dict["_atom_site.label_alt_id"], altlocid(at) == ' ' ? "." : string(altlocid(at)))
    push!(atom_dict["_atom_site.label_comp_id"], res_name)
    push!(atom_dict["_atom_site.label_asym_id"], "?")
    push!(atom_dict["_atom_site.label_entity_id"], "?")
    push!(atom_dict["_atom_site.label_seq_id"], "?")
    push!(atom_dict["_atom_site.pdbx_PDB_ins_code"], inscode(at) == ' ' ? "?" : string(inscode(at)))
    push!(atom_dict["_atom_site.Cartn_x"], pyfmt(coordspec, round(x(at), digits=3)))
    push!(atom_dict["_atom_site.Cartn_y"], pyfmt(coordspec, round(y(at), digits=3)))
    push!(atom_dict["_atom_site.Cartn_z"], pyfmt(coordspec, round(z(at), digits=3)))
    push!(atom_dict["_atom_site.occupancy"], pyfmt(floatspec, occupancy(at)))
    push!(atom_dict["_atom_site.B_iso_or_equiv"], pyfmt(floatspec, tempfactor(at)))
    push!(atom_dict["_atom_site.pdbx_formal_charge"], length(charge(at)) == 0 ? "?" : charge(at))
    push!(atom_dict["_atom_site.auth_seq_id"], res_n)
    push!(atom_dict["_atom_site.auth_comp_id"], res_name)
    push!(atom_dict["_atom_site.auth_asym_id"], chain_id)
    push!(atom_dict["_atom_site.auth_atom_id"], atomname(at))
    push!(atom_dict["_atom_site.pdbx_PDB_model_num"], model_n)
end
