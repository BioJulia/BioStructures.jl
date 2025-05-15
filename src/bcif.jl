export BCIFDict, BCIFFormat
using LinearAlgebra
import MsgPack


"""
    read(input::IO, ::Type{BCIFFormat}, structure_name::AbstractString="", remove_disorder::Bool=false, read_std_atoms::Bool=true, read_het_atoms::Bool=true, run_dssp::Bool=false, run_stride::Bool=false)

A function to read a binary CIF file from MolStar and extract the list of attributes and their compressed bytes.
"""
function Base.read(io::IO,
            ::Type{BCIFFormat};
            structure_name::AbstractString="",
            remove_disorder::Bool=false,
            read_std_atoms::Bool=true,
            read_het_atoms::Bool=true,
            run_dssp::Bool=false,
            run_stride::Bool=false)
    
    bcif_dict = MsgPack.unpack(io, BCIFDict)
    
    return MolecularStructure(
        bcif_dict; 
        structure_name=structure_name,
        remove_disorder=remove_disorder,
        read_std_atoms=read_std_atoms,
        read_het_atoms=read_het_atoms,
        run_dssp=run_dssp,
        run_stride=run_stride
    )
end
function Base.read(filename::AbstractString, ::Type{BCIFFormat}; kwargs...)
    open(filename) do io
        read(io, BCIFFormat; kwargs...)
    end
end

struct BCIFDict <: AbstractDict{String, Any}
    dict::Dict{String,Any}
    
    function BCIFDict(dict::Dict{String, Dict{String}})
        new(convert(Dict{String,Any}, dict))
    end
end

Base.keys(mmcif_dict::BCIFDict) = keys(mmcif_dict.dict)
Base.values(mmcif_dict::BCIFDict) = values(mmcif_dict.dict)
Base.haskey(mmcif_dict::BCIFDict, key) = haskey(mmcif_dict.dict, key)
Base.get(mmcif_dict::BCIFDict, key, default) = get(mmcif_dict.dict, key, default)
Base.length(mmcif_dict::BCIFDict) = length(mmcif_dict.dict)
Base.iterate(mmcif_dict::BCIFDict) = iterate(mmcif_dict.dict)
Base.iterate(mmcif_dict::BCIFDict, i) = iterate(mmcif_dict.dict, i)

function MolecularStructure(bcif_dict::BCIFDict; 
    structure_name::AbstractString="",
    remove_disorder::Bool=false,
    read_std_atoms::Bool=true,
    read_het_atoms::Bool=true,
    run_dssp::Bool=false,
    run_stride::Bool=false)

    struc = MolecularStructure(structure_name)

    for i in 1:length(bcif_dict["_atom_site"]["id"])
        model_n = bcif_dict["_atom_site"]["pdbx_PDB_model_num"][i]
        if !haskey(models(struc), model_n)
            struc[model_n] = Model(model_n, struc)
        end
        unsafe_addatomtomodel!(struc[model_n], AtomRecord(bcif_dict, i))
    end
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

function AtomRecord(d::BCIFDict, i::Int) 
    d = d["_atom_site"]
    return AtomRecord(
        d["group_PDB"][i] == "HETATM",
        d["id"][i],
        d["auth_atom_id"][i],
        d["label_alt_id"][i] == "" ? ' ' : d["label_alt_id"][i][1],
        d["auth_comp_id"][i],
        string(d["auth_asym_id"][i]),
        d["auth_seq_id"][i],
        d["pdbx_PDB_ins_code"][i] == "" ? ' ' : d["pdbx_PDB_ins_code"][i][1],
        [
            d["Cartn_x"][i],
            d["Cartn_y"][i],
            d["Cartn_z"][i]
        ],
        d["occupancy"][i],
        d["B_iso_or_equiv"][i],
        d["type_symbol"][i],
        string(d["pdbx_formal_charge"][i]),
    )
end

# Utility functions for encoding/decoding

function MsgPack.from_msgpack(::Type{BCIFDict}, data::Dict{String, Any})
    categories = data["dataBlocks"][1]["categories"]
    return BCIFDict(reduce(merge, [Dict(category["name"] => columns_to_dict(category["columns"])) for category in categories]))
end

function columns_to_dict(columns::Vector{Any})
    tasks = map(columns) do column
        Threads.@spawn Dict(column["name"] => decode_column(column))
    end
    return reduce(merge, fetch.(tasks))
end

function encode_stepwise(data, encodings)
    for encoding in encodings
        data = encode(encoding, data)
    end
    return data
end

function decode_stepwise(data, encodings)
    for encoding in reverse(encodings)
        data = decode(encoding, data)
    end
    return data
end

function deserialize_numeric_encoding(content::Any)
    if isa(content, Vector)
        return [deserialize_numeric_encoding(item) for item in content]
    end

    if isa(content, Encoding)
        return content
    end
    kind = content["kind"]

    # if byte convert to integer
    for (key, value) in content
        content[key] = value isa UInt8 ? Int32(value) : value
    end
    params = content

    if kind == "ByteArray"
        return ByteArrayEncoding(INT_TO_TYPE[params["type"]])
    elseif kind == "FixedPoint"
        return FixedPointEncoding(params["factor"]; srcType=INT_TO_TYPE[params["srcType"]])
    elseif kind == "IntervalQuantization"
        return IntervalQuantizationEncoding(params["min"], params["max"], params["numSteps"]; srcType=INT_TO_TYPE[params["srcType"]])
    elseif kind == "RunLength"
        return RunLengthEncoding(srcSize=params["srcSize"], srcType=INT_TO_TYPE[params["srcType"]])
    elseif kind == "Delta"
        # Pass the actual integer type code, not the Type object
        return DeltaEncoding(srcType=INT_TO_TYPE[params["srcType"]], origin=Int32(params["origin"]))
    elseif kind == "IntegerPacking"
        return IntegerPackingEncoding(params["byteCount"], srcSize=params["srcSize"], isUnsigned=params["isUnsigned"])
    end
end


function decode_column(column::Dict)
    column_data = column["data"]
    encodings = []

    # collect the encodings. If it's a string encoding then it should be a single encoding
    # that contains it's own dataEncoding and offsetEncoding which also need to be handled
    for enc in column_data["encoding"]
        if enc["kind"] == "StringArray"
            push!(encodings, StringArrayEncoding(
                stringData=enc["stringData"],
                dataEncoding=deserialize_numeric_encoding(enc["dataEncoding"]),
                offsetEncoding=deserialize_numeric_encoding(enc["offsetEncoding"]),
                offsets=enc["offsets"]
            ))
        else
            push!(encodings, deserialize_numeric_encoding(enc))
        end
    end

    return decode_stepwise(column_data["data"], encodings)
end


# Below are the encoding and decoding types for BCIF format

# Data types defined for the BCIF encoding by are indicated by integer values
# there are not well discussed in the official spec, had to ask about it explicitly
# https://github.com/molstar/BinaryCIF/issues/4
const INT_TO_TYPE = Dict(
    1 => Int8,
    2 => Int16,
    3 => Int32,
    4 => UInt8,
    5 => UInt16,
    6 => UInt32,
    32 => Float32,
    33 => Float64
)

const TYPES_TO_INT = IdDict{Any, Int}(t => i for (i, t) in INT_TO_TYPE)

const EncodingDataTypes = Union{values(INT_TO_TYPE)...}

# Mapping from Julia types to TypeCode

# Safe casting function
# function safe_cast(array, dtype)
#     if eltype(array) == dtype
#         return array
#     end

#     if dtype <: Integer && !(eltype(array) <: Integer)
#         throw(ArgumentError("Cannot cast floating point to integer"))
#     end

#     if dtype <: Integer
#         type_min, type_max = typemin(dtype), typemax(dtype)
#         if any(x -> x < type_min || x > type_max, array)
#             throw(ArgumentError("Integer values do not fit into the given dtype"))
#         end
#     end

#     return convert(Array{dtype}, array)
# end

# Abstract encoding type
abstract type Encoding end

# ByteArrayEncoding
mutable struct ByteArrayEncoding <: Encoding
    type::DataType
end

# Add a constructor that handles the type lookup from INT_TO_TYPE
function ByteArrayEncoding(type_code::Integer)
    if haskey(INT_TO_TYPE, type_code)
        return ByteArrayEncoding(INT_TO_TYPE[type_code])
    else
        throw(ArgumentError("Invalid type code: $type_code"))
    end
end

# function encode(enc::ByteArrayEncoding, data)
#     if enc.type === nothing
#         enc.type = TYPE_TO_TYPE_CODE[eltype(data)]
#     end
#     return reinterpret(UInt8, safe_cast(data, TYPE_CODE_TO_TYPE[enc.type]))
# end

function decode(enc::ByteArrayEncoding, data)
    return reinterpret(enc.type, data)
end

# FixedPointEncoding
mutable struct FixedPointEncoding <: Encoding
    factor::Float64
    srcType::DataType

    function FixedPointEncoding(factor; srcType::DataType)
        new(factor, srcType)
    end
end

function encode(enc::FixedPointEncoding, data)
    return round.(Int32, data .* enc.factor)
end

function decode(enc::FixedPointEncoding, data)
    return convert(Array{enc.srcType}, data ./ enc.factor)
end

# IntervalQuantizationEncoding
mutable struct IntervalQuantizationEncoding <: Encoding
    min::Float64
    max::Float64
    numSteps::Int
    srcType::DataType
end

# function encode(enc::IntervalQuantizationEncoding, data)
#     # Convert to normalized values between 0 and numSteps-1
#     normalized = (data .- enc.min) ./ (enc.max - enc.min) .* (enc.numSteps - 1)
#     # Clamp to valid range and convert to integers
#     indices = clamp.(round.(Int32, normalized), 0, enc.numSteps - 1)
#     return indices
# end

function decode(enc::IntervalQuantizationEncoding, data)
    # Convert indices back to values in the original range
    normalized = data ./ Float64(enc.numSteps - 1)
    output = normalized .* (enc.max - enc.min) .+ enc.min
    return convert(Array{enc.srcType}, output)
end

# RunLengthEncoding
mutable struct RunLengthEncoding <: Encoding
    srcSize::Int
    srcType::DataType

    function RunLengthEncoding(; srcSize, srcType)
        new(srcSize, srcType)
    end
end

# function encode(enc::RunLengthEncoding, data)
#     if enc.srcType === nothing
#         enc.srcType = TYPE_TO_TYPE_CODE[eltype(data)]
#     end
#     if enc.srcSize === nothing
#         enc.srcSize = length(data)
#     elseif enc.srcSize != length(data)
#         throw(ArgumentError("Given source size does not match actual data size"))
#     end

#     # Pessimistic allocation - worst case is run length of 1 for every element
#     output = zeros(Int32, length(data) * 2)
#     j = 1
#     val = data[1]
#     run_length = 0

#     for i in 1:length(data)
#         curr_val = data[i]
#         if curr_val == val
#             run_length += 1
#         else
#             # New element -> Write element with run-length
#             output[j] = val
#             output[j+1] = run_length
#             j += 2
#             val = curr_val
#             run_length = 1
#         end
#     end

#     # Write last element
#     output[j] = val
#     output[j+1] = run_length
#     j += 2

#     # Trim to correct size
#     return output[1:j-1]
# end

function decode(enc::RunLengthEncoding, data)
    if length(data) % 2 != 0
        throw(ArgumentError("Invalid run-length encoded data"))
    end

    length_output = 0
    if enc.srcSize === nothing
        # Determine length of output array by summing run lengths
        for i in 2:2:length(data)
            length_output += data[i]
        end
    else
        length_output = enc.srcSize
    end

    output = zeros(enc.srcType, length_output)
    j = 1

    for i in 1:2:length(data)
        value = data[i]
        repeat_count = data[i+1]
        output[j:j+repeat_count-1] .= value
        j += repeat_count
    end

    return output
end

# DeltaEncoding
mutable struct DeltaEncoding <: Encoding
    srcType::DataType
    origin::Int32

    # Constructor for Type parameter
    function DeltaEncoding(; srcType::Type, origin::Int32=0)
        new(srcType, origin)
    end
end

function decode(enc::DeltaEncoding, data)
    output = cumsum(data)
    output = convert(Array{enc.srcType}, output)
    output .+= enc.origin
    return output
end

# IntegerPackingEncoding
mutable struct IntegerPackingEncoding <: Encoding
    byteCount::Int
    srcSize::Int
    isUnsigned::Bool

    function IntegerPackingEncoding(byteCount; srcSize, isUnsigned=false)
        new(byteCount, srcSize, isUnsigned)
    end
end

function determine_packed_dtype(enc::IntegerPackingEncoding)
    if enc.byteCount == 1
        return enc.isUnsigned ? UInt8 : Int8
    elseif enc.byteCount == 2
        return enc.isUnsigned ? UInt16 : Int16
    else
        throw(ArgumentError("Unsupported byte count"))
    end
end

# function encode(enc::IntegerPackingEncoding, data)
#     if enc.srcSize === nothing
#         enc.srcSize = length(data)
#     elseif enc.srcSize != length(data)
#         throw(ArgumentError("Given source size does not match actual data size"))
#     end

#     data = convert(Array{Int32}, data)
#     packed_type = determine_packed_dtype(enc)
#     min_val = typemin(packed_type)
#     max_val = typemax(packed_type)

#     # Get length of output array by summing up required length of each element
#     length_output = 0
#     for num in data
#         if num < 0
#             if min_val == 0
#                 throw(ArgumentError("Cannot pack negative numbers into unsigned type"))
#             end
#             # Required packed length is number of times min_val needs to be repeated + 1
#             length_output += div(num, min_val) + 1
#         elseif num > 0
#             length_output += div(num, max_val) + 1
#         else
#             # num = 0
#             length_output += 1
#         end
#     end

#     # Fill output
#     output = zeros(packed_type, length_output)
#     j = 1

#     for i in 1:length(data)
#         remainder = data[i]
#         if remainder < 0
#             if min_val == 0
#                 throw(ArgumentError("Cannot pack negative numbers into unsigned type"))
#             end
#             while remainder <= min_val
#                 remainder -= min_val
#                 output[j] = min_val
#                 j += 1
#             end
#         elseif remainder > 0
#             while remainder >= max_val
#                 remainder -= max_val
#                 output[j] = max_val
#                 j += 1
#             end
#         end
#         output[j] = remainder
#         j += 1
#     end

#     return output
# end

function decode(enc::IntegerPackingEncoding, data)
    packed_type = determine_packed_dtype(enc)
    min_val = typemin(packed_type)
    max_val = typemax(packed_type)

    # For unsigned integers, do not check lower bound (is always 0)
    # -> Set lower bound to value that is never reached
    if min_val == 0
        min_val = -1
    end

    output = zeros(Int32, enc.srcSize)
    j = 1
    unpacked_val = 0

    for i in 1:length(data)
        packed_val = data[i]
        if packed_val == max_val || packed_val == min_val
            unpacked_val += packed_val
        else
            unpacked_val += packed_val
            output[j] = unpacked_val
            unpacked_val = 0
            j += 1
        end
    end

    return output
end

# StringArrayEncoding
mutable struct StringArrayEncoding <: Encoding
    stringData::String
    dataEncoding::Vector{Encoding}
    offsetEncoding::Vector{Encoding}
    offsets::Vector{UInt8}

    function StringArrayEncoding(; stringData=nothing, dataEncoding=nothing, offsetEncoding=nothing, offsets=nothing)
        if dataEncoding === nothing
            dataEncoding = [ByteArrayEncoding(Int32)]
        end
        if offsetEncoding === nothing
            offsetEncoding = [ByteArrayEncoding(Int32)]
        end
        new(stringData, dataEncoding, offsetEncoding, offsets)
    end
end

# function encode(enc::StringArrayEncoding, data)
#     if !(eltype(data) <: AbstractString)
#         throw(ArgumentError("Data must be of string type"))
#     end

#     if enc.stringData === nothing
#         # Get unique stringData
#         enc.stringData = unique(data)
#         check_present = false
#     else
#         check_present = true
#     end

#     # Sort stringData for binary search
#     sorted_indices = sortperm(enc.stringData)
#     sorted_strings = enc.stringData[sorted_indices]

#     # Find indices of each string in data
#     indices = zeros(Int32, length(data))
#     for i in 1:length(data)
#         idx = searchsortedfirst(sorted_strings, data[i])
#         if idx <= length(sorted_strings) && sorted_strings[idx] == data[i]
#             indices[i] = sorted_indices[idx]
#         else
#             if check_present
#                 throw(ArgumentError("Data contains stringData not present in 'stringData'"))
#             end
#         end
#     end

#     # Apply encodings
#     encoded_data = indices
#     for encoding in enc.dataEncoding
#         encoded_data = encode(encoding, encoded_data)
#     end

#     return encoded_data
# end

function decode(enc::StringArrayEncoding, data)
    # Apply decodings in reverse order
    indices = decode_stepwise(data, enc.dataEncoding) .+ 1
    offsets = decode_stepwise(enc.offsets, enc.offsetEncoding)

    substrings = Vector{String}()

    # break up the string into the substrings that are individual occurrences
    for (i, offset) in enumerate(offsets[1:end-1])
        start_i = offsets[i] + 1
        end_i = offsets[i+1]
        push!(substrings, String(enc.stringData[start_i:end_i]))
    end

    return substrings[indices]
end
