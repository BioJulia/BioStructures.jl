export decode_column
using LinearAlgebra
import MsgPack


"""
    read(input::IO, ::Type{BCIFFormat}, structure_name::AbstractString="", remove_disorder::Bool=false, read_std_atoms::Bool=true, read_het_atoms::Bool=true, run_dssp::Bool=false, run_stride::Bool=false)

A function to read a binary CIF file from MolStar and extract the list of attributes and their compressed bytes.
"""
function Base.read(input::IO,
    ::Type{BCIFFormat},
    structure_name::AbstractString="",
    remove_disorder::Bool=false,
    read_std_atoms::Bool=true,
    read_het_atoms::Bool=true,
    run_dssp::Bool=false,
    run_stride::Bool=false)

    file = MsgPack.unpack(read(input))

    # currently just looking for the first data block
    categories = file["dataBlocks"][1]["categories"]
    atom_site = get_category(categories, "_atom_site")
    columns = atom_site[1]["columns"]
    tasks = map(columns) do column
        Threads.@spawn decode_column(column)
    end

    bcif_dict = BCIFDict(Dict(columns[i]["name"] => result for (i, result) in enumerate(fetch.(tasks))))
    # return bcif_dict
    struc = MolecularStructure(structure_name)
    struc[1] = Model(1, struc)
    for i in 1:length(bcif_dict["id"])
        unsafe_addatomtomodel!(struc[1], AtomRecord(bcif_dict, i))
    end
    fixlists!(struc)
    return struc
end

BCIFArrayTypes = Union{Vector{String},Vector{Int32},Vector{Float64}}

struct BCIFDict <: AbstractDict{String,BCIFArrayTypes}
    dict::Dict{String,BCIFArrayTypes}
end

function BCIFDict(dict::Dict{String,BCIFArrayTypes})
    new(dict)
end

Base.keys(mmcif_dict::BCIFDict) = keys(mmcif_dict.dict)
Base.values(mmcif_dict::BCIFDict) = values(mmcif_dict.dict)
Base.haskey(mmcif_dict::BCIFDict, key) = haskey(mmcif_dict.dict, key)
Base.get(mmcif_dict::BCIFDict, key, default) = get(mmcif_dict.dict, key, default)
Base.length(mmcif_dict::BCIFDict) = length(mmcif_dict.dict)
Base.iterate(mmcif_dict::BCIFDict) = iterate(mmcif_dict.dict)
Base.iterate(mmcif_dict::BCIFDict, i) = iterate(mmcif_dict.dict, i)

AtomRecord = AtomRecord(d::BCIFDict, i::Int) = AtomRecord(
    d["group_PDB"][i] == "HETATM",
    d["id"][i],
    d["auth_atom_id"][i],
    'A',# d["label_alt_id"][i],
    d["auth_comp_id"][i],
    d["auth_asym_id"][i],
    d["auth_seq_id"][i],
    'A', # d["pdbx_PDB_ins_code"][i],
    [
        d["Cartn_x"][i],
        d["Cartn_y"][i],
        d["Cartn_z"][i]
    ],
    d["occupancy"][i],
    d["B_iso_or_equiv"][i],
    d["type_symbol"][i],
    d["pdbx_formal_charge"][i],
)



function get_category(cats::Vector{Any}, name::String)
    idx = findall(getindex.(cats, "name") .== name)

    if isnothing(idx)
        throw(ArgumentError("Category $name not found"))
    end
    if length(idx) > 1
        throw(ArgumentError("Multiple categories with name $name found"))
    end

    return cats[idx]
end

# Enum for type codes
@enum TypeCode begin
    INT8 = 1
    INT16 = 2
    INT32 = 3
    UINT8 = 4
    UINT16 = 5
    UINT32 = 6
    FLOAT32 = 32
    FLOAT64 = 33
end

# Mapping from TypeCode to Julia types
const TYPE_CODE_TO_TYPE = Dict(
    INT8 => Int8,
    INT16 => Int16,
    INT32 => Int32,
    UINT8 => UInt8,
    UINT16 => UInt16,
    UINT32 => UInt32,
    FLOAT32 => Float32,
    FLOAT64 => Float64
)

# Mapping from Julia types to TypeCode
const TYPE_TO_TYPE_CODE = Dict(value => key for (key, value) in TYPE_CODE_TO_TYPE)

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


# Safe casting function
function safe_cast(array, dtype)
    if eltype(array) == dtype
        return array
    end

    if dtype <: Integer && !(eltype(array) <: Integer)
        throw(ArgumentError("Cannot cast floating point to integer"))
    end

    if dtype <: Integer
        type_min, type_max = typemin(dtype), typemax(dtype)
        if any(x -> x < type_min || x > type_max, array)
            throw(ArgumentError("Integer values do not fit into the given dtype"))
        end
    end

    return convert(Array{dtype}, array)
end

# Abstract encoding type
abstract type Encoding end

# ByteArrayEncoding
mutable struct ByteArrayEncoding <: Encoding
    type::Union{TypeCode,Nothing}

    function ByteArrayEncoding(type=nothing)
        if type !== nothing
            type = type isa TypeCode ? type : TYPE_TO_TYPE_CODE[type]
        end
        new(type)
    end
end

function encode(enc::ByteArrayEncoding, data)
    if enc.type === nothing
        enc.type = TYPE_TO_TYPE_CODE[eltype(data)]
    end
    return reinterpret(UInt8, safe_cast(data, TYPE_CODE_TO_TYPE[enc.type]))
end

function decode(enc::ByteArrayEncoding, data)
    return reinterpret(TYPE_CODE_TO_TYPE[enc.type], data)
end

# FixedPointEncoding
mutable struct FixedPointEncoding <: Encoding
    factor::Float64
    srcType::TypeCode

    function FixedPointEncoding(factor; srcType=FLOAT32)
        srcType = srcType isa TypeCode ? srcType : TYPE_TO_TYPE_CODE[srcType]
        if !(srcType in (FLOAT32, FLOAT64))
            throw(ArgumentError("Only floating point types are supported"))
        end
        new(factor, srcType)
    end
end

function encode(enc::FixedPointEncoding, data)
    return round.(Int32, data .* enc.factor)
end

function decode(enc::FixedPointEncoding, data)
    return convert(Array{TYPE_CODE_TO_TYPE[enc.srcType]}, data ./ enc.factor)
end

# IntervalQuantizationEncoding
mutable struct IntervalQuantizationEncoding <: Encoding
    min::Float64
    max::Float64
    numSteps::Int
    srcType::TypeCode

    function IntervalQuantizationEncoding(min, max, numSteps; srcType=FLOAT32)
        srcType = srcType isa TypeCode ? srcType : TYPE_TO_TYPE_CODE[srcType]
        new(min, max, numSteps, srcType)
    end
end

function encode(enc::IntervalQuantizationEncoding, data)
    # Convert to normalized values between 0 and numSteps-1
    normalized = (data .- enc.min) ./ (enc.max - enc.min) .* (enc.numSteps - 1)
    # Clamp to valid range and convert to integers
    indices = clamp.(round.(Int32, normalized), 0, enc.numSteps - 1)
    return indices
end

function decode(enc::IntervalQuantizationEncoding, data)
    # Convert indices back to values in the original range
    normalized = data ./ Float64(enc.numSteps - 1)
    output = normalized .* (enc.max - enc.min) .+ enc.min
    return convert(Array{TYPE_CODE_TO_TYPE[enc.srcType]}, output)
end

# RunLengthEncoding
mutable struct RunLengthEncoding <: Encoding
    srcSize::Union{Int,Nothing}
    srcType::Union{TypeCode,Nothing}

    function RunLengthEncoding(; srcSize=nothing, srcType=nothing)
        if srcType !== nothing
            srcType = srcType isa TypeCode ? srcType : TYPE_TO_TYPE_CODE[srcType]
        end
        new(srcSize, srcType)
    end
end

function encode(enc::RunLengthEncoding, data)
    if enc.srcType === nothing
        enc.srcType = TYPE_TO_TYPE_CODE[eltype(data)]
    end
    if enc.srcSize === nothing
        enc.srcSize = length(data)
    elseif enc.srcSize != length(data)
        throw(ArgumentError("Given source size does not match actual data size"))
    end

    # Pessimistic allocation - worst case is run length of 1 for every element
    output = zeros(Int32, length(data) * 2)
    j = 1
    val = data[1]
    run_length = 0

    for i in 1:length(data)
        curr_val = data[i]
        if curr_val == val
            run_length += 1
        else
            # New element -> Write element with run-length
            output[j] = val
            output[j+1] = run_length
            j += 2
            val = curr_val
            run_length = 1
        end
    end

    # Write last element
    output[j] = val
    output[j+1] = run_length
    j += 2

    # Trim to correct size
    return output[1:j-1]
end

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

    output = zeros(TYPE_CODE_TO_TYPE[enc.srcType], length_output)
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
    srcType::Union{TypeCode,Nothing}
    origin::Int

    function DeltaEncoding(; srcType=nothing, origin=0)
        if srcType !== nothing
            srcType = srcType isa TypeCode ? srcType : TYPE_TO_TYPE_CODE[srcType]
        end
        new(srcType, origin)
    end
end

function encode(enc::DeltaEncoding, data)
    if enc.srcType === nothing
        enc.srcType = TYPE_TO_TYPE_CODE[eltype(data)]
    end

    data = data .- enc.origin
    diffs = vcat([0], diff(data))
    return convert(Array{Int32}, diffs)
end

function decode(enc::DeltaEncoding, data)
    output = cumsum(data)
    output = convert(Array{TYPE_CODE_TO_TYPE[enc.srcType]}, output)
    output .+= enc.origin
    return output
end

# IntegerPackingEncoding
mutable struct IntegerPackingEncoding <: Encoding
    byteCount::Int
    srcSize::Union{Int,Nothing}
    isUnsigned::Bool

    function IntegerPackingEncoding(byteCount; srcSize=nothing, isUnsigned=false)
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

function encode(enc::IntegerPackingEncoding, data)
    if enc.srcSize === nothing
        enc.srcSize = length(data)
    elseif enc.srcSize != length(data)
        throw(ArgumentError("Given source size does not match actual data size"))
    end

    data = convert(Array{Int32}, data)
    packed_type = determine_packed_dtype(enc)
    min_val = typemin(packed_type)
    max_val = typemax(packed_type)

    # Get length of output array by summing up required length of each element
    length_output = 0
    for num in data
        if num < 0
            if min_val == 0
                throw(ArgumentError("Cannot pack negative numbers into unsigned type"))
            end
            # Required packed length is number of times min_val needs to be repeated + 1
            length_output += div(num, min_val) + 1
        elseif num > 0
            length_output += div(num, max_val) + 1
        else
            # num = 0
            length_output += 1
        end
    end

    # Fill output
    output = zeros(packed_type, length_output)
    j = 1

    for i in 1:length(data)
        remainder = data[i]
        if remainder < 0
            if min_val == 0
                throw(ArgumentError("Cannot pack negative numbers into unsigned type"))
            end
            while remainder <= min_val
                remainder -= min_val
                output[j] = min_val
                j += 1
            end
        elseif remainder > 0
            while remainder >= max_val
                remainder -= max_val
                output[j] = max_val
                j += 1
            end
        end
        output[j] = remainder
        j += 1
    end

    return output
end

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
            dataEncoding = [ByteArrayEncoding(INT32)]
        end
        if offsetEncoding === nothing
            offsetEncoding = [ByteArrayEncoding(INT32)]
        end
        new(stringData, dataEncoding, offsetEncoding, offsets)
    end
end

function encode(enc::StringArrayEncoding, data)
    if !(eltype(data) <: AbstractString)
        throw(ArgumentError("Data must be of string type"))
    end

    if enc.stringData === nothing
        # Get unique stringData
        enc.stringData = unique(data)
        check_present = false
    else
        check_present = true
    end

    # Sort stringData for binary search
    sorted_indices = sortperm(enc.stringData)
    sorted_strings = enc.stringData[sorted_indices]

    # Find indices of each string in data
    indices = zeros(Int32, length(data))
    for i in 1:length(data)
        idx = searchsortedfirst(sorted_strings, data[i])
        if idx <= length(sorted_strings) && sorted_strings[idx] == data[i]
            indices[i] = sorted_indices[idx]
        else
            if check_present
                throw(ArgumentError("Data contains stringData not present in 'stringData'"))
            end
        end
    end

    # Apply encodings
    encoded_data = indices
    for encoding in enc.dataEncoding
        encoded_data = encode(encoding, encoded_data)
    end

    return encoded_data
end

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

# Utility functions for encoding/decoding
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

function deserialize_encoding(content::Any)
    if isa(content, Vector)
        return [deserialize_encoding(item) for item in content]
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

    # Handle nested encodings
    if haskey(params, "data_encoding")
        params["data_encoding"] = deserialize_encoding(params["data_encoding"])
    end

    if haskey(params, "offsetEncoding")
        params["offsetEncoding"] = deserialize_encoding(params["offsetEncoding"])
    end

    encoding_constructors = Dict(
        "ByteArray" => () -> ByteArrayEncoding(INT_TO_TYPE[get(params, "type", nothing)]),
        "FixedPoint" => () -> FixedPointEncoding(params["factor"]; srcType=INT_TO_TYPE[get(params, "srcType", FLOAT32)]),
        "StringArray" => () -> StringArrayEncoding(
            stringData=get(params, "stringData", nothing),
            dataEncoding=get(params, "dataEncoding", nothing),
            offsetEncoding=get(params, "offsetEncoding", nothing),
            offsets=get(params, "offsets", nothing)
        ),
        "IntervalQuantization" => () -> IntervalQuantizationEncoding(params["min"], params["max"], params["numSteps"];
            srcType=INT_TO_TYPE[get(params, "srcType", 32)]),
        "RunLength" => () -> RunLengthEncoding(srcSize=get(params, "srcSize", nothing),
            srcType=INT_TO_TYPE[get(params, "srcType", nothing)]),
        "Delta" => () -> DeltaEncoding(srcType=INT_TO_TYPE[get(params, "srcType", nothing)],
            origin=get(params, "origin", 0)),
        "IntegerPacking" => () -> IntegerPackingEncoding(params["byteCount"],
            srcSize=get(params, "srcSize", nothing),
            isUnsigned=get(params, "isUnsigned", false))
    )

    if haskey(encoding_constructors, kind)
        return encoding_constructors[kind]()
    else
        error("Unknown encoding kind: $kind")
    end
end


function decode_column(column::Dict)
    data = column["data"]
    encodings = []

    # Handle the encoding array properly
    for enc in data["encoding"]
        if haskey(enc, "dataEncoding")
            if haskey(enc, "offsetEncoding")
                push!(encodings, StringArrayEncoding(
                    stringData=enc["stringData"],
                    dataEncoding=deserialize_encoding(enc["dataEncoding"]),
                    offsetEncoding=deserialize_encoding(enc["offsetEncoding"]),
                    offsets=enc["offsets"]
                ))
            else
                push!(encodings, deserialize_encoding(enc["dataEncoding"]))
            end
        else
            push!(encodings, deserialize_encoding(enc))
        end
    end

    # Flatten the encodings if needed
    flat_encodings = []
    for enc in encodings
        if enc isa Vector
            append!(flat_encodings, enc)
        else
            push!(flat_encodings, enc)
        end
    end


    # return flat_encodings

    decoded = decode_stepwise(data["data"], flat_encodings)


end
