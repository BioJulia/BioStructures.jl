export
    MMTFDict,
    writemmtf

"""
A MMTF dictionary.
Keys are field names as a `String` and values are various types.
"""
struct MMTFDict
    dict::Dict{String, Any}
end

# Create an empty MMTF dictionary
# Matches the decoded form of a MMTF file using MMTF.jl
# Encoding and decoding this Dict gives an identical Dict
MMTFDict() = MMTFDict(Dict{String, Any}(
        "altLocList"         => Char[],
        "atomIdList"         => Int32[],
        "bFactorList"        => Float32[],
        "bioAssemblyList"    => Any[],
        "bondAtomList"       => Int32[],
        "bondOrderList"      => Int8[],
        "chainIdList"        => String[],
        "chainNameList"      => String[],
        "chainsPerModel"     => Any[],
        "depositionDate"     => "",
        "entityList"         => Any[],
        "experimentalMethods"=> Any[],
        "groupIdList"        => Int32[],
        "groupList"          => Any[],
        "groupsPerChain"     => Any[],
        "groupTypeList"      => Int32[],
        "insCodeList"        => Char[],
        "mmtfProducer"       => "",
        "mmtfVersion"        => "",
        "ncsOperatorList"    => Any[],
        "numAtoms"           => 0,
        "numBonds"           => 0,
        "numChains"          => 0,
        "numGroups"          => 0,
        "numModels"          => 0,
        "occupancyList"      => Float32[],
        "releaseDate"        => "",
        "resolution"         => 0.0,
        "rFree"              => "",
        "rWork"              => "",
        "secStructList"      => Int8[],
        "sequenceIndexList"  => Int32[],
        "spaceGroup"         => "",
        "structureId"        => "",
        "title"              => "",
        "unitCell"           => Any[],
        "xCoordList"         => Float32[],
        "yCoordList"         => Float32[],
        "zCoordList"         => Float32[]
    ))

MMTFDict(filepath::AbstractString) = MMTFDict(parsemmtf(filepath))

Base.getindex(mmtf_dict::MMTFDict, field::AbstractString) = mmtf_dict.dict[field]

function Base.setindex!(mmtf_dict::MMTFDict,
                    val,
                    field::AbstractString)
    mmtf_dict.dict[field] = val
    return mmtf_dict
end

Base.keys(mmtf_dict::MMTFDict) = keys(mmtf_dict.dict)
Base.values(mmtf_dict::MMTFDict) = values(mmtf_dict.dict)
Base.haskey(mmtf_dict::MMTFDict, key) = haskey(mmtf_dict.dict, key)

function Base.show(io::IO, mmtf_dict::MMTFDict)
    print(io, "MMTF dictionary with $(length(keys(mmtf_dict))) fields")
end


function Base.read(input::IO,
            ::Type{MMTF};
            structure_name::AbstractString="",
            remove_disorder::Bool=false,
            read_std_atoms::Bool=true,
            read_het_atoms::Bool=true)
    d = parsemmtf(input)
    struc = ProteinStructure(structure_name)
    # Extract hetero atom information from entity list
    hets = trues(length(d["chainIdList"]))
    for e in d["entityList"]
        if e["type"] == "polymer"
            for i in e["chainIndexList"]
                # 0-based indexing in MMTF
                hets[i + 1] = false
            end
        end
    end
    model_i = 0
    chain_i = 0
    group_i = 0
    atom_i = 0
    for modelchaincount in d["chainsPerModel"]
        model_i += 1
        struc[model_i] = Model(model_i, struc)
        for ci in 1:modelchaincount
            chain_i += 1
            if (!read_std_atoms && !hets[chain_i]) || (!read_het_atoms && hets[chain_i])
                continue
            end
            for gi in 1:d["groupsPerChain"][chain_i]
                group_i += 1
                # 0-based indexing in MMTF
                group = d["groupList"][d["groupTypeList"][group_i] + 1]
                for ai in 1:length(group["atomNameList"])
                    atom_i += 1
                    unsafe_addatomtomodel!(
                        struc[model_i],
                        AtomRecord(
                            hets[chain_i],
                            d["atomIdList"][atom_i],
                            group["atomNameList"][ai],
                            d["altLocList"][atom_i] == '\0' ? ' ' : d["altLocList"][atom_i],
                            group["groupName"],
                            d["chainNameList"][chain_i],
                            d["groupIdList"][group_i],
                            d["insCodeList"][group_i] == '\0' ? ' ' : d["insCodeList"][group_i],
                            [
                                d["xCoordList"][atom_i],
                                d["yCoordList"][atom_i],
                                d["zCoordList"][atom_i],
                            ],
                            d["occupancyList"][atom_i],
                            d["bFactorList"][atom_i],
                            group["elementList"][ai],
                            # Add + to positive charges to match PDB convention
                            group["formalChargeList"][ai] > 0 ? "+$(group["formalChargeList"][ai])" :
                                        string(group["formalChargeList"][ai])
                        );
                        remove_disorder=remove_disorder)
                end
            end
        end
    end
    fixlists!(struc)
    return struc
end


function writemmtf(output::Union{AbstractString, IO},
                d::MMTFDict;
                gzip::Bool=false)
    writemmtf(d.dict, output; gzip=gzip)
    return
end

function writemmtf(filepath::AbstractString,
                el::StructuralElementOrList,
                atom_selectors::Function...;
                gzip::Bool=false)
    open(filepath, "w") do output
        writemmtf(output, el, atom_selectors...; gzip=gzip)
    end
end

generatechainid(i::Integer) = string(Char(64 + i))

function writemmtf(output::IO,
                el::StructuralElementOrList,
                atom_selectors::Function...;
                gzip::Bool=false)
    d = MMTFDict()
    # Index of standard and hetero atom entities; 0 means not yet added
    at_entity_i, het_entity_i = 0, 0
    group_names = String[]
    for mod in collectmodels(el)
        push!(d["chainsPerModel"], countchains(mod))
        chain_i = 0
        for ch in mod
            # MMTF splits chains up by molecules so we determine chain splits
            #   at the residue level
            prev_resname = ""
            prev_het = true
            group_count = 0
            for res in ch
                # Determine whether we have changed entity
                # ATOM blocks, and hetero molecules with the same name, are
                #   treated as the same entity
                if ishetero(res) != prev_het || (prev_het && resname(res) != prev_resname)
                    println(res)
                    chain_i += 1
                    push!(d["chainIdList"], generatechainid(chain_i))
                    push!(d["chainNameList"], chainid(ch))
                    # Add the groupsPerChain for the previous chain
                    if group_count > 0
                        push!(d["groupsPerChain"], group_count)
                        group_count = 0
                    end
                    # Check whether we need to add the entity for ATOM/HETATM
                    if (ishetero(res) && het_entity_i == 0) || (!ishetero(res) && at_entity_i == 0)
                        push!(d["entityList"], Dict{Any, Any}(
                            "chainIndexList"=> Any[],
                            "description"   => "",
                            "sequence"      => "",
                            "type"          => ishetero(res) ? "non-polymer" : "polymer"
                        ))
                        if ishetero(res)
                            het_entity_i = length(d["entityList"])
                        else
                            at_entity_i = length(d["entityList"])
                        end
                    entity_i = ishetero(res) ? het_entity_i : at_entity_i
                    push!(d["entityList"][entity_i]["chainIndexList"],
                        length(d["chainIdList"]) - 1)
                    end
                end
                group_count += 1

                if !(resname(res) in group_names)
                    push!(group_names, resname(res))
                    push!(d["groupList"], Dict{Any, Any}(
                        "groupName"       => resname(res),
                        "bondAtomList"    => Any[],
                        "elementList"     => Any[element(at) for at in res],
                        "formalChargeList"=> Any[parse(Int64, charge(at)) for at in res],
                        "singleLetterCode"=> "",
                        "chemCompType"    => "",
                        "atomNameList"    => Any[atomnames(res)...],
                        "bondOrderList"   => Any[]
                    ))
                end
                push!(d["groupIdList"], resnumber(res))
                push!(d["groupTypeList"], findfirst(isequal(resname(res)), group_names) - 1)
                push!(d["insCodeList"], inscode(res) == ' ' ? '\0' : inscode(res))
                for outer_at in collectatoms(res, atom_selectors...)
                    for at in outer_at
                        push!(d["altLocList"], altlocid(at) == ' ' ? '\0' : altlocid(at))
                        push!(d["atomIdList"], serial(at))
                        push!(d["bFactorList"], tempfactor(at))
                        push!(d["occupancyList"], occupancy(at))
                        push!(d["xCoordList"], x(at))
                        push!(d["yCoordList"], y(at))
                        push!(d["zCoordList"], z(at))
                    end
                end
                prev_resname = resname(res)
                prev_het = ishetero(res)
            end
            if group_count > 0
                push!(d["groupsPerChain"], group_count)
            end
        end
    end
    writemmtf(output, d; gzip=gzip)
end
