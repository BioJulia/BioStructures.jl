export
    MMTFDict,
    generatechainid,
    writemmtf

"""
    MMTFDict(filepath; gzip=false)
    MMTFDict(io; gzip=false)
    MMTFDict()

A Macromolecular Transmission Format (MMTF) dictionary.

Can be accessed using similar functions to a standard `Dict`.
Keys are field names as a `String` and values are various types.
To directly access the underlying dictionary of `MMTFDict` `d`, use
`d.dict`.
Call `MMTFDict` with a filepath or stream to read the dictionary from that
source.
The keyword argument `gzip` (default `false`) determines if the file is gzipped.
"""
struct MMTFDict <: AbstractDict{String, Any}
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
        "zCoordList"         => Float32[]))

MMTFDict(filepath::AbstractString; gzip::Bool=false) = MMTFDict(parsemmtf(filepath; gzip=gzip))
MMTFDict(io::IO; gzip::Bool=false) = MMTFDict(parsemmtf(io; gzip=gzip))

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
Base.get(mmtf_dict::MMTFDict, key, default) = get(mmtf_dict.dict, key, default)
Base.length(mmtf_dict::MMTFDict) = length(mmtf_dict.dict)
Base.iterate(mmtf_dict::MMTFDict) = iterate(mmtf_dict.dict)
Base.iterate(mmtf_dict::MMTFDict, i) = iterate(mmtf_dict.dict, i)

function Base.show(io::IO, mmtf_dict::MMTFDict)
    print(io, "MMTF dictionary with $(length(mmtf_dict)) fields")
end

function Base.read(input::IO,
            ::Type{MMTF};
            structure_name::AbstractString="",
            remove_disorder::Bool=false,
            read_std_atoms::Bool=true,
            read_het_atoms::Bool=true,
            gzip::Bool=false,
            run_dssp::Bool=false,
            run_stride::Bool=false,)
    d = MMTFDict(parsemmtf(input; gzip=gzip))
    ProteinStructure(d;
                     structure_name=structure_name,
                     remove_disorder=remove_disorder,
                     read_std_atoms=read_std_atoms,
                     read_het_atoms=read_het_atoms,
                     run_dssp=run_dssp,
                     run_stride=run_stride)
end

function ProteinStructure(d::MMTFDict;
            structure_name::AbstractString="",
            remove_disorder::Bool=false,
            read_std_atoms::Bool=true,
            read_het_atoms::Bool=true,
            run_dssp::Bool=false,
            run_stride::Bool=false)
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
            for gi in 1:d["groupsPerChain"][chain_i]
                group_i += 1
                # 0-based indexing in MMTF
                group = d["groupList"][d["groupTypeList"][group_i] + 1]
                for ai in 1:length(group["atomNameList"])
                    atom_i += 1
                    if (read_std_atoms || hets[chain_i]) && (read_het_atoms || !hets[chain_i])
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
    end
    if length(d["atomIdList"]) != atom_i
        throw(ArgumentError("Discrepancy between atom count ($(length(d["atomIdList"]))) and atoms read in ($atom_i)"))
    end
    # Remove any models that were not added to
    for model in struc
        if countchains(model) == 0
            delete!(models(struc), modelnumber(model))
        end
    end
    fixlists!(struc)

    # Run DSSP and STRIDE if required
    if run_dssp
        rundssp!(struc)
    end

    if run_stride
        runstride!(struc)
    end

    return struc
end

"""
    generatechainid(entity_id)

Convert a positive `Integer` into a chain ID.

Goes A to Z, then AA to ZA, AB to ZB etc.
This is in line with Protein Data Bank (PDB) conventions.
"""
function generatechainid(entity_id::Integer)
    entity_id > 0 || throw(ArgumentError("Entity ID $entity_id is not positive"))
    divisor = entity_id
    out_string = ""
    while divisor > 0
        mod = (divisor - 1) % 26
        out_string *= Char(65 + mod)
        divisor = Int(floor((divisor - mod) / 26))
    end
    return out_string
end

"""
    writemmtf(output, element, atom_selectors...; gzip=false)
    writemmtf(output, mmtf_dict; gzip=false)

Write a `StructuralElementOrList` or a `MMTFDict` to a MMTF file or output
stream.

Atom selector functions can be given as additional arguments - only atoms
that return `true` from all the functions are retained.
The keyword argument `expand_disordered` (default `true`) determines whether to
return all copies of disordered residues and atoms.
The keyword argument `gzip` (default `false`) determines if the file should be
gzipped.
"""
function writemmtf(output::Union{AbstractString, IO},
                d::MMTFDict;
                gzip::Bool=false)
    writemmtf(d.dict, output; gzip=gzip)
    return
end

function writemmtf(output::Union{AbstractString, IO},
                el::StructuralElementOrList,
                atom_selectors::Function...;
                expand_disordered::Bool=true,
                gzip::Bool=false)
    loop_el = isa(el, ProteinStructure) ? collectmodels(el) : el
    ats = sort(sort(sort(collectatoms(loop_el, atom_selectors...;
                            expand_disordered=expand_disordered)),
                            by=residue), by=model)
    d = MMTFDict()
    if length(ats) == 0
        return writemmtf(output, d; gzip=gzip)
    end
    last_model = model(first(ats))
    last_chain = chain(first(ats))
    last_res = residue(first(ats))
    chain_i = 0
    group_count = 0
    sequence = ""
    prev_resname = ""
    prev_het = true
    for (i, at) in enumerate(ats)
        if model(at) != last_model
            push!(d["chainsPerModel"], chain_i)
            chain_i = 0
            last_model = model(at)
        end

        if chain(at) != last_chain || i == 1
            # MMTF splits chains up by molecules so we determine chain splits
            #   at the residue level
            # The below just ensures we have a new entity for each chain
            prev_resname = ""
            prev_het = true
            last_chain = chain(at)
        end

        res = residue(at)
        if residue(at) != last_res || i == 1
            # Determine whether we have changed entity
            # ATOM blocks, and hetero molecules with the same name, are
            #   treated as the same entity
            if ishetero(res) != prev_het || (prev_het && resname(res) != prev_resname) || group_count == 0
                chain_i += 1
                push!(d["chainIdList"], generatechainid(chain_i))
                push!(d["chainNameList"], chainid(res))
                # Add the groupsPerChain and sequence for the previous chain
                if group_count > 0
                    push!(d["groupsPerChain"], group_count)
                    d["entityList"][end]["sequence"] = sequence
                    group_count = 0
                    sequence = ""
                end
                # Checking for similar entities is non-trivial so we treat
                #   each molecule as a separate entity
                push!(d["entityList"], Dict{Any, Any}(
                    "chainIndexList"=> Any[length(d["chainIdList"]) - 1],
                    "description"   => "",
                    "sequence"      => "", # This is changed later
                    "type"          => ishetero(res) ? "non-polymer" : "polymer"
                ))
            end
            if !ishetero(res)
                sequence *= string(LongAA(res; gaps=false))
            end
            group_count += 1

            # Look for an existing group with the correct residue name and
            #   atom names present
            # Look forwards in the loop to get atom names
            ats_res = AbstractAtom[at]
            for j in (i + 1):length(ats)
                if residue(ats[j]) == res
                    push!(ats_res, ats[j])
                else
                    break
                end
            end
            at_names = atomname.(ats_res)
            group_i = 0
            for (gi, group) in enumerate(d["groupList"])
                if group["groupName"] == resname(res) && group["atomNameList"] == at_names
                    group_i = gi
                    break
                end
            end

            if group_i == 0
                push!(d["groupList"], Dict{Any, Any}(
                    "groupName"       => resname(res),
                    "bondAtomList"    => Any[],
                    "elementList"     => Any[element(at) for at in ats_res],
                    # MMTF specifies missing charges as zero
                    "formalChargeList"=> Any[charge(at) == "" ? 0 : parse(Int64, charge(at)) for at in ats_res],
                    "singleLetterCode"=> "",
                    "chemCompType"    => "",
                    "atomNameList"    => Any[at_names...],
                    "bondOrderList"   => Any[]
                ))
            end

            push!(d["groupIdList"], resnumber(res))
            push!(d["groupTypeList"], group_i == 0 ? length(d["groupList"]) - 1 : group_i - 1)
            push!(d["insCodeList"], inscode(res) == ' ' ? '\0' : inscode(res))
            push!(d["secStructList"], -1)
            push!(d["sequenceIndexList"], ishetero(res) ? -1 : length(sequence) - 1)

            prev_resname = resname(res)
            prev_het = ishetero(res)
            last_res = res
        end

        push!(d["altLocList"], altlocid(at) == ' ' ? '\0' : altlocid(at))
        push!(d["atomIdList"], serial(at))
        push!(d["bFactorList"], tempfactor(at))
        push!(d["occupancyList"], occupancy(at))
        push!(d["xCoordList"], x(at))
        push!(d["yCoordList"], y(at))
        push!(d["zCoordList"], z(at))
    end

    push!(d["chainsPerModel"], chain_i)
    push!(d["groupsPerChain"], group_count)
    d["entityList"][end]["sequence"] = sequence

    d["numModels"] = length(d["chainsPerModel"])
    d["numChains"] = length(d["chainIdList"])
    d["numGroups"] = length(d["groupIdList"])
    d["numAtoms"] = length(d["atomIdList"])
    d["structureId"] = structurename(first(el))
    d["mmtfVersion"] = "1.0.0"
    d["mmtfProducer"] = "BioStructures.jl"

    return writemmtf(output, d; gzip=gzip)
end
