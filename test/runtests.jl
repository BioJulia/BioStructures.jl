module TestBioStructures

using Aqua
using BioAlignments
using BioSequences
import BioCore # Imported to avoid clash with BioGenerics distance
using CodecZlib
using DSSP_jll
using DataFrames
using Format
using Graphs
import MMTF # Imported to avoid clash with writemmtf
using MetaGraphs
using RecipesBase
using STRIDE_jll

using Downloads
using LinearAlgebra
using Test

using BioStructures
using BioStructures:
    x,
    x!,
    y,
    y!,
    z,
    z!,
    fixlists!,
    parseserial,
    parseatomname,
    parsealtloc,
    parseresname,
    parsechainid,
    parseresnumber,
    parseinscode,
    parsecoordx,
    parsecoordy,
    parsecoordz,
    parseoccupancy,
    parsetempfac,
    parseelement,
    parsecharge,
    spacestring,
    coordspec,
    floatspec,
    checkchainerror,
    splitline,
    tokenizecif,
    tokenizecifstructure,
    formatmmcifcol,
    requiresnewline,
    requiresquote,
    pdb_download_prefix

# Get the path to BioFmtSpecimens and download it if required
fmtdir = BioCore.Testing.get_bio_fmt_specimens("master", false)

# Access files in BioFmtSpecimens to test against
testfilepath(path::AbstractString...) = joinpath(fmtdir, path...)

# All writing is done to one temporary file and temporary directory which are removed at the end
const temp_filename = tempname()
const temp_dir = mktempdir()

# Gzip a file
function gzip_file(infile, outfile)
    open(infile, "r") do in
        open(outfile, "w") do out
            gz = GzipCompressorStream(out)
            write(gz, in)
            close(gz)
        end
    end
end

# countlines with optional gzip
function countlines_gzip(filename::AbstractString; gzip=false)
    if gzip
        open(GzipDecompressorStream, filename) do f
            return countlines(f)
        end
    else
        return countlines(filename)
    end
end

getparent(a::Atom) = a.residue
getparent(r::Residue) = r.chain
getparent(c::Chain) = c.model
getparent(m::Model) = m.structure
getchildren(::Atom) = ()
getchildren(r::Residue) = values(r.atoms)
getchildren(c::Chain) = values(c.residues)
getchildren(m::Model) = values(m.chains)
getchildren(s::MolecularStructure) = values(s.models)

function testparent(children, parent)
    for child in children
        if isa(child, DisorderedAtom)
            for a in values(child.alt_loc_ids)
                @test getparent(a) === parent
            end
        elseif isa(child, DisorderedResidue)
            for r in values(child.names)
                @test getparent(r) === parent
                testparent(getchildren(r), r)
            end
        else
            @test getparent(child) === parent
            testparent(getchildren(child), child)
        end
    end
end

Aqua.test_all(BioStructures; ambiguities=(recursive=false))

# This is the only test set that requires an internet connection
@testset "PDB interface" begin
    @test length(pdbentrylist()) > 100000

    # This may be empty on a given date so we just check it has the correct type
    @test isa(pdbstatuslist("$pdb_download_prefix/data/status/latest/added.pdb"), Vector{String})
    # Invalid URL
    @test_throws RequestError pdbstatuslist("$pdb_download_prefix/data/status/latest/dummy.pdb")

    addedlist, modifiedlist, obsoletelist = pdbrecentchanges()

    @test length(pdbobsoletelist()) > 3600

    # Invalid PDB ID format
    @test_throws ArgumentError downloadpdb("1a df")
    # Valid PDB ID format but PDB does not exist
    @test_throws RequestError downloadpdb("no1e", dir=temp_dir)
    # Invalid file format
    @test_throws TypeError downloadpdb("1alw", dir=temp_dir, format=String)
    # Biological assembly not available in PDBXML and MMTF
    @test_throws ArgumentError downloadpdb("1alw", dir=temp_dir, format=PDBXMLFormat, ba_number=1)
    # Invalid BA number for this PDB entry
    @test_throws RequestError downloadpdb("1alw", dir=temp_dir, format=MMCIFFormat, ba_number=10)
    # Test if downloadpdb returns the path to the downloaded file
    @test isfile(downloadpdb("1crn", dir=temp_dir))

    # PDB format
    downloadpdb("1alw", dir=temp_dir, format=PDBFormat)
    pdbpath = joinpath(temp_dir, "1ALW.$(pdbextension[PDBFormat])")
    @test isfile(pdbpath) && filesize(pdbpath) > 0
    # PDBXML format
    downloadpdb("1alw", dir=temp_dir, format=PDBXMLFormat)
    pdbpath = joinpath(temp_dir, "1ALW.$(pdbextension[PDBXMLFormat])")
    @test isfile(pdbpath) && filesize(pdbpath) > 0
    # mmCIF format
    downloadpdb("1alw", dir=temp_dir, format=MMCIFFormat)
    pdbpath = joinpath(temp_dir, "1ALW.$(pdbextension[MMCIFFormat])")
    @test isfile(pdbpath) && filesize(pdbpath) > 0
    @test isfile(pdbpath) && filesize(pdbpath) > 0
    # Obsolete PDB
    downloadpdb("116l", dir=temp_dir, format=PDBFormat, obsolete=true)
    pdbpath = joinpath(temp_dir, "obsolete", "116L.$(pdbextension[PDBFormat])")
    @test isfile(pdbpath) && filesize(pdbpath) > 0
    # Biological assembly - PDB format
    downloadpdb("1alw", dir=temp_dir, format=PDBFormat, ba_number=1)
    pdbpath = joinpath(temp_dir, "1ALW_ba1.$(pdbextension[PDBFormat])")
    @test isfile(pdbpath) && filesize(pdbpath) > 0
    # Biological assembly - mmCIF format
    downloadpdb("5a9z", dir=temp_dir, format=MMCIFFormat, ba_number=1)
    pdbpath = joinpath(temp_dir, "5A9Z_ba1.$(pdbextension[MMCIFFormat])")
    @test isfile(pdbpath) && filesize(pdbpath) > 0
    # Download multiple PDB files
    pdbidlist = ["1ent", "1en2"]
    downloadpdb(pdbidlist, dir=temp_dir, format=PDBFormat)
    for pdbid in pdbidlist
        pdbpath = joinpath(temp_dir, "$(uppercase(pdbid)).$(pdbextension[PDBFormat])")
        @test isfile(pdbpath) && filesize(pdbpath) > 0
    end

    # Test function as first argument to downloadpdb
    @test downloadpdb(countlines, "1alw") == 3584
    @test downloadpdb("1alw") do fp
        s = read(fp, PDBFormat)
        return countresidues(s, standardselector)
    end == 346
    @test downloadpdb(["169l", "2lzm"]) do fps
        n_chains = Int[]
        for fp in fps
            s = read(fp, PDBFormat)
            push!(n_chains, countchains(s))
        end
        return n_chains
    end == [5, 1]

    # Test retrievepdb
    struc = retrievepdb("1AKE", dir=temp_dir, structure_name="New name")
    @test structurename(struc) == "New name"
    @test countatoms(struc) == 3804

    struc = retrievepdb("1AKE", dir=temp_dir, obsolete=true, read_het_atoms=false)
    @test countatoms(struc) == 3312
    @test serial(collectatoms(struc)[2000]) == 2005
    @test sum(ishetero, collectatoms(struc)) == 0

    struc = retrievepdb("1AKE", dir=temp_dir, ba_number=1, read_het_atoms=false, read_std_atoms=false)
    @test countatoms(struc) == 0
    @test countresidues(struc) == 0
    @test countchains(struc) == 0
    @test countmodels(struc) == 0
end

@testset "Types" begin
    # Test constructors and indexing
    struc = MolecularStructure("Test structure")
    struc[1] = Model(1, struc)
    mo = struc[1]
    @test isa(mo, Model)
    struc[3] = Model(3, struc)
    struc[1]['A'] = Chain('A', mo)
    ch = struc[1]['A']
    @test isa(ch, Chain)
    struc['B'] = Chain('B', mo)
    struc['A'][10] = Residue("ALA", 10, ' ', false, ch)
    res = struc['A'][10]
    @test isa(res, Residue)
    struc['A']["H_20A"] = DisorderedResidue(Dict(
        " VA" => Residue(" VA", 20, 'A', true, ch),
        "ILE" => Residue("ILE", 20, 'A', true, ch)
    ), " VA")
    dis_res = struc['A']["H_20A"]
    @test isa(dis_res, DisorderedResidue)
    a = Atom(100, " CA ", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res)
    struc['A'][10][" CA "] = a
    at = struc['A'][10]["CA"]
    @test isa(at, Atom)
    struc['A'][10][" CB "] = DisorderedAtom(Dict(
        'A' => Atom(200, " CB ", 'A', [10.0, 20.0, 30.0], 0.6, 20.0, " C", "1+", res),
        'B' => Atom(201, " CB ", 'B', [11.0, 21.0, 31.0], 0.4, 30.0, " C", "1+", res)
    ), 'A')
    dis_at = struc['A'][10]["CB"]
    @test isa(dis_at, DisorderedAtom)
    struc['A']["H_20A"][" CG "] = Atom(
        300, " CG ", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", defaultresidue(dis_res))
    disorderedres(dis_res, "ILE")[" O  "] = Atom(
        400, " O  ", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " O", "  ", disorderedres(dis_res, "ILE"))
    fixlists!(struc)

    # copy doesn't share memory
    testparent(getchildren(struc), struc)
    struc_copy = copy(struc)
    testparent(getchildren(struc_copy), struc_copy)
    struc_copy['A'][10]["CA"].coords[2] = 100
    @test struc_copy['A'][10]["CA"].coords[2] == 100
    @test a.coords[2] == 2
    @test struc['A'][10]["CA"].coords[2] == 2
    # Intermediate copies preserve parenting up to the node of the copy
    testparent(getchildren(mo), mo)
    mo_copy = copy(mo)
    testparent(getchildren(mo_copy), mo_copy)
    testparent(getchildren(ch), ch)
    ch_copy = copy(ch)
    testparent(getchildren(ch_copy), ch_copy)
    testparent(getchildren(res), res)
    res_copy = copy(res)
    testparent(getchildren(res_copy), res_copy)
    @test copy(dis_res) isa DisorderedResidue
    @test copy(dis_at) isa DisorderedAtom

    # Test alternate constructors
    MolecularStructure("struc", Dict(1 => Model()))
    MolecularStructure()
    mmcif_dict = MMCIFDict(testfilepath("mmCIF", "1AKE.cif"))
    struc_mmcif_1ake = MolecularStructure(mmcif_dict)
    @test countatoms(struc_mmcif_1ake) == 3804
    mmtf_dict = MMTFDict(testfilepath("MMTF", "1AKE.mmtf"))
    struc_mmtf_1ake = MolecularStructure(mmtf_dict)
    @test countatoms(struc_mmtf_1ake) == 3804

    Model(1, Dict("A" => Chain('A')), MolecularStructure())
    Model(1, MolecularStructure())
    Model(1)
    Model()

    Chain("A", ["1"], Dict("1" => Residue("ALA", 1, ' ', false, Chain('A'))), Model())
    Chain("A", Model())
    Chain('A', Model())
    Chain("A")
    Chain('A')

    Residue("ALA", 1, ' ', false, ["CA"], Dict("CA" =>
        Atom(1, "CA", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "  ", "  ",
        Residue("ALA", 1, ' ', false, Chain('A')))), Chain('A'), '-')
    Residue("ALA", 1, ' ', false, Chain('A'))

    # Test show
    show(devnull, at)
    show(devnull, dis_at)
    show(devnull, res)
    show(devnull, dis_res)
    show(devnull, ch)
    show(devnull, mo)
    show(devnull, struc)
    show(devnull, Model())
    show(devnull, MolecularStructure())

    # Test getters/setters
    @test serial(at) == 100
    @test serial(dis_at) == 200

    @test atomname(at) == "CA"
    @test atomname(dis_at) == "CB"
    @test atomname(at, strip=false) == " CA "
    @test atomname(dis_at, strip=false) == " CB "

    @test altlocid(at) == ' '
    @test altlocid(dis_at) == 'A'

    @test x(at) == 1.0
    @test x(dis_at) == 10.0

    x!(at, 5.0)
    @test coords(at) == [5.0, 2.0, 3.0]
    # Only the coordinates on the default atom are changed
    x!(dis_at, 12.0)
    @test coords(dis_at) == [12.0, 20.0, 30.0]
    @test coords(dis_at['B']) == [11.0, 21.0, 31.0]

    @test y(at) == 2.0
    @test y(dis_at) == 20.0

    y!(at, 6.0)
    @test coords(at) == [5.0, 6.0, 3.0]
    y!(dis_at, 22.0)
    @test coords(dis_at) == [12.0, 22.0, 30.0]
    @test coords(dis_at['B']) == [11.0, 21.0, 31.0]

    @test z(at) == 3.0
    @test z(dis_at) == 30.0

    z!(at, 7.0)
    @test coords(at) == [5.0, 6.0, 7.0]
    z!(dis_at, 32.0)
    @test coords(dis_at) == [12.0, 22.0, 32.0]
    @test coords(dis_at['B']) == [11.0, 21.0, 31.0]

    coords!(at, [40.0, 50.0, 60.0])
    @test coords(at) == [40.0, 50.0, 60.0]
    coords!(dis_at, [40.0, 50.0, 60.0])
    @test coords(dis_at) == [40.0, 50.0, 60.0]
    @test coords(dis_at['B']) == [11.0, 21.0, 31.0]
    @test_throws ArgumentError coords!(at, [40.0, 50.0, 60.0, 70.0])
    x!(dis_at['A'], 100.0)
    @test coords(dis_at) == [100.0, 50.0, 60.0]
    @test coords(dis_at['B']) == [11.0, 21.0, 31.0]
    x!(dis_at['B'], 110.0)
    @test coords(dis_at) == [100.0, 50.0, 60.0]
    @test coords(dis_at['B']) == [110.0, 21.0, 31.0]

    @test occupancy(at) == 1.0
    @test occupancy(dis_at) == 0.6

    @test tempfactor(at) == 10.0
    @test tempfactor(dis_at) == 20.0

    @test element(at) == "C"
    @test element(dis_at) == "C"
    @test element(at, strip=false) == " C"
    @test element(dis_at, strip=false) == " C"

    @test charge(at) == ""
    @test charge(dis_at) == "1+"
    @test charge(at, strip=false) == "  "
    @test charge(dis_at, strip=false) == "1+"

    @test residue(at) == res
    @test residue(dis_at) == res
    @test residue(res) == res
    @test residue(dis_res) == dis_res

    @test !ishetero(res)
    @test !ishetero(at)
    @test !ishetero(dis_at)
    @test ishetero(dis_res)

    @test !isdisorderedatom(at)
    @test isdisorderedatom(dis_at)

    @test defaultaltlocid(dis_at) == 'A'

    dis_at_mod = DisorderedAtom(dis_at, 'B')
    @test defaultaltlocid(dis_at_mod) == 'B'
    @test serial(dis_at_mod) == 201
    @test_throws ArgumentError DisorderedAtom(dis_at, 'C')

    @test isa(defaultatom(dis_at), Atom)
    @test serial(defaultatom(dis_at)) == 200

    @test altlocids(at) == [' ']
    @test altlocids(dis_at) == ['A', 'B']

    @test atomid(at) == ("10:A", "ALA", "CA")
    @test atomid(dis_at) == ("10:A", "ALA", "CB")

    @test resname(res) == "ALA"
    @test resname(at) == "ALA"
    @test resname(dis_at) == "ALA"
    @test resname(dis_res) == "VA"
    @test resname(dis_res, strip=false) == " VA"

    @test resnumber(res) == 10
    @test resnumber(at) == 10
    @test resnumber(dis_at) == 10
    @test resnumber(dis_res) == 20

    @test inscode(res) == ' '
    @test inscode(at) == ' '
    @test inscode(dis_at) == ' '
    @test inscode(dis_res) == 'A'

    @test resid(at) == "10"
    @test resid(dis_at) == "10"
    @test resid(at, full=true) == "10:A"
    @test resid(dis_at, full=true) == "10:A"
    @test resid(res) == "10"
    @test resid(res, full=true) == "10:A"
    @test resid(dis_res) == "H_20A"
    @test resid(dis_res, full=true) == "H_20A:A"

    @test atomnames(res) == ["CA", "CB"]
    @test atomnames(dis_res) == ["CG"]
    @test atomnames(res, strip=false) == [" CA ", " CB "]
    @test atomnames(dis_res, strip=false) == [" CG "]

    @test isa(atoms(res), Dict{String, AbstractAtom})
    @test length(atoms(res)) == 2
    @test serial(atoms(res)[" CA "]) == 100
    @test isa(atoms(dis_res), Dict{String, AbstractAtom})
    @test length(atoms(dis_res)) == 1
    @test serial(atoms(dis_res)[" CG "]) == 300

    @test !isdisorderedres(res)
    @test isdisorderedres(dis_res)

    @test isa(disorderedres(dis_res, "ILE"), AbstractResidue)
    @test resname(disorderedres(dis_res, "ILE")) == "ILE"

    @test defaultresname(dis_res) == " VA"

    @test isa(defaultresidue(dis_res), Residue)
    @test resname(defaultresidue(dis_res)) == "VA"

    @test resnames(res) == ["ALA"]
    @test resnames(dis_res) == [" VA", "ILE"]

    dis_res_mod = DisorderedResidue(dis_res, "ILE")
    @test defaultresname(dis_res_mod) == "ILE"
    @test atomnames(dis_res_mod) == ["O"]
    @test_throws ArgumentError DisorderedResidue(dis_res, "SER")

    @test chain(at) == ch
    @test chain(dis_at) == ch
    @test chain(res) == ch
    @test chain(dis_res) == ch
    @test chain(ch) == ch

    @test chainid(at) == "A"
    @test chainid(dis_at) == "A"
    @test chainid(res) == "A"
    @test chainid(dis_res) == "A"
    @test chainid(ch) == "A"

    # Test modifying the chain ID
    chainid!(ch, "C")
    @test chainid(at) == "C"
    @test chainid(dis_at) == "C"
    @test chainid(res) == "C"
    @test chainid(dis_res) == "C"
    @test chainid(ch) == "C"
    @test mo["C"] == ch
    @test_throws KeyError mo["A"]

    @test chainids(mo) == ["B", "C"]
    chainid!(ch, "A")

    # Move one of the residues to a new chain and StructuralElements below it
    #   should identify with the new chain
    chainid!(res, "C")
    @test chainid(res) == "C"
    @test chainid(at) == "C"
    @test chainid(dis_at) == "C"
    @test chainid(dis_res) == "A"
    @test chainid(ch) == "A"

    @test chainids(mo) == ["A", "B", "C"]

    # Emptying a chain of residues by moving its residues deletes the chain
    chainid!(res, "A")

    @test chainids(mo) == ["A", "B"]

    # Reassigning a chain ID to one that already exists errors
    @test_throws PDBConsistencyError chainid!(ch, "B")
    struc['B'][10] = Residue("ALA", 10, ' ', false, struc['B'])
    # Reassigning a residue with a number to a chain that already has one of that number errors
    @test_throws PDBConsistencyError chainid!(struc['B'][10], "A")

    error = PDBConsistencyError("message")
    showerror(devnull, error)

    @test resids(ch) == ["10", "H_20A"]

    @test isa(residues(ch), Dict{String, AbstractResidue})
    @test length(residues(ch)) == 2
    @test serial(residues(ch)["10"]["CA"]) == 100

    @test model(at) == mo
    @test model(dis_at) == mo
    @test model(res) == mo
    @test model(dis_res) == mo
    @test model(ch) == mo
    @test model(mo) == mo

    @test modelnumber(at) == 1
    @test modelnumber(dis_at) == 1
    @test modelnumber(res) == 1
    @test modelnumber(dis_res) == 1
    @test modelnumber(ch) == 1
    @test modelnumber(mo) == 1

    @test chainids(mo) == ["A", "B"]
    @test chainids(struc) == ["A", "B"]
    @test chainids(MolecularStructure()) == String[]
    @test chainids(Model()) == String[]

    @test isa(chains(mo), Dict{String, Chain})
    @test length(chains(mo)) == 2
    @test resname(chains(mo)["A"]["H_20A"]) == "VA"
    @test isa(chains(Model()), Dict{String, Chain})
    @test isa(chains(struc), Dict{String, Chain})
    @test length(chains(struc)) == 2
    @test isa(chains(MolecularStructure()), Dict{String, Chain})

    @test structure(at) == struc
    @test structure(dis_at) == struc
    @test structure(res) == struc
    @test structure(dis_res) == struc
    @test structure(ch) == struc
    @test structure(mo) == struc
    @test structure(struc) == struc

    @test structurename(at) == "Test structure"
    @test structurename(dis_at) == "Test structure"
    @test structurename(res) == "Test structure"
    @test structurename(dis_res) == "Test structure"
    @test structurename(ch) == "Test structure"
    @test structurename(mo) == "Test structure"
    @test structurename(struc) == "Test structure"

    @test modelnumbers(struc) == [1, 3]

    @test isa(models(struc), Dict{Int, Model})
    @test length(models(struc)) == 2
    @test resids(models(struc)[1]['A']) == ["10", "H_20A"]

    @test isa(defaultmodel(struc), Model)
    @test modelnumber(defaultmodel(struc)) == 1

    # Test iteration and collecting of elements
    for at_col in (collect(at), eltype(at)[i for i in at])
        @test isa(at_col, Vector{Atom})
        @test length(at_col) == 1
        @test serial(at_col[1]) == 100
    end

    for dis_at_col in (collect(dis_at), eltype(dis_at)[i for i in dis_at])
        @test isa(dis_at_col, Vector{Atom})
        @test length(dis_at_col) == 2
        @test serial(dis_at_col[2]) == 201
    end

    for res_col in (collect(res), eltype(res)[i for i in res])
        @test isa(res_col, Vector{AbstractAtom})
        @test length(res_col) == 2
        @test serial(res_col[2]) == 200
    end

    for dis_res_col in (collect(dis_res), eltype(dis_res)[i for i in dis_res])
        @test isa(dis_res_col, Vector{AbstractAtom})
        @test length(dis_res_col) == 1
        @test serial(dis_res_col[1]) == 300
    end

    for ch_col in (collect(ch), eltype(ch)[i for i in ch])
        @test isa(ch_col, Vector{AbstractResidue})
        @test length(ch_col) == 2
        @test resname(ch_col[2]) == "VA"
    end

    for mo_col in (collect(mo), eltype(mo)[i for i in mo])
        @test isa(mo_col, Vector{Chain})
        @test length(mo_col) == 2
        @test chainid(mo_col[2]) == "B"
    end

    for struc_col in (collect(struc), eltype(struc)[i for i in struc])
        @test isa(struc_col, Vector{Model})
        @test length(struc_col) == 2
        @test modelnumber(struc_col[2]) == 3
    end

    # Test element indices
    @test isa(dis_at['A'], Atom)
    @test serial(dis_at['A']) == 200
    @test serial(dis_at['B']) == 201
    @test_throws KeyError dis_at['C']
    @test_throws MethodError dis_at["A"]

    @test isa(res["CA"], AbstractAtom)
    @test serial(res["CA"]) == 100
    @test serial(res["CB"]) == 200
    @test_throws KeyError res["N"]

    @test isa(dis_res["CG"], AbstractAtom)
    @test serial(dis_res["CG"]) == 300
    @test_throws KeyError dis_res["N"]

    @test isa(ch["10"], AbstractResidue)
    @test isa(ch[10], AbstractResidue)
    @test_throws KeyError ch["H_10"]
    @test serial(ch["10"]["CA"]) == 100
    @test serial(ch[10]["CA"]) == 100
    @test resname(ch["H_20A"]) == "VA"
    @test_throws KeyError ch["11"]
    @test_throws KeyError ch[11]

    @test isa(mo["A"], Chain)
    @test isa(mo['A'], Chain)
    @test resname(mo['A']["H_20A"]) == "VA"
    @test_throws KeyError mo['C']

    @test isa(struc[1], Model)
    @test isa(struc[3], Model)
    @test isa(struc["A"], Chain)
    @test isa(struc['A'], Chain)
    @test_throws KeyError struc[2]
    @test_throws KeyError struc['C']
    @test ishetero(struc[1]['A']["H_20A"])

    # Test selector functions
    ch_a = Chain('A')
    ch_a["10"] = Residue("ALA", 10, ' ', false, ch_a)
    res_a = ch_a["10"]
    ch_a["H_11"] = Residue("MG", 11, ' ', true, ch_a)
    res_b = ch_a["H_11"]
    ch_a["H_100"] = Residue("HOH", 100, ' ', true, ch_a)
    res_c = ch_a["H_100"]
    ch_a["10"]["CA"] = Atom(
        100, " CA ", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res_a)
    at_a = ch_a["10"]["CA"]
    ch_a["H_11"]["MG"] = Atom(
        110, "MG", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res_b)
    at_b = ch_a["H_11"]["MG"]

    @test standardselector(at_a)
    @test !standardselector(at_b)
    @test standardselector(res_a)
    @test !standardselector(res_b)
    @test !heteroselector(at_a)
    @test heteroselector(at_b)
    @test !heteroselector(res_a)
    @test heteroselector(res_b)
    @test atomnameselector(at_a, Set(["CA", "N", "C"]))
    @test atomnameselector(at_a, ["CA", "N", "C"])
    @test !atomnameselector(at_b, Set(["CA", "N", "C"]))
    @test !atomnameselector(at_a, ["CA", "N", "C"], strip=false)
    @test calphaselector(at_a)
    @test !calphaselector(at_b)
    @test !calphaselector(Atom(
        100, "CA", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res_b))
    @test !cbetaselector(at_a)
    @test cbetaselector(Atom(
        100, "CB", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " H", "  ", res_a))
    @test backboneselector(at_a)
    @test !backboneselector(at_b)
    @test !backboneselector(Atom(
        100, "CA", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res_b))
    @test heavyatomselector(at_a)
    @test !heavyatomselector(at_b)
    @test !heavyatomselector(Atom(
        100, "H1", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " H", "  ", res_a))
    @test resnameselector(at_a, Set(["ALA"]))
    @test resnameselector(at_a, ["ALA"])
    @test !resnameselector(at_b, Set(["ALA"]))
    @test resnameselector(res_a, Set(["ALA"]))
    @test !resnameselector(res_b, Set(["ALA"]))
    @test !waterselector(res_a)
    @test waterselector(res_c)
    @test !waterselector(at_a)
    @test notwaterselector(res_a)
    @test !notwaterselector(res_c)
    @test notwaterselector(at_a)
    @test !disorderselector(at_a)
    @test disorderselector(dis_at)
    @test !disorderselector(res_a)
    @test disorderselector(dis_res)
    @test hydrogenselector(Atom(
        100, "H", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " H", "  ", res_a))
    @test !hydrogenselector(Atom(
        100, "H", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res_a))
    @test hydrogenselector(Atom(
        100, "H1", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "  ", "  ", res_a))
    @test hydrogenselector(Atom(
        100, "1H", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "  ", "  ", res_a))
    @test !hydrogenselector(Atom(
        100, "NH1", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "  ", "  ", res_a))
    @test allselector(at_a)
    @test allselector(res_a)

    # Further tests for structural element ordering
    # Order when looping over a DisorderedAtom is the atom serial
    res = Residue("ALA", 1, ' ', false, Chain('A'))
    dis_at_ord = DisorderedAtom(Dict(
        'A' => Atom(102, "CA", 'A', [1.0, 2.0, 3.0], 0.3, 10.0, "C", "", res),
        'B' => Atom(101, "CA", 'B', [1.0, 2.0, 3.0], 0.4, 10.0, "C", "", res),
        'C' => Atom(100, "CA", 'C', [1.0, 2.0, 3.0], 0.3, 10.0, "C", "", res),
    ), 'B')
    @test altlocids(dis_at_ord) == ['C', 'B', 'A']

    # Order when sorting an atom list is the atom serial
    at_list_ord = AbstractAtom[
        DisorderedAtom(Dict(
            'A' => Atom(100, "CA", 'A', [1.0, 2.0, 3.0], 0.4, 10.0, "C", "", res),
            'B' => Atom(104, "CA", 'B', [1.0, 2.0, 3.0], 0.6, 10.0, "C", "", res)
        ), 'B'),
        Atom(102, "CB", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", "", res),
        Atom(103, "CG", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, "C", "", res)
    ]
    @test atomname.(sort(at_list_ord)) == ["CB", "CG", "CA"]
    sort!(at_list_ord)
    @test atomname.(at_list_ord) == ["CB", "CG", "CA"]

    # Order when sorting a residue list is chain ID, then residue number,
    #   then insertion code, then stdres/hetres
    res_ord = AbstractResidue[
        Residue("ALA", 201, 'A', false, Chain('A')),
        Residue("ALA", 203, ' ', false, Chain('A')),
        Residue("ALA", 200, ' ', true,  Chain('A')),
        Residue("ALA", 201, 'B', false, Chain('A')),
        Residue("ALA", 202, ' ', false, Chain('A')),
        Residue("ALA", 300, ' ', false, Chain('B')),
        Residue("ALA", 201, ' ', true,  Chain('A')),
        Residue("ALA", 201, ' ', false, Chain('A')),
        Residue("ALA", 201, 'A', true,  Chain('A')),
        Residue("ALA", 100, ' ', false, Chain('B')),
        Residue("ALA", 203, ' ', true,  Chain('A')),
        Residue("ALA", 200, ' ', false, Chain('A')),
    ]

    # Test sequentialresidues
    @test sequentialresidues(res_ord[3], res_ord[7])
    @test !sequentialresidues(res_ord[7], res_ord[3])
    @test sequentialresidues(res_ord[7], res_ord[9])

    @test resid.(res_ord) == [
        "201A", "203", "H_200", "201B", "202", "300", "H_201", "201", "H_201A",
        "100", "H_203", "200"]
    @test resid.(res_ord; full=true) == [
        "201A:A", "203:A", "H_200:A", "201B:A", "202:A", "300:B", "H_201:A",
        "201:A", "H_201A:A", "100:B", "H_203:A", "200:A"]
    @test resid.(sort(res_ord); full=true) == [
        "200:A", "H_200:A", "201:A", "H_201:A", "201A:A", "H_201A:A", "201B:A",
        "202:A", "203:A", "H_203:A", "100:B", "300:B"]
    sort!(res_ord)
    @test resid.(res_ord; full=true) == [
        "200:A", "H_200:A", "201:A", "H_201:A", "201A:A", "H_201A:A", "201B:A",
        "202:A", "203:A", "H_203:A", "100:B", "300:B"]

    # Order of listing residue names in a DisorderedResidue is default then alphabetical
    dis_res_ord = DisorderedResidue(Dict(
        "THR" => Residue("THR", 201, ' ', false, Chain('A')),
        "ALA" => Residue("ALA", 201, ' ', false, Chain('A')),
        "ILE" => Residue("ILE", 201, ' ', false, Chain('A')),
        "SER" => Residue("SER", 201, ' ', false, Chain('A')),
        "VAL" => Residue("VAL", 201, ' ', false, Chain('A'))
    ), "SER")
    @test defaultresname(dis_res_ord) == "SER"
    @test resnames(dis_res_ord) == ["SER", "ALA", "ILE", "THR", "VAL"]

    # Order when sorting chain IDs is character ordering with the empty chain ID at the end
    mo_ord = Model(1, Dict(
        "AB"  => Chain("AB"),
        "A"   => Chain("A"),
        " "   => Chain(" "),
        "1"   => Chain("1"),
        "AAA" => Chain("AAA"),
        "BC"  => Chain("BC"),
        "a"   => Chain("a"),
        "X"   => Chain("X"),
        "AC"  => Chain("AC"),
    ), MolecularStructure())
    @test chainids(mo_ord) == ["1", "A", "X", "a", "AB", "AC", "BC", "AAA", " "]

    # Test sequence extraction
    @test threeletter_to_aa["ALA"] == AA_A
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat)
    seq = LongAA(struc['B'])
    @test seq == LongAA(
        "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNG" *
        "FLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQ" *
        "EETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILGX-------------------------" *
        "--------------------------------------------------------------------------------" *
        "--------------------------------------------------------------------------------" *
        "--------------------------------------------------------------------------------" *
        "---------------------------------------------------XX---------------------------" *
        "----------------------------------------XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" *
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" *
        "XXXXXXXXXXXXXXX"
    )
    seq = LongAA(struc['B'], gaps=false)
    @test seq == LongAA(
        "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNG" *
        "FLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQ" *
        "EETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILGXXXXXXXXXXXXXXXXXXXXXXXXXX" *
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" *
        "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    )
    seq = LongAA(struc['B'], standardselector)
    @test seq == LongAA(
        "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNG" *
        "FLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQ" *
        "EETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG"
    )
    seq = LongAA(struc, standardselector)
    @test seq == LongAA(
        "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNG" *
        "FLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQ" *
        "EETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG" *
        "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNG" *
        "FLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQ" *
        "EETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG"
    )
    seq = LongAA(AbstractResidue[
        Residue("VAL", 20, 'B', true, Chain('B')),
        Residue("ALA", 10, 'A', false, Chain('A')),
    ])
    @test seq == LongAA("VA")

    # Test pairwise alignment
    res = collectresidues(struc["A"], standardselector)
    alres = pairalign(res[1:40], res[21:60])
    al = alignment(alres)
    @test score(alres) == 40
    @test count_matches(al) == 20
    @test count_insertions(al) == 20
    @test length(al) == 60
    alres = pairalign(res[1:40], res[21:60], aligntype=LocalAlignment())
    al = alignment(alres)
    @test score(alres) == 100
    @test count_matches(al) == 20
    @test count_insertions(al) == 0
    @test length(al) == 20
    alres = pairalign(res[1:40], res[21:60],
            scoremodel=AffineGapScoreModel(BLOSUM62, gap_open=-50, gap_extend=-1))
    al = alignment(alres)
    @test score(alres) == -29
    @test count_matches(al) == 5
    @test count_insertions(al) == 0
    @test length(al) == 40

    # Test DataFrame constructor
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat)
    df = DataFrame(collectatoms(struc))
    @test size(df) == (3816, 17)
    @test first(df)[:x] == 26.981
    @test size(describe(df), 1) == 17
    @test isapprox(sum(df.tempfactor), 165121.99)
    @test first(sort(df, :x))[:x] == -7.668
    df = DataFrame(collectatoms(struc), expand_disordered=false)
    @test size(df) == (3804, 17)
    df = DataFrame(collectresidues(struc))
    @test size(df) == (808, 8)
    @test first(df)[:resname] == "MET"
    @test size(describe(df), 1) == 8
    @test sum(df.countatoms) == 3816
    @test first(sort(df, :resnumber, rev=true))[:resnumber] == 735
    df = DataFrame(collectresidues(struc), expand_disordered=false)
    @test size(df) == (808, 8)
    @test sum(df.countatoms) == 3804
    struc = read(testfilepath("PDB", "1EN2.pdb"), PDBFormat)
    # Residues need expanding first or you miss disordered residues
    df = DataFrame(collectatoms(collectresidues(struc; expand_disordered=true)))
    @test size(df) == (819, 17)
    df = DataFrame(collectatoms(struc), expand_disordered=false)
    @test size(df) == (754, 17)
    df = DataFrame(collectresidues(struc))
    @test size(df) == (171, 8)
    df = DataFrame(collectresidues(struc), expand_disordered=false)
    @test size(df) == (166, 8)
end

@testset "Selection syntax" begin
    struc = retrievepdb("4YC6", dir=temp_dir, run_dssp=true)
    @test length(collectatoms(struc, sel"all")) == 12271
    @test length(collectatoms(struc, sel"name CA")) == 1420
    @test length(collectatoms(struc, sel"name CA", sel"x > 0")) == 312
    @test length(collectatoms(struc, calphaselector, sel"x > 0")) == 312
    sel = collectatoms(struc, sel"index = 13")
    @test length(sel) == 1
    @test serial(sel[1]) == 13
    @test length(collectatoms(struc, sel"index > 1 and index < 13")) == 11
    @test length(collectatoms(struc, sel"index >= 1 and index <= 13")) == 13
    @test length(collectatoms(struc, sel"beta > 100.0")) == 4807
    @test length(collectatoms(struc, sel"protein")) == 11632
    @test length(collectatoms(struc, sel"water")) == 639
    @test length(collectatoms(struc, sel"resname GLY")) == 320
    @test length(collectatoms(struc, sel"protein and resnum = 2")) == 36
    @test length(collectatoms(struc, sel"neutral")) == 8300
    @test length(collectatoms(struc, sel"charged")) == 3332
    @test length(collectatoms(struc, sel"sidechain")) == 5952
    @test length(collectatoms(struc, sel"acidic")) == 1604
    @test length(collectatoms(struc, sel"basic")) == 1728
    @test length(collectatoms(struc, sel"hydrophobic")) == 3604
    @test length(collectatoms(struc, sel"not hydrophobic")) == 8667
    @test length(collectatoms(struc, sel"aliphatic")) == 3276
    @test length(collectatoms(struc, sel"aromatic")) == 2340
    @test length(collectatoms(struc, sel"polar")) == 6544
    @test length(collectatoms(struc, sel"nonpolar")) == 5088
    @test length(collectatoms(struc, sel"backbone")) == 5680
    @test length(collectatoms(struc, sel"element H")) == 0
    @test length(collectatoms(struc, sel"name CA or element S")) == 1464
    @test length(collectatoms(struc, sel"disordered")) == 68
    @test length(collectatoms(struc, sel"sscode E")) == 2448
    @test length(collectatoms(struc, sel"helix")) == 4047
    # Check interpolation support
    ss_type = "helix"
    @test length(collectatoms(struc, sel"$ss_type")) == 4047
    sel_chains = ('A', 'B')
    @test length(collectresidues(struc, sel"chain $(first(sel_chains)) or chain $(last(sel_chains))")) == 544

    @test length(collectresidues(struc, sel"chain A or chain B")) == 544
    @test length(collectresidues(struc, sel"standard")) == 1420
    @test length(collectresidues(struc, sel"coil")) == 642
    @test length(collectchains(struc, sel"chain B")) == 1
    @test length(collectchains(struc, sel"chain I")) == 0
    @test length(collectmodels(struc, sel"model 1")) == 1
    @test length(collectmodels(struc, sel"model 2")) == 0

    @test_throws ArgumentError collectatoms(struc, BioStructures.Select("abc")) # Invalid selection syntax
    @test_throws ArgumentError collectatoms(struc, BioStructures.Select("index = A")) # Invalid value type
    @test_throws ArgumentError collectatoms(struc, BioStructures.Select("resnum C"))

    # Test show method for @sel_str
    buff = IOBuffer()
    show(buff, MIME"text/plain"(), sel"name CA and resnum 1")
    @test String(take!(buff)) == """Select("name CA and resnum 1")"""
end

@testset "PDB reading" begin
    # Test parsing functions
    line = "ATOM    591  C   GLY A  80      29.876  54.131  35.806  1.00 40.97           C1+"
    @test parseserial(line) == 591
    @test parseatomname(line) == " C  "
    @test parsealtloc(line) == ' '
    @test parseresname(line) == "GLY"
    @test parsechainid(line) == "A"
    @test parseresnumber(line) == 80
    @test parseinscode(line) == ' '
    @test parsecoordx(line) == 29.876
    @test parsecoordy(line) == 54.131
    @test parsecoordz(line) == 35.806
    @test parseoccupancy(line) == 1.0
    @test parsetempfac(line) == 40.97
    @test parseelement(line) == " C"
    @test parsecharge(line) == "1+"

    line_short = "ATOM    591  C"
    @test_throws PDBParseError("line too short", 37, line_short) AtomRecord(line_short, 37)
    @test_throws PDBParseError    parseserial("ATOM         C   GLY A  80      29.876  54.131  35.806  1.00 40.97           C1+")
    @test_throws PDBParseError parseresnumber("ATOM    591  C   GLY A          29.876  54.131  35.806  1.00 40.97           C1+")
    @test_throws PDBParseError    parsecoordx("ATOM    591  C   GLY A  80      xxxxxx  54.131  35.806  1.00 40.97           C1+")
    @test_throws PDBParseError    parsecoordy("ATOM    591  C   GLY A  80      29.876  xxxxxx  35.806  1.00 40.97           C1+")
    @test_throws PDBParseError    parsecoordz("ATOM    591  C   GLY A  80      29.876  54.131  xxxxxx  1.00 40.97           C1+")
    line_medium = "ATOM    591  C   GLY A  80      29.876  54.131  35.806"
    rec = AtomRecord(line_medium, 55)
    @test rec.occupancy == 1.0
    @test rec.temp_factor == 0.0
    @test rec.element == "  "
    @test rec.charge == "  "

    # Test AtomRecord constructor
    line_a = "ATOM    669  CA  ILE A  90      31.743  33.110  31.221  1.00 25.76           C  "
    line_b = "HETATM 3474  O  B XX A 334A      8.802  62.000   8.672  1.00 39.15           O1-"
    line_c = "ATOM    669  CA  ILE A  90      xxxxxx  33.110  31.221  1.00 25.76           C  "
    line_d = "ATOM    669  CA  ILE A  90      31.743   "
    line_e = "REMARK   1 REFERENCE 1                                                          "
    at_rec = AtomRecord(line_a, 10)
    show(devnull, at_rec)
    @test !at_rec.het_atom
    @test at_rec.serial == 669
    @test at_rec.atom_name == " CA "
    @test at_rec.alt_loc_id == ' '
    @test at_rec.res_name == "ILE"
    @test at_rec.chain_id == "A"
    @test at_rec.res_number == 90
    @test at_rec.ins_code == ' '
    @test at_rec.coords == [31.743, 33.110, 31.221]
    @test at_rec.occupancy == 1.00
    @test at_rec.temp_factor == 25.76
    @test at_rec.element == " C"
    @test at_rec.charge == "  "
    at_rec = AtomRecord(line_b)
    @test at_rec.het_atom
    @test at_rec.serial == 3474
    @test at_rec.atom_name == " O  "
    @test at_rec.alt_loc_id == 'B'
    @test at_rec.res_name == " XX"
    @test at_rec.chain_id == "A"
    @test at_rec.res_number == 334
    @test at_rec.ins_code == 'A'
    @test at_rec.coords == [8.802, 62.0, 8.672]
    @test at_rec.occupancy == 1.00
    @test at_rec.temp_factor == 39.15
    @test at_rec.element == " O"
    @test at_rec.charge == "1-"
    @test_throws PDBParseError AtomRecord(line_c)
    @test_throws PDBParseError AtomRecord(line_d)

    # Test parsing 1AKE (multiple chains, disordered atoms)
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat)
    @test structurename(struc) == "1AKE.pdb"
    @test countmodels(struc) == 1
    @test modelnumbers(struc) == [1]
    @test countchains(struc) == 2
    @test countchains(struc[1]) == 2
    @test chainids(struc) == ["A", "B"]
    @test chainids(struc[1]) == ["A", "B"]
    @test resname(struc['A'][10]) == "GLY"
    @test !isdisorderedres(struc['A'][10])
    @test serial(struc['A'][200]["NZ"]) == 1555
    @test !isdisorderedatom(struc['A'][200]["NZ"])
    @test isdisorderedatom(struc['A'][167]["CD"])
    @test altlocids(struc['A'][167]["CD"]) == ['A', 'B']
    @test x(struc['A'][167]["CD"]) == 24.502
    @test x(struc['A'][167]["CD"]['A']) == 24.502
    @test x(struc['A'][167]["CD"]['B']) == 24.69
    mos = collect(struc)
    @test modelnumber(mos[1]) == 1
    chs = collect(mos[1])
    @test chainid.(chs) == ["A", "B"]
    res = collect(chs[1])
    @test length(res) == 456
    @test resid(res[20]) == "20"
    ats = collect(res[20])
    @test length(ats) == 8
    @test atomname.(ats) == ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"]
    @test chainid(struc[1][begin]) == "A"
    @test chainid(struc[1][end  ]) == "B"
    @test resnumber(struc[1]["A"][begin]) == 1
    @test resnumber(struc[1]["A"][end  ]) == 543
    @test serial(struc[1]["A"][10][begin]) == 68
    @test serial(struc[1]["A"][10][end  ]) == 71
    @test serial(struc[1]["A"][167]["CD"][begin]) == 1288
    @test serial(struc[1]["A"][167]["CD"][end  ]) == 1289

    # Test choosedefaultaltlocid
    res = Residue("ALA", 1, ' ', false, Chain('A'))
    at_a = Atom(100, "CA", 'A', [1.0, 2.0, 3.0], 0.4, 10.0, "C", "", res)
    at_b = Atom(101, "CA", 'B', [1.0, 2.0, 3.0], 0.6, 10.0, "C", "", res)
    @test choosedefaultaltlocid(at_a, at_b) == 'B'
    @test choosedefaultaltlocid(at_b, at_a) == 'B'
    at_a = Atom(100, "CA", 'A', [1.0, 2.0, 3.0], 0.5, 10.0, "C", "", res)
    at_b = Atom(101, "CA", 'B', [1.0, 2.0, 3.0], 0.5, 10.0, "C", "", res)
    @test choosedefaultaltlocid(at_a, at_b) == 'A'
    @test choosedefaultaltlocid(at_b, at_a) == 'A'

    # Test applyselectors
    ats = collectatoms(struc)
    # Not providing any selector functions just returns the input list
    ats_min = applyselectors(ats)
    @test length(ats_min) == length(ats)
    @test serial.(ats_min) == serial.(ats)
    applyselectors!(ats_min)
    @test length(ats_min) == length(ats)
    @test serial.(ats_min) == serial.(ats)
    ats_min = applyselectors(ats, standardselector)
    @test length(ats_min) == 3312
    @test serial(ats_min[2000]) == 2006
    applyselectors!(ats, standardselector)
    @test length(ats) == 3312
    @test serial(ats[2000]) == 2006
    ats = collectatoms(struc)
    ats_min = applyselectors(ats, standardselector, disorderselector)
    @test length(ats_min) == 5
    @test serial(ats_min[4]) == 1294
    applyselectors!(ats, standardselector, disorderselector)
    @test length(ats) == 5
    @test serial(ats[4]) == 1294

    res = collectresidues(struc)
    res_min = applyselectors(res)
    @test length(res_min) == length(res)
    @test resid.(res_min; full=true) == resid.(res; full=true)
    applyselectors!(res_min)
    @test length(res_min) == length(res)
    @test resid.(res_min; full=true) == resid.(res; full=true)
    res_min = applyselectors(res, waterselector)
    @test length(res_min) == 378
    @test resid(res_min[300], full=true) == "H_657:B"
    applyselectors!(res, waterselector)
    @test length(res) == 378
    @test resid(res[300], full=true) == "H_657:B"
    res = collectresidues(struc)
    # Test anonymous selector function
    res_min = applyselectors(res, standardselector, res -> chainid(res) == "A")
    @test length(res_min) == 214
    @test resid(res_min[200], full=true) == "200:A"
    applyselectors!(res, standardselector, res -> chainid(res) == "A")
    @test length(res) == 214
    @test resid(res[200], full=true) == "200:A"

    # Test parsing options
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat, structure_name="New name")
    @test structurename(struc) == "New name"
    @test countatoms(struc) == 3804

    struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat, read_het_atoms=false)
    @test countatoms(struc) == 3312
    @test serial(collectatoms(struc)[2000]) == 2006
    @test sum(ishetero, collectatoms(struc)) == 0

    struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat, read_std_atoms=false)
    @test countatoms(struc) == 492
    @test serial(collectatoms(struc)[400]) == 3726
    @test sum(ishetero, collectatoms(struc)) == 492

    struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat, read_het_atoms=false, read_std_atoms=false)
    @test countatoms(struc) == 0
    @test countresidues(struc) == 0
    @test countchains(struc) == 0
    @test countmodels(struc) == 0

    struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat, remove_disorder=true)
    @test countatoms(struc) == 3804
    @test sum(isdisorderedatom, collectatoms(struc)) == 0
    @test tempfactor(struc['A'][167]["NE"]) == 23.32

    # Test parsing from stream
    open(testfilepath("PDB", "1AKE.pdb")) do file
        struc = read(file, PDBFormat)
        @test countatoms(struc) == 3804
        @test countresidues(struc) == 808
    end

    # Test parsing 1EN2 (disordered residue)
    struc = read(testfilepath("PDB", "1EN2.pdb"), PDBFormat)
    @test modelnumbers(struc) == [1]
    @test chainids(struc[1]) == ["A"]
    @test serial(struc['A'][48]["CA"]) == 394
    @test isdisorderedres(struc['A'][10])
    @test defaultresname(struc['A'][10]) == "SER"
    @test resname(disorderedres(struc['A'][10], "SER")) == "SER"
    @test resname(disorderedres(struc['A'][10], "GLY")) == "GLY"
    # Atoms in a disordered residue are not necessarily disordered
    @test !isdisorderedatom(struc['A'][10]["CA"])
    @test altlocid(struc['A'][10]["CA"]) == 'A'
    @test isdisorderedres(struc['A'][16])
    @test defaultresname(struc['A'][16]) == "ARG"
    @test resname(disorderedres(struc['A'][16], "ARG")) == "ARG"
    @test resname(disorderedres(struc['A'][16], "TRP")) == "TRP"
    @test !isdisorderedatom(disorderedres(struc['A'][16], "TRP")["CA"])
    @test altlocid(disorderedres(struc['A'][16], "TRP")["CA"]) == 'C'
    @test isdisorderedatom(struc['A'][16]["CA"])
    @test altlocids(struc['A'][16]["CA"]) == ['A', 'B']
    @test defaultaltlocid(struc['A'][16]["CA"]) == 'A'
    @test occupancy(struc['A'][16]["CA"]) == 0.22
    ats = collectatoms(DisorderedResidue[struc['A'][10], struc['A'][16]])
    @test length(ats) == 17
    @test isa(ats, Vector{AbstractAtom})
    @test serial(ats[10]) == 113
    @test isa(ats[10], DisorderedAtom)
    @test countatoms(DisorderedResidue[struc['A'][10], struc['A'][16]]) == 17
    res = collectresidues(DisorderedResidue[struc['A'][16], struc['A'][10]])
    @test length(res) == 2
    @test isa(res, Vector{DisorderedResidue})
    @test resnumber(res[1]) == 16
    @test countresidues(DisorderedResidue[struc['A'][16], struc['A'][10]]) == 2
    sort!(res)
    @test resnumber(res[1]) == 10
    @test serial(struc[1]["A"][10][begin]) == 57
    @test serial(struc[1]["A"][10][end  ]) == 62

    # Test parsing 1SSU (multiple models)
    struc = read(testfilepath("PDB", "1SSU.pdb"), PDBFormat)
    # Test countmodels
    @test countmodels(struc) == 20
    @test modelnumbers(struc) == collect(1:20)
    @test countchains(struc) == 1
    @test countchains(struc[5]) == 1
    @test chainids(struc) == ["A"]
    @test chainids(struc[5]) == ["A"]
    @test serial(struc[10]['A'][40]["HB2"]) == 574
    @test countatoms(struc) == 756
    @test countatoms.([struc[i] for i in modelnumbers(struc)]) == 756 * ones(Int, 20)
    @test countatoms(struc, hydrogenselector) == 357
    ats = collectatoms(Model[struc[5], struc[10]])
    @test length(ats) == 1512
    @test z(ats[20]) == -14.782
    @test z(ats[1000]) == -3.367
    @test countatoms(Model[struc[5], struc[10]]) == 1512
    res = collectresidues(Model[struc[5], struc[10]])
    @test length(res) == 102
    @test y(res[10]["O"]) == -1.612
    @test y(res[100]["O"]) == -13.184
    @test countresidues(Model[struc[5], struc[10]]) == 102
    chs = collectchains(Model[struc[5], struc[10]])
    @test length(chs) == 2
    @test chainid.(chs) == ["A", "A"]
    @test z(chs[2][5]["CA"]) == -5.667
    @test countchains(Model[struc[5], struc[10]]) == 2
    mos = collectmodels(Model[struc[10], struc[5]])
    @test length(mos) == 2
    @test modelnumber.(mos) == [10, 5]
    @test z(mos[2]['A'][5]["CA"]) == -5.837
    @test countmodels(Model[struc[10], struc[5]]) == 2
    @test modelnumber(struc[begin]) == 1
    @test modelnumber(struc[end  ]) == 20

    # Test collectatoms
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat)
    ats = collectatoms(struc)
    @test length(ats) == 3804
    @test isa(ats, Vector{AbstractAtom})
    @test isa(ats[70], Atom)
    @test isa(ats[1290], DisorderedAtom)
    @test serial(ats[1660]) == 3323
    ats = collectatoms(struc, heteroselector)
    @test length(ats) == 492
    @test serial(ats[80]) == 3463
    ats = collectatoms(struc, disorderselector)
    @test length(ats) == 12
    @test serial(ats[10]) == 3338
    @test all(isa.(ats, DisorderedAtom))
    ats = collectatoms(struc, standardselector, disorderselector)
    @test length(ats) == 5
    @test serial(ats[4]) == 1294
    ats = collectatoms(struc[1])
    @test length(ats) == 3804
    @test serial(ats[1660]) == 3323
    ats = collectatoms(struc['A'])
    @test length(ats) == 1954
    @test serial(ats[240]) == 240
    ats = collectatoms(struc['A'][50])
    @test length(ats) == 9
    @test serial(ats[4]) == 358
    ats = collectatoms(struc['A'][50]["CA"])
    @test length(ats) == 1
    @test isa(ats, Vector{AbstractAtom})
    @test serial(ats[1]) == 356
    ats = collectatoms(struc['A'][167]["CZ"])
    @test length(ats) == 1
    @test isa(ats, Vector{AbstractAtom})
    @test isa(ats[1], DisorderedAtom)
    @test serial(ats[1]) == 1292
    ats = collectatoms(Chain[struc['B'], struc['A']])
    @test length(ats) == 3804
    @test serial(ats[5]) == 1667
    ats = collectatoms(Residue[struc['A'][51], struc['A'][50]])
    @test length(ats) == 17
    @test serial(ats[5]) == 368
    ats = collectatoms(AbstractResidue[struc['A'][51], struc['A'][50]])
    @test length(ats) == 17
    @test serial(ats[5]) == 368
    ats = collectatoms(Atom[struc['A'][51]["CA"], struc['A'][50]["CA"]])
    @test length(ats) == 2
    @test isa(ats, Vector{Atom})
    @test serial(ats[2]) == 356
    ats = collectatoms(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]])
    @test length(ats) == 2
    @test isa(ats, Vector{DisorderedAtom})
    @test serial(ats[2]) == 1288
    ats = collectatoms(AbstractAtom[struc['A'][50]["CA"], struc['A'][167]["CZ"]])
    @test length(ats) == 2
    @test isa(ats, Vector{AbstractAtom})
    @test serial(ats[2]) == 1292
    ats = collectatoms(struc, expand_disordered=true)
    @test length(ats) == 3816
    @test isa(ats, Vector{AbstractAtom})
    @test all(isa.(ats, Atom))
    ats = collectatoms(struc, standardselector, expand_disordered=true)
    @test length(ats) == 3317
    ats = collectatoms(struc, disorderselector, expand_disordered=true)
    @test length(ats) == 24
    @test isa(ats, Vector{AbstractAtom})
    @test all(isa.(ats, Atom))

    # Test countatoms
    @test countatoms(struc) == 3804
    @test countatoms(struc[1]) == 3804
    @test countatoms(struc['A']) == 1954
    @test countatoms(struc['A'][50]) == 9
    @test countatoms(struc['A'][50]["CA"]) == 1
    @test countatoms([struc['A'], struc['B']]) == 3804
    @test countatoms(AbstractResidue[struc['A'][50], struc['A'][51]]) == 17
    @test countatoms(Residue[struc['A'][50], struc['A'][51]]) == 17
    @test countatoms(collectatoms(struc['A'])) == 1954
    @test countatoms(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]]) == 2
    @test countatoms(Atom[struc['A'][51]["CA"], struc['A'][50]["CA"]]) == 2

    @test countatoms(struc['A'], standardselector) == 1656
    @test countatoms(struc['A'], !standardselector) == 298
    @test countatoms(struc['A'], heteroselector) == 298
    @test countatoms(struc['A'], standardselector, disorderselector) == 5
    @test countatoms(struc, expand_disordered=true) == 3816
    @test countatoms(struc, standardselector, expand_disordered=true) == 3317

    @test countatoms(MolecularStructure()) == 0
    @test countatoms(Model()) == 0
    @test countatoms(Chain('X')) == 0
    @test countatoms(Residue("ALA", 100, ' ', false, Chain('A'))) == 0

    # Test collectresidues
    res = collectresidues(struc)
    @test length(res) == 808
    @test isa(res, Vector{AbstractResidue})
    @test isa(res[50], Residue)
    @test resnumber(res[220]) == 305
    res = collectresidues(struc, heteroselector)
    @test length(res) == 380
    @test resnumber(res[370]) == 725
    res = collectresidues(struc, standardselector, res -> chainid(res) == "A")
    @test length(res) == 214
    @test resnumber(res[200]) == 200
    res = collectresidues(struc[1])
    @test length(res) == 808
    @test resnumber(res[220]) == 305
    res = collectresidues(struc['A'])
    @test length(res) == 456
    @test resnumber(res[220]) == 305
    res = collectresidues(struc['A'][50])
    @test length(res) == 1
    @test isa(res, Vector{AbstractResidue})
    @test resnumber(res[1]) == 50
    res = collectresidues(struc['A'][50]["CA"])
    @test length(res) == 1
    @test isa(res, Vector{AbstractResidue})
    @test isa(res[1], Residue)
    @test resnumber(res[1]) == 50
    res = collectresidues(struc['A'][167]["CZ"])
    @test length(res) == 1
    @test isa(res, Vector{AbstractResidue})
    @test isa(res[1], Residue)
    @test resnumber(res[1]) == 167
    res = collectresidues(Chain[struc['B'], struc['A']])
    @test length(res) == 808
    @test resid(res[5], full=true) == "5:B"
    res = collectresidues(Residue[struc['A'][51], struc['A'][50]])
    @test length(res) == 2
    @test isa(res, Vector{Residue})
    @test resnumber(res[1]) == 51
    res = collectresidues(AbstractResidue[struc['A'][51], struc['A'][50]])
    @test length(res) == 2
    @test isa(res, Vector{AbstractResidue})
    @test resnumber(res[1]) == 51
    res = collectresidues(Atom[struc['A'][51]["CA"], struc['A'][50]["CA"]])
    @test length(res) == 2
    @test isa(res, Vector{AbstractResidue})
    @test resnumber(res[2]) == 50
    res = collectresidues(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]])
    @test length(res) == 1
    @test isa(res, Vector{AbstractResidue})
    @test atomnames(res[1]) == ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"]
    res = collectresidues(AbstractAtom[struc['A'][50]["CA"], struc['A'][167]["CZ"]])
    @test length(res) == 2
    @test isa(res, Vector{AbstractResidue})
    @test atomnames(res[1]) == ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"]
    struc = read(testfilepath("PDB", "1EN2.pdb"), PDBFormat)
    res = collectresidues(struc)
    @test length(res) == 166
    @test isa(res[10], DisorderedResidue)
    res = collectresidues(struc, expand_disordered=true)
    @test length(res) == 171
    @test isa(res, Vector{AbstractResidue})
    @test all(isa.(res, Residue))
    res = collectresidues(struc, standardselector, expand_disordered=true)
    @test length(res) == 90
    res = collectresidues(struc, disorderselector)
    @test length(res) == 5
    @test isa(res, Vector{AbstractResidue})
    @test all(isa.(res, DisorderedResidue))
    res = collectresidues(struc, disorderselector, expand_disordered=true)
    @test length(res) == 10
    @test isa(res, Vector{AbstractResidue})
    @test all(isa.(res, Residue))

    # Test countresidues
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat)
    @test countresidues(struc) == 808
    @test countresidues(struc[1]) == 808
    @test countresidues(struc['A']) == 456
    @test countresidues(struc['A'][50]) == 1
    @test countresidues(struc['A'][50]["CA"]) == 1
    @test countresidues([struc['A'], struc['B']]) == 808
    @test countresidues(AbstractResidue[struc['A'][50], struc['A'][51]]) == 2
    @test countresidues(Residue[struc['A'][50], struc['A'][51]]) == 2
    @test countresidues(collectatoms(struc['A'])) == 456
    @test countresidues(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]]) == 1
    @test countresidues(Atom[struc['A'][51]["CA"], struc['A'][50]["CA"]]) == 2

    @test countresidues(struc['A'], standardselector) == 214
    @test countresidues(struc['A'], heteroselector) == 242
    @test countresidues(struc, standardselector, res -> chainid(res) == "A") == 214
    struc = read(testfilepath("PDB", "1EN2.pdb"), PDBFormat)
    @test countresidues(struc, expand_disordered=true) == 171
    @test countresidues(struc, standardselector, expand_disordered=true) == 90

    @test countresidues(MolecularStructure()) == 0
    @test countresidues(Model()) == 0
    @test countresidues(Chain('X')) == 0
    @test countresidues(Residue("ALA", 100, ' ', false, Chain('A'))) == 1

    # Test collectchains
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat)
    chs = collectchains(struc)
    @test length(chs) == 2
    @test isa(chs, Vector{Chain})
    @test chainid(chs[2]) == "B"
    chs = collectchains(struc, ch -> chainid(ch) == "B")
    @test length(chs) == 1
    @test chainid(chs[1]) == "B"
    chs = collectchains(struc[1])
    @test length(chs) == 2
    @test chainid(chs[2]) == "B"
    chs = collectchains(struc['A'])
    @test length(chs) == 1
    @test chainid(chs[1]) == "A"
    chs = collectchains(struc['A'][50])
    @test length(chs) == 1
    @test chainid(chs[1]) == "A"
    chs = collectchains(struc['A'][50]["CA"])
    @test length(chs) == 1
    @test chainid(chs[1]) == "A"
    chs = collectchains(Chain[struc['B'], struc['A']])
    @test length(chs) == 2
    @test chainid(chs[2]) == "A"
    chs = collectchains(Residue[struc['A'][51], struc['B'][50]])
    @test length(chs) == 2
    @test chainid(chs[2]) == "B"
    chs = collectchains(AbstractResidue[struc['A'][51], struc['B'][50]])
    @test length(chs) == 2
    @test chainid(chs[2]) == "B"
    chs = collectchains(Atom[struc['B'][51]["CA"], struc['A'][50]["CA"]])
    @test length(chs) == 2
    @test chainid(chs[2]) == "A"
    chs = collectchains(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]])
    @test length(chs) == 1
    @test chainid(chs[1]) == "A"
    chs = collectchains(AbstractAtom[struc['A'][50]["CA"], struc['A'][167]["CZ"]])
    @test length(chs) == 1
    @test chainid(chs[1]) == "A"

    # Test countchains
    @test countchains(struc) == 2
    @test countchains(struc[1]) == 2
    @test countchains(struc['A']) == 1
    @test countchains(struc['A'][50]) == 1
    @test countchains(struc['A'][50]["CA"]) == 1
    @test countchains([struc['A'], struc['B']]) == 2
    @test countchains(AbstractResidue[struc['A'][50], struc['B'][51]]) == 2
    @test countchains(Residue[struc['A'][50], struc['B'][51]]) == 2
    @test countchains(collectatoms(struc['A'])) == 1
    @test countchains(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]]) == 1
    @test countchains(Atom[struc['A'][51]["CA"], struc['B'][50]["CA"]]) == 2

    @test countchains(struc, ch -> chainid(ch) == "B") == 1

    @test countchains(MolecularStructure()) == 0
    @test countchains(Model()) == 0
    @test countchains(Chain('X')) == 1
    @test countchains(Residue("ALA", 100, ' ', false, Chain('A'))) == 1

    # Test collectmodels
    struc_1SSU = read(testfilepath("PDB", "1SSU.pdb"), PDBFormat)
    mos = collectmodels(struc_1SSU)
    @test length(mos) == 20
    @test isa(mos, Vector{Model})
    @test modelnumber.(mos) == collect(1:20)
    mos = collectmodels(struc_1SSU, mo -> modelnumber(mo) < 4)
    @test length(mos) == 3
    @test modelnumber(mos[2]) == 2
    mos = collectmodels(struc_1SSU[10])
    @test length(mos) == 1
    @test modelnumber(mos[1]) == 10
    mos = collectmodels(struc_1SSU['A'])
    @test length(mos) == 1
    @test modelnumber(mos[1]) == 1
    mos = collectmodels(struc_1SSU[7]['A'][50])
    @test length(mos) == 1
    @test modelnumber(mos[1]) == 7
    mos = collectmodels(struc_1SSU['A'][50]["CA"])
    @test length(mos) == 1
    @test modelnumber(mos[1]) == 1
    mos = collectmodels(Chain[struc['B'], struc['A']])
    @test length(mos) == 1
    @test modelnumber(mos[1]) == 1
    mos = collectmodels(Residue[struc['A'][51], struc['B'][50]])
    @test length(mos) == 1
    @test modelnumber(mos[1]) == 1
    mos = collectmodels(AbstractResidue[struc_1SSU[12]['A'][51], struc_1SSU['A'][50]])
    @test length(mos) == 2
    @test modelnumber(mos[2]) == 1
    mos = collectmodels(Atom[struc_1SSU['A'][51]["CA"], struc_1SSU['A'][50]["CA"]])
    @test length(mos) == 1
    @test modelnumber(mos[1]) == 1
    mos = collectmodels(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]])
    @test length(mos) == 1
    @test modelnumber(mos[1]) == 1
    mos = collectmodels(AbstractAtom[struc['A'][50]["CA"], struc['A'][167]["CZ"]])
    @test length(mos) == 1
    @test modelnumber(mos[1]) == 1

    # Test countmodels
    @test countmodels(struc_1SSU) == 20
    @test countmodels(struc_1SSU[10]) == 1
    @test countmodels(struc_1SSU['A']) == 1
    @test countmodels(struc_1SSU['A'][50]) == 1
    @test countmodels(struc_1SSU[10]['A'][50]["CA"]) == 1
    @test countmodels([struc['A'], struc['B']]) == 1
    @test countmodels(AbstractResidue[struc_1SSU[1]['A'][50], struc_1SSU[2]['A'][51]]) == 2
    @test countmodels(Residue[struc['A'][50], struc['B'][51]]) == 1
    @test countmodels(collectatoms(struc_1SSU)) == 1
    @test countmodels(collectatoms(collect(struc_1SSU))) == 20
    @test countmodels(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]]) == 1
    @test countmodels(Atom[struc_1SSU[7]['A'][51]["CA"], struc_1SSU['A'][50]["CA"]]) == 2

    @test countmodels(struc_1SSU, mo -> modelnumber(mo) < 4) == 3

    @test countmodels(MolecularStructure()) == 0
    @test countmodels(Model()) == 1
    @test countmodels(Chain('X')) == 1
    @test countmodels(Residue("ALA", 100, ' ', false, Chain('A'))) == 1

    # Test parser error handling
    error = PDBParseError("message", 10, "line")
    showerror(devnull, error)
    # Missing coordinate (blank string)
    @test_throws PDBParseError read(testfilepath("PDB", "1AKE_err_a.pdb"), PDBFormat)
    # Missing chain ID (line ends early)
    @test_throws PDBParseError read(testfilepath("PDB", "1AKE_err_b.pdb"), PDBFormat)
    # Bad MODEL record
    @test_throws PDBParseError read(testfilepath("PDB", "1SSU_err.pdb"), PDBFormat)
    # Truncated MODEL record from ASTRAL/SCOP95
    struc = read(testfilepath("PDB", "d9pcya_.ent"), PDBFormat)
    @test isa(struc, MolecularStructure)
    @test countmodels(struc) == 16
    # Duplicate atom names in same residue
    @test_throws ErrorException read(testfilepath("PDB", "1AKE_err_c.pdb"), PDBFormat)
    # Non-existent file
    @test_throws SystemError read(testfilepath("PDB", "non_existent_file.pdb"), PDBFormat)

    # Test parsing empty file
    struc = read(IOBuffer(""), PDBFormat)
    @test isa(struc, MolecularStructure)
    @test countmodels(struc) == 0
end

@testset "PDB writing" begin
    # Test spacestring
    @test spacestring(1.5, 5) == "  1.5"
    @test spacestring("A", 3) == "  A"
    @test spacestring('A', 3) == "  A"
    @test_throws ArgumentError spacestring(1.456789, 5)
    @test_throws ArgumentError spacestring("ABCDEF", 3)

    # Test spaceatomname
    res = Residue("ALA", 1, ' ', false, Chain('A'))
    @test spaceatomname(Atom(1, " CA ", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " C", "  ", res)) == " CA "
    @test spaceatomname(Atom(1, "N",    ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " N", "  ", res)) == " N  "
    @test spaceatomname(Atom(1, "N",    ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "  ", "  ", res)) == " N  "
    @test spaceatomname(Atom(1, "CA",   ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " C", "  ", res)) == " CA "
    @test spaceatomname(Atom(1, "NE1",  ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " N", "  ", res)) == " NE1"
    @test spaceatomname(Atom(1, "2HD1", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " H", "  ", res)) == "2HD1"
    @test spaceatomname(Atom(1, "HH11", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " H", "  ", res)) == "HH11"
    @test spaceatomname(Atom(1, "1H",   ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " H", "  ", res)) == "1H  "
    @test spaceatomname(Atom(1, "MG",   ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "MG", "  ", res)) == "MG  "
    @test spaceatomname(Atom(1, "MG",   ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "  ", "  ", res)) == " MG "
    @test_throws ArgumentError spaceatomname(Atom(1, "11H",   ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " H", "  ", res))
    @test_throws ArgumentError spaceatomname(Atom(1, "11H11", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " H", "  ", res))
    @test_throws ArgumentError spaceatomname(Atom(1, "1MG",   ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "MG", "  ", res))

    # Test output formatting
    @test pyfmt(coordspec, 5.0) == "5.000"
    @test pyfmt(floatspec, -10.0) == "-10.00"

    # Test pdbline
    ch_a = Chain('A')
    ch_a["1"] = Residue("ALA", 1, ' ', false, ch_a)
    ch_a["1"][" N  "] = Atom(10, " N  ", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " N", "  ", ch_a["1"])
    line = pdbline(ch_a["1"][" N  "])
    @test line == "ATOM     10  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  "
    ch_b = Chain('B')
    ch_b["H_20"] = Residue("X", 20, ' ', true, ch_b)
    ch_b["H_20"]["C"] = Atom(101, "C", 'A', [10.5, 20.12345, -5.1227], 0.50, 50.126, "C", "1+", ch_b["H_20"])
    line = pdbline(ch_b["H_20"]["C"])
    @test line == "HETATM  101  C  A  X B  20      10.500  20.123  -5.123  0.50 50.13           C1+"
    ch_b["H_20"]["11H11"] = Atom(1, "11H11", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " H", "  ", ch_b["H_20"])
    @test_throws ArgumentError pdbline(ch_b["H_20"]["11H11"])
    ch_b["H_20"]["H1"] = Atom(1, "H1", ' ', [-1000.123, 0.0, 0.0], 1.0, 0.0, " H", "  ", ch_b["H_20"])
    @test_throws ArgumentError pdbline(ch_b["H_20"]["H1"])

    struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat)
    @test pdbline(struc["A"][167]["NH1"]) == "ATOM   1294  NH1AARG A 167      24.181  40.144  13.699  0.50 27.31           N  "

    line_a = "ATOM    669  CA  ILE A  90      31.743  33.110  31.221  1.00 25.76           C  "
    line_b = "HETATM 3474  O  B XX A 334A      8.802  62.000   8.672  1.00 39.15           O1-"
    at_rec = AtomRecord(line_a)
    @test pdbline(at_rec) == "ATOM    669  CA  ILE A  90      31.743  33.110  31.221  1.00 25.76           C  "
    at_rec = AtomRecord(line_b)
    @test pdbline(at_rec) == "HETATM 3474  O  B XX A 334A      8.802  62.000   8.672  1.00 39.15           O1-"

    # Test writepdb
    struc = read(testfilepath("PDB", "1SSU.pdb"), PDBFormat)
    writepdb(temp_filename, struc)
    @test countlines(temp_filename) == 15160
    struc_written = read(temp_filename, PDBFormat)
    @test isa(struc_written, MolecularStructure)
    @test modelnumbers(struc_written) == collect(1:20)
    @test countatoms(struc_written) == 756
    @test z(struc_written[4]['A']["30"]["OG"]) == -2.177
    @test atomnames(struc_written[15]['A']["39"]) == [
        "N", "CA", "C", "O", "CB", "SG", "H", "HA", "HB2", "HB3"]

    # Test writing to stream
    open(temp_filename, "w") do file
        writepdb(file, struc)
    end
    @test countlines(temp_filename) == 15160
    struc_written = read(temp_filename, PDBFormat)
    @test modelnumbers(struc_written) == collect(1:20)
    @test countatoms(struc_written) == 756
    @test z(struc_written[4]['A']["30"]["OG"]) == -2.177
    @test atomnames(struc_written[15]['A']["39"]) == [
        "N", "CA", "C", "O", "CB", "SG", "H", "HA", "HB2", "HB3"]

    # Test selectors
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat)
    writepdb(temp_filename, struc, heteroselector)
    @test countlines(temp_filename) == 499
    struc_written = read(temp_filename, PDBFormat)
    @test modelnumbers(struc_written) == [1]
    @test countatoms(struc_written) == 492
    @test chainids(struc_written) == ["A", "B"]
    @test tempfactor(struc_written['B']["H_705"]["O"]) == 64.17
    writepdb(temp_filename, collectatoms(struc, standardselector,
                                            disorderselector))
    @test countlines(temp_filename) == 10
    struc_written = read(temp_filename, PDBFormat)
    @test countatoms(struc_written) == 5
    @test sum(isdisorderedatom, collectatoms(struc_written)) == 5
    @test defaultaltlocid(struc_written['A'][167]["NH1"]) == 'A'
    writepdb(temp_filename, struc)
    @test countlines(temp_filename) == 3816
    writepdb(temp_filename, struc, expand_disordered=false)
    @test countlines(temp_filename) == 3804
    struc_written = read(temp_filename, PDBFormat)
    @test !any(isdisorderedatom.(collectatoms(struc_written)))

    # Test writing different element types
    writepdb(temp_filename, struc[1])
    @test countlines(temp_filename) == 3816
    struc_written = read(temp_filename, PDBFormat)
    @test modelnumbers(struc_written) == [1]
    @test countatoms(struc_written) == 3804
    writepdb(temp_filename, struc['A'])
    @test countlines(temp_filename) == 1966
    struc_written = read(temp_filename, PDBFormat)
    @test chainids(struc_written) == ["A"]
    writepdb(temp_filename, struc['A'][50])
    @test countlines(temp_filename) == 9
    struc_written = read(temp_filename, PDBFormat)
    @test chainids(struc_written) == ["A"]
    @test countresidues(struc_written) == 1
    @test atomnames(struc_written['A'][50]) == [
        "N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"]
    writepdb(temp_filename, struc['A'][50]["CA"])
    @test countlines(temp_filename) == 1
    struc_written = read(temp_filename, PDBFormat)
    @test countatoms(struc_written) == 1
    @test !isdisorderedatom(collectatoms(struc_written)[1])
    @test serial(collectatoms(struc_written)[1]) == 356
    writepdb(temp_filename, struc['A'][167]["CZ"])
    @test countlines(temp_filename) == 2
    struc_written = read(temp_filename, PDBFormat)
    @test countatoms(struc_written) == 1
    @test isdisorderedatom(collectatoms(struc_written)[1])
    @test length(collectatoms(struc_written)[1]) == 2
    @test tempfactor(collectatoms(struc_written)[1]) == 16.77
    writepdb(temp_filename, Chain[struc['A'], struc['B']])
    @test countlines(temp_filename) == 3816
    struc_written = read(temp_filename, PDBFormat)
    @test chainids(struc_written) == ["A", "B"]
    @test countatoms(struc_written['A']) == 1954
    @test countatoms(struc_written['B']) == 1850
    @test altlocids(struc_written['A']["H_215"]["O1G"]) == ['A', 'B']
    writepdb(temp_filename, AbstractResidue[struc['A'][51], struc['A'][50]])
    @test countlines(temp_filename) == 17
    struc_written = read(temp_filename, PDBFormat)
    @test countresidues(struc_written) == 2
    @test resid.(collectresidues(struc_written)) == ["50", "51"]
    @test countatoms(struc_written) == 17
    writepdb(temp_filename, AbstractAtom[struc['A'][51]["CA"], struc['A'][50]["CA"]])
    @test countlines(temp_filename) == 2
    struc_written = read(temp_filename, PDBFormat)
    @test countatoms(struc_written) == 2
    @test !ishetero(struc_written['A'][51]["CA"])
    writepdb(temp_filename, struc['A'][51]["N"], calphaselector)
    @test countlines(temp_filename) == 0
    struc_written = read(temp_filename, PDBFormat)
    @test countatoms(struc_written) == 0
    writepdb(temp_filename, MolecularStructure())
    @test countlines(temp_filename) == 0
    struc_written = read(temp_filename, PDBFormat)
    @test countatoms(struc_written) == 0

    # Test multiple model writing
    struc = read(testfilepath("PDB", "1SSU.pdb"), PDBFormat)
    writepdb(temp_filename, Model[struc[10], struc[5]])
    @test countlines(temp_filename) == 1516
    struc_written = read(temp_filename, PDBFormat)
    @test modelnumbers(struc_written) == [5, 10]
    @test modelnumber(defaultmodel(struc_written)) == 5
    @test countatoms(struc_written[5]) == 756
    @test countatoms(struc_written[10]) == 756
    @test_throws KeyError struc_written[1]

    # Test disordered residue writing
    struc = read(testfilepath("PDB", "1EN2.pdb"), PDBFormat)
    writepdb(temp_filename, struc)
    @test countlines(temp_filename) == 819
    struc_written = read(temp_filename, PDBFormat)
    @test countatoms(struc_written) == 754
    @test isa(struc_written['A'][15], Residue)
    @test isa(struc_written['A'][16], DisorderedResidue)
    @test defaultresname(struc_written['A'][16]) == "ARG"
    @test isa(struc_written['A'][16]["N"], DisorderedAtom)
    @test defaultaltlocid(struc_written['A'][16]["N"]) == 'A'
    @test isa(disorderedres(struc_written['A'][16], "TRP")["N"], Atom)
    @test countatoms(struc_written['A'][16]) == 11
    @test countatoms(disorderedres(struc_written['A'][16], "TRP")) == 14
    writepdb(temp_filename, AbstractResidue[struc['A'][16], struc['A'][10]])
    @test countlines(temp_filename) == 46
    struc_written = read(temp_filename, PDBFormat)
    @test countresidues(struc_written) == 2
    @test isa(struc_written['A'][10], DisorderedResidue)
    @test isa(struc_written['A'][16], DisorderedResidue)
    @test defaultresname(struc_written['A'][10]) == "SER"
    @test isa(disorderedres(struc_written['A'][10], "GLY")["O"], Atom)
    @test altlocid(disorderedres(struc_written['A'][10], "GLY")["O"]) == 'B'
    @test countatoms(struc_written['A'][10]) == 6
    @test countatoms(struc_written['A'][16]) == 11
    writepdb(temp_filename, struc, expand_disordered=false)
    @test countlines(temp_filename) == 754
    struc_written = read(temp_filename, PDBFormat)
    @test !any(isdisorderedres.(collectresidues(struc_written)))
    @test !any(isdisorderedatom.(collectatoms(struc_written)))
    @test countresidues(struc_written, expand_disordered=true) == 166
    @test countatoms(struc_written, expand_disordered=true) == 754

    @test_throws ArgumentError writepdb(temp_filename, Atom(
        1, "11H11", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " H", "  ", res))

    checkchainerror(Chain("A"))
    @test_throws ArgumentError checkchainerror(Chain("AA"))
    res = Residue("ALA", 10, ' ', false, [], Dict(), Chain("AA"), '-')
    res["CA"] = Atom(100, " CA ", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res)
    push!(res.atom_list, "CA")
    @test_throws ArgumentError writepdb(temp_filename, res)
end

@testset "mmCIF" begin
    # Test mmCIF dictionary
    dic = MMCIFDict()
    dic = MMCIFDict(Dict())
    show(IOBuffer(), MIME("text/plain"), MMCIFDict())
    show(IOBuffer(), MIME("text/plain"), MMCIFDict(Dict("a" => ["b"])))

    mmcif_1ake = testfilepath("mmCIF", "1AKE.cif")
    gzip_file(mmcif_1ake, temp_filename)
    for dic in (MMCIFDict(mmcif_1ake), MMCIFDict(temp_filename; gzip=true))
        @test isa(dic.dict, Dict{String, Vector{String}})
        @test dic["_pdbx_database_status.recvd_initial_deposition_date"] == ["1991-11-08"]
        @test dic["_audit_author.name"] == ["Mueller, C.W.", "Schulz, G.E."]
        @test length(dic["_atom_site.group_PDB"]) == 3816
        dic["_pdbx_database_status.recvd_initial_deposition_date"] = ["changed"]
        @test dic["_pdbx_database_status.recvd_initial_deposition_date"] == ["changed"]
        @test length(keys(dic)) == 610
        @test length(values(dic)) == 610
        @test haskey(dic, "_cell.entry_id")
        @test !haskey(dic, "nokey")
        @test get(dic, "_cell.entry_id", ["default"]) == ["1AKE"]
        @test get(dic, "nokey", ["default"]) == ["default"]
        @test ismissing(get(dic, "nokey", missing))
        show(devnull, dic)
    end

    multiline_str = """
        data_test
        _test_value
        ;first line
            second line
        third line
        ;
        """
    dic = MMCIFDict(IOBuffer(multiline_str))
    @test dic["data_"] == ["test"]
    @test dic["_test_value"] == ["first line\n    second line\nthird line"]

    gz = GzipCompressorStream(IOBuffer(multiline_str))
    dic = MMCIFDict(gz; gzip=true)
    @test dic["data_"] == ["test"]
    @test dic["_test_value"] == ["first line\n    second line\nthird line"]
    close(gz)

    comment_str = """
        data_test
        _test_single foo # Ignore this comment
        loop_
        _test_loop_1
        _test_loop_2
        a b # Ignore this comment
        c d
        """
    dic = MMCIFDict(IOBuffer(comment_str))
    @test dic["_test_single"] == ["foo"]
    @test dic["_test_loop_1"] == ["a", "c"]
    @test dic["_test_loop_2"] == ["b", "d"]

    quote_str = """
        data_1MOM
        loop_
        _struct_conf.conf_type_id
        _struct_conf.id
        _struct_conf.pdbx_PDB_helix_id
        _struct_conf.beg_label_comp_id
        _struct_conf.beg_label_asym_id
        _struct_conf.beg_label_seq_id
        _struct_conf.pdbx_beg_PDB_ins_code
        _struct_conf.end_label_comp_id
        _struct_conf.end_label_asym_id
        _struct_conf.end_label_seq_id
        _struct_conf.pdbx_end_PDB_ins_code
        _struct_conf.beg_auth_comp_id
        _struct_conf.beg_auth_asym_id
        _struct_conf.beg_auth_seq_id
        _struct_conf.end_auth_comp_id
        _struct_conf.end_auth_asym_id
        _struct_conf.end_auth_seq_id
        _struct_conf.pdbx_PDB_helix_class
        _struct_conf.details
        _struct_conf.pdbx_PDB_helix_length
        HELX_P HELX_P1  A     PRO A 11  ? ASN A 23  ? PRO A 11  ASN A 23  1 ?                         13
        HELX_P HELX_P2  "A'"  ALA A 44  ? ARG A 46  ? ALA A 44  ARG A 46  5 ?                         3
        HELX_P HELX_P3  B     PRO A 86  ? SER A 92  ? PRO A 86  SER A 92  1 ?                         7
        HELX_P HELX_P4  C     TYR A 111 ? ALA A 118 ? TYR A 111 ALA A 118 1 ?                         8
        HELX_P HELX_P5  "B'"  ARG A 122 ? LYS A 124 ? ARG A 122 LYS A 124 5 ?                         3
        HELX_P HELX_P6  D     LEU A 129 ? LEU A 139 ? LEU A 129 LEU A 139 1 ?                         11
        HELX_P HELX_P7  E     SER A 144 ? ILE A 155 ? SER A 144 ILE A 155 1 ?                         12
        HELX_P HELX_P8  "C'"  THR A 158 ? ARG A 163 ? THR A 158 ARG A 163 1 ?                         6
        HELX_P HELX_P9  F     LYS A 165 ? GLU A 173 ? LYS A 165 GLU A 173 1 ?                         9
        HELX_P HELX_P10 G     LEU A 183 ? SER A 191 ? LEU A 183 SER A 191 1 ?                         9
        HELX_P HELX_P11 H     TRP A 192 ? LEU A 201 ? TRP A 192 LEU A 201 1 ?                         10
        HELX_P HELX_P12 "D'"  LYS A 231 ? ASN A 236 ? LYS A 231 ASN A 236 3 ?                         6
        HELX_P HELX_P13 "E'"  THR A 243 ? ILE A 246 ? THR A 243 ILE A 246 5 ?                         4
        TURN_P TURN_P1  'A'"' ASP A 1   ? ILE A 34  ? ASP A 1   ILE A 34  5
        ;TYPE I'
        ;
        ?
        TURN_P TURN_P2  BC    ASP A 1   ? GLY A 57  ? ASP A 1   GLY A 57  5 'TYPE I'                  ?
        TURN_P TURN_P3  CD    ASP A 1   ? ASN A 68  ? ASP A 1   ASN A 68  5 'TYPE I (H-BOND O65-N69)' ?
        TURN_P TURN_P4  DE    ASP A 1   ? THR A 79  ? ASP A 1   THR A 79  5
        ;TYPE II'
        ;
        ?
        """
    dic = MMCIFDict(IOBuffer(quote_str))
    @test dic["_struct_conf.pdbx_PDB_helix_id"] == [
        "A", "A'", "B", "C", "B'", "D", "E", "C'",
        "F", "G", "H", "D'", "E'", "A'\"", "BC", "CD", "DE"
    ]
    @test dic["_struct_conf.details"] == [
        "?", "?", "?", "?", "?", "?", "?", "?", "?", "?", "?", "?", "?",
        "TYPE I'", "TYPE I", "TYPE I (H-BOND O65-N69)", "TYPE II'"
    ]

    underscore_str = """
        data_4Q9R
        loop_
        _pdbx_audit_revision_item.ordinal
        _pdbx_audit_revision_item.revision_ordinal
        _pdbx_audit_revision_item.data_content_type
        _pdbx_audit_revision_item.item
        1  5 'Structure model' '_atom_site.B_iso_or_equiv'
        2  5 'Structure model' '_atom_site.Cartn_x'
        3  5 'Structure model' '_atom_site.Cartn_y'
        4 5 'Structure model' '_atom_site.Cartn_z'
        """
    dic = MMCIFDict(IOBuffer(underscore_str))
    @test length(keys(dic)) == 5
    @test dic["_pdbx_audit_revision_item.item"] == [
        "_atom_site.B_iso_or_equiv", "_atom_site.Cartn_x",
        "_atom_site.Cartn_y", "_atom_site.Cartn_z"
    ]

    keyerr_str = """
        data_test
        loop_
        bad.key
        bad.key2
        1 2
        3 4
        """
    @test_throws ArgumentError MMCIFDict(IOBuffer(keyerr_str))

    # Test splitline
    @test splitline("foo bar") == ["foo", "bar"]
    @test splitline("  foo bar  ") == ["foo", "bar"]
    @test splitline("'foo' bar") == ["foo", "bar"]
    @test splitline("foo \"bar\"") == ["foo", "bar"]
    @test splitline("foo 'bar a' b") == ["foo", "bar a", "b"]
    @test splitline("foo 'bar'a' b") == ["foo", "bar'a", "b"]
    @test splitline("foo \"bar' a\" b") == ["foo", "bar' a", "b"]
    @test splitline("foo '' b") == ["foo", "", "b"]
    @test splitline("foo #bar") == ["foo"]
    @test splitline("foo b#ar") == ["foo", "b#ar"]
    @test_throws ArgumentError splitline("foo 'bar")
    @test_throws ArgumentError splitline("foo 'ba'r  ")
    @test_throws ArgumentError splitline("foo \"bar'")
    @test_throws ArgumentError splitline("foo b'ar'")

    # Test tokenizecif and tokenizecifstructure
    open(testfilepath("mmCIF", "1AKE.cif")) do f
        tokens = tokenizecif(f)
        @test length(tokens) == 93983
        @test tokens[90000] == "HOH"
    end

    open(testfilepath("mmCIF", "1AKE.cif")) do f
        tokens = tokenizecifstructure(f)
        @test length(tokens) == 80158
        @test tokens[1] == "loop_"
        @test tokens[10] == "_atom_site.label_seq_id"
        @test tokens[80089] == "77.69"
    end

    # Test parsing empty file
    dic = MMCIFDict(IOBuffer(""))
    @test isa(dic, MMCIFDict)
    @test length(keys(dic)) == 0
    @test length(values(dic)) == 0

    # Test AtomRecord
    at_rec = AtomRecord(MMCIFDict(testfilepath("mmCIF", "1AKE.cif")), 5)
    show(devnull, at_rec)
    @test !at_rec.het_atom
    @test at_rec.serial == 5
    @test at_rec.atom_name == "CB"
    @test at_rec.alt_loc_id == ' '
    @test at_rec.res_name == "MET"
    @test at_rec.chain_id == "A"
    @test at_rec.res_number == 1
    @test at_rec.ins_code == ' '
    @test at_rec.coords == [24.677, 53.310, 39.580]
    @test at_rec.occupancy == 1.00
    @test at_rec.temp_factor == 38.06
    @test at_rec.element == "C"
    @test at_rec.charge == "  "

    # Test parsing 1AKE
    struc = read(testfilepath("mmCIF", "1AKE.cif"), MMCIFFormat)
    @test structurename(struc) == "1AKE.cif"
    @test countmodels(struc) == 1
    @test modelnumbers(struc) == [1]
    @test countchains(struc) == 2
    @test countchains(struc[1]) == 2
    @test chainids(struc) == ["A", "B"]
    @test chainids(struc[1]) == ["A", "B"]
    @test resname(struc['A'][10]) == "GLY"
    @test !isdisorderedres(struc['A'][10])
    @test serial(struc['A'][200]["NZ"]) == 1555
    @test !isdisorderedatom(struc['A'][200]["NZ"])
    @test isdisorderedatom(struc['A'][167]["CD"])
    @test altlocids(struc['A'][167]["CD"]) == ['A', 'B']
    @test x(struc['A'][167]["CD"]) == 24.502
    @test x(struc['A'][167]["CD"]['A']) == 24.502
    @test x(struc['A'][167]["CD"]['B']) == 24.69
    mos = collect(struc)
    @test modelnumber(mos[1]) == 1
    chs = collect(mos[1])
    @test chainid.(chs) == ["A", "B"]
    res = collect(chs[1])
    @test length(res) == 456
    @test resid(res[20]) == "20"
    ats = collect(res[20])
    @test length(ats) == 8
    @test atomname.(ats) == ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"]

    # Test parsing options
    struc = read(testfilepath("mmCIF", "1AKE.cif"), MMCIFFormat, structure_name="New name")
    @test structurename(struc) == "New name"
    @test countatoms(struc) == 3804

    struc = read(testfilepath("mmCIF", "1AKE.cif"), MMCIFFormat, read_het_atoms=false)
    @test countatoms(struc) == 3312
    # Different to the PDB file due to the lack of TER label serial
    @test serial(collectatoms(struc)[2000]) == 2005
    @test sum(ishetero, collectatoms(struc)) == 0

    struc = read(testfilepath("mmCIF", "1AKE.cif"), MMCIFFormat, read_std_atoms=false)
    @test countatoms(struc) == 492
    @test serial(collectatoms(struc)[400]) == 3724
    @test sum(ishetero, collectatoms(struc)) == 492

    struc = read(testfilepath("mmCIF", "1AKE.cif"), MMCIFFormat, read_het_atoms=false,
                read_std_atoms=false)
    @test countatoms(struc) == 0
    @test countresidues(struc) == 0
    @test countchains(struc) == 0
    @test countmodels(struc) == 0

    struc = read(testfilepath("mmCIF", "1AKE.cif"), MMCIFFormat, remove_disorder=true)
    @test countatoms(struc) == 3804
    @test sum(isdisorderedatom, collectatoms(struc)) == 0
    @test tempfactor(struc['A'][167]["NE"]) == 23.32

    # Test parsing from stream
    open(testfilepath("mmCIF", "1AKE.cif")) do file
        struc = read(file, MMCIFFormat)
        @test countatoms(struc) == 3804
        @test countresidues(struc) == 808
    end

    # Test parsing 1EN2
    struc = read(testfilepath("mmCIF", "1EN2.cif"), MMCIFFormat)
    @test modelnumbers(struc) == [1]
    @test chainids(struc[1]) == ["A"]
    @test serial(struc['A'][48]["CA"]) == 394
    @test isdisorderedres(struc['A'][10])
    @test defaultresname(struc['A'][10]) == "SER"
    @test resname(disorderedres(struc['A'][10], "SER")) == "SER"
    @test resname(disorderedres(struc['A'][10], "GLY")) == "GLY"
    @test !isdisorderedatom(struc['A'][10]["CA"])
    @test altlocid(struc['A'][10]["CA"]) == 'A'
    @test isdisorderedres(struc['A'][16])
    @test defaultresname(struc['A'][16]) == "ARG"
    @test resname(disorderedres(struc['A'][16], "ARG")) == "ARG"
    @test resname(disorderedres(struc['A'][16], "TRP")) == "TRP"
    @test !isdisorderedatom(disorderedres(struc['A'][16], "TRP")["CA"])
    @test altlocid(disorderedres(struc['A'][16], "TRP")["CA"]) == 'C'
    @test isdisorderedatom(struc['A'][16]["CA"])
    @test altlocids(struc['A'][16]["CA"]) == ['A', 'B']
    @test defaultaltlocid(struc['A'][16]["CA"]) == 'A'
    @test occupancy(struc['A'][16]["CA"]) == 0.22
    ats = collectatoms(DisorderedResidue[struc['A'][10], struc['A'][16]])
    @test length(ats) == 17
    @test isa(ats, Vector{AbstractAtom})
    @test serial(ats[10]) == 113
    @test isa(ats[10], DisorderedAtom)
    @test countatoms(DisorderedResidue[struc['A'][10], struc['A'][16]]) == 17
    res = collectresidues(DisorderedResidue[struc['A'][16], struc['A'][10]])
    @test length(res) == 2
    @test isa(res, Vector{DisorderedResidue})
    @test resnumber(res[1]) == 16
    @test countresidues(DisorderedResidue[struc['A'][16], struc['A'][10]]) == 2

    # Test parsing 1SSU
    struc = read(testfilepath("mmCIF", "1SSU.cif"), MMCIFFormat)
    @test countmodels(struc) == 20
    @test modelnumbers(struc) == collect(1:20)
    @test countchains(struc) == 1
    @test countchains(struc[5]) == 1
    @test chainids(struc) == ["A"]
    @test chainids(struc[5]) == ["A"]
    # This is different to the PDB file as in the mmCIF file the atom serial
    #   does not reset for new models
    @test serial(struc[10]['A'][40]["HB2"]) == 7378
    @test countatoms(struc) == 756
    @test countatoms.([struc[i] for i in modelnumbers(struc)]) == 756 * ones(Int, 20)
    @test countatoms(struc, hydrogenselector) == 357
    ats = collectatoms(Model[struc[5], struc[10]])
    @test length(ats) == 1512
    # Due to the changes in atom serial this is also different to the PDB file
    @test z(ats[20]) == -14.782
    @test z(ats[1000]) == -3.367
    @test countatoms(Model[struc[5], struc[10]]) == 1512
    res = collectresidues(Model[struc[5], struc[10]])
    @test length(res) == 102
    @test y(res[10]["O"]) == -1.612
    @test y(res[100]["O"]) == -13.184
    @test countresidues(Model[struc[5], struc[10]]) == 102
    chs = collectchains(Model[struc[5], struc[10]])
    @test length(chs) == 2
    @test chainid.(chs) == ["A", "A"]
    @test z(chs[2][5]["CA"]) == -5.667
    @test countchains(Model[struc[5], struc[10]]) == 2
    mos = collectmodels(Model[struc[10], struc[5]])
    @test length(mos) == 2
    @test modelnumber.(mos) == [10, 5]
    @test z(mos[2]['A'][5]["CA"]) == -5.837
    @test countmodels(Model[struc[10], struc[5]]) == 2

    # Test parsing multi-character chain IDs
    multichar_str = """
        data_test
        loop_
        _atom_site.group_PDB
        _atom_site.id
        _atom_site.type_symbol
        _atom_site.label_atom_id
        _atom_site.label_alt_id
        _atom_site.label_comp_id
        _atom_site.label_asym_id
        _atom_site.label_entity_id
        _atom_site.label_seq_id
        _atom_site.pdbx_PDB_ins_code
        _atom_site.Cartn_x
        _atom_site.Cartn_y
        _atom_site.Cartn_z
        _atom_site.occupancy
        _atom_site.B_iso_or_equiv
        _atom_site.pdbx_formal_charge
        _atom_site.auth_seq_id
        _atom_site.auth_comp_id
        _atom_site.auth_asym_id
        _atom_site.auth_atom_id
        _atom_site.pdbx_PDB_model_num
        ATOM   1    N N   . MET A 1 1   ? 26.981 53.977  40.085 1.00 40.83  ? 1   MET A1 N   1
        ATOM   2    C CA  . MET A 1 1   ? 26.091 52.849  39.889 1.00 37.14  ? 1   MET A1 CA  1
        ATOM   3    C C   . MET A 1 1   ? 26.679 52.163  38.675 1.00 30.15  ? 1   MET A2 C   1
        ATOM   4    O O   . MET A 1 1   ? 27.020 52.865  37.715 1.00 27.59  ? 1   MET A2 O   1
        """
    struc = read(IOBuffer(multichar_str), MMCIFFormat)
    @test chainids(struc) == ["A1", "A2"]

    # Test parsing multi-line construct to structure
    multlinestruc_str = """
        data_test
        loop_
        _atom_site.group_PDB
        _atom_site.id
        _atom_site.type_symbol
        _atom_site.label_atom_id
        _atom_site.label_alt_id
        _atom_site.label_comp_id
        _atom_site.label_asym_id
        _atom_site.label_entity_id
        _atom_site.label_seq_id
        _atom_site.pdbx_PDB_ins_code
        _atom_site.Cartn_x
        _atom_site.Cartn_y
        _atom_site.Cartn_z
        _atom_site.occupancy
        _atom_site.B_iso_or_equiv
        _atom_site.pdbx_formal_charge
        _atom_site.auth_seq_id
        _atom_site.auth_comp_id
        _atom_site.auth_asym_id
        _atom_site.auth_atom_id
        _atom_site.pdbx_PDB_model_num
        ATOM   1    N N   . MET A 1 1   ? 26.981 53.977  40.085 1.00 40.83  ? 1   MET A N   1
        ATOM   2    C CA  . MET A 1 1   ? 26.091
        ;52.849
        ;
        39.889 1.00 37.14  ? 1   MET A CA  1
        ATOM   3    C C   . MET A 1 1   ? 26.679 52.163  38.675 1.00 30.15  ? 1   MET A C   1
        ATOM   4    O O   . MET A 1 1   ? 27.020 52.865  37.715 1.00 27.59  ? 1   MET A O   1
        """
    struc = read(IOBuffer(multlinestruc_str), MMCIFFormat)
    @test coords(struc['A'][1]["CA"]) == [26.091, 52.849, 39.889]
    @test serial(struc['A'][1]["O"]) == 4

    # Test files that should not parse
    @test_throws Exception read(testfilepath("mmCIF", "1AKE_err.cif"), MMCIFFormat)
    @test_throws ErrorException read(testfilepath("mmCIF", "1EN2_err.cif"), MMCIFFormat)

    # Test parsing empty file
    struc = read(IOBuffer(""), MMCIFFormat)
    @test isa(struc, MolecularStructure)
    @test countmodels(struc) == 0

    # Test formatting
    @test !requiresnewline("foobar")
    @test requiresnewline("foo\nbar")
    @test !requiresnewline("foo' bar")
    @test requiresnewline("foo' b\" ar")

    @test !requiresquote("foobar")
    @test requiresquote("foo bar")
    @test requiresquote("_foobar")
    @test requiresquote("data_foobar")
    @test requiresquote("global_")

    @test formatmmcifcol("foo") == "foo"
    @test formatmmcifcol("foo\nbar", 5) == "\n;foo\nbar\n;\n"
    @test formatmmcifcol("_foobar", 10) == "'_foobar' "
    @test formatmmcifcol("foo' bar", 12) == "\"foo' bar\"  "

    # Test writemmcif
    #
    # Note: We test everything for the uncompressed as well as the
    # gzip-compressed case.  Because of this, every call to the mmcif
    # functions has a trailing `; kwargs...`.  We also have to use a
    # custom `countlines_gzip(filename; gzip=false)` function that
    # handles the case that `filename` is compressed.
    for kwargs in [(), (gzip=true,)]
        dic_one = MMCIFDict(testfilepath("mmCIF", "1AKE.cif"))
        writemmcif(temp_filename, dic_one; kwargs...)
        dic_two = MMCIFDict(temp_filename; kwargs...)
        @test length(keys(dic_one)) == length(keys(dic_two))
        @test all([haskey(dic_two, k) for k in keys(dic_one)])
        @test all([dic_one[k] == dic_two[k] for k in keys(dic_one)])

        writemmcif(temp_filename, MMCIFDict(Dict("key.key" => ["value"])); kwargs...)
        @test_throws ArgumentError writemmcif(temp_filename, MMCIFDict(Dict("key" => ["value"]));
                                              kwargs...)
        @test_throws ArgumentError writemmcif(temp_filename, MMCIFDict(Dict("key.key.key" => ["value"]));
                                              kwargs...)
        @test_throws ArgumentError writemmcif(temp_filename, MMCIFDict(Dict(
            "key.one" => ["value"],
            "key.two" => ["value1", "value2"],
        )); kwargs...)

        struc = read(testfilepath("PDB", "1SSU.pdb"), PDBFormat)
        writemmcif(temp_filename, struc; kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 15145
        dic_written = MMCIFDict(temp_filename; kwargs...)
        @test length(keys(dic_written)) == 22
        @test length(dic_written["_atom_site.group_PDB"]) == 15120
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test isa(struc_written, MolecularStructure)
        @test modelnumbers(struc_written) == collect(1:20)
        @test countatoms(struc_written) == 756
        @test z(struc_written[4]['A']["30"]["OG"]) == -2.177
        @test atomnames(struc_written[15]['A']["39"]) == [
            "N", "CA", "C", "O", "CB", "SG", "H", "HA", "HB2", "HB3"]

        writemmcif(temp_filename, MMCIFDict(); kwargs...)
        dic_written = MMCIFDict(temp_filename; kwargs...)
        @test length(keys(dic_written)) == 1
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test countatoms(struc_written) == 0

        # Test writing to stream
        open(temp_filename, "w") do file
            writemmcif(file, struc; kwargs...)
        end
        @test countlines_gzip(temp_filename; kwargs...) == 15145

        # Test selectors
        struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat)
        writemmcif(temp_filename, struc, heteroselector; kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 524
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test modelnumbers(struc_written) == [1]
        @test countatoms(struc_written) == 492
        @test chainids(struc_written) == ["A", "B"]
        @test tempfactor(struc_written['B']["H_705"]["O"]) == 64.17
        writemmcif(temp_filename, collectatoms(struc, standardselector,
                                               disorderselector);
                   kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 35
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test countatoms(struc_written) == 5
        @test sum(isdisorderedatom, collectatoms(struc_written)) == 5
        @test defaultaltlocid(struc_written['A'][167]["NH1"]) == 'A'
        writemmcif(temp_filename, struc; kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 3841
        writemmcif(temp_filename, struc, expand_disordered=false;
                   kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 3829
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test !any(isdisorderedatom.(collectatoms(struc_written)))

        # Test writing different element types
        writemmcif(temp_filename, struc[1]; kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 3841
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test modelnumbers(struc_written) == [1]
        @test countatoms(struc_written) == 3804
        writemmcif(temp_filename, struc['A']; kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 1991
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test chainids(struc_written) == ["A"]
        writemmcif(temp_filename, struc['A'][50]; kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 34
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test chainids(struc_written) == ["A"]
        @test countresidues(struc_written) == 1
        @test atomnames(struc_written['A'][50]) == [
            "N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"]
        writemmcif(temp_filename, struc['A'][50]["CA"]; kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 24
        writemmcif(temp_filename, struc['A'][167]["CZ"]; kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 27
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test countatoms(struc_written) == 1
        @test isdisorderedatom(collectatoms(struc_written)[1])
        @test length(collectatoms(struc_written)[1]) == 2
        @test tempfactor(collectatoms(struc_written)[1]) == 16.77
        writemmcif(temp_filename, Chain[struc['A'], struc['B']];
                   kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 3841
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test chainids(struc_written) == ["A", "B"]
        @test countatoms(struc_written['A']) == 1954
        @test countatoms(struc_written['B']) == 1850
        @test altlocids(struc_written['A']["H_215"]["O1G"]) == ['A', 'B']
        writemmcif(temp_filename, AbstractResidue[struc['A'][51], struc['A'][50]]; kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 42
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test countresidues(struc_written) == 2
        @test resid.(collectresidues(struc_written)) == ["50", "51"]
        @test countatoms(struc_written) == 17
        writemmcif(temp_filename, AbstractAtom[struc['A'][51]["CA"], struc['A'][50]["CA"]]; kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 27
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test countatoms(struc_written) == 2
        @test !ishetero(struc_written['A'][51]["CA"])
        writemmcif(temp_filename, struc['A'][51]["N"], calphaselector;
                   kwargs...)
        dic_written = MMCIFDict(temp_filename; kwargs...)
        @test length(keys(dic_written)) == 1
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test countatoms(struc_written) == 0
        writemmcif(temp_filename, MolecularStructure(); kwargs...)
        dic_written = MMCIFDict(temp_filename; kwargs...)
        @test length(keys(dic_written)) == 1
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test countatoms(struc_written) == 0

        # Test multiple model writing
        struc = read(testfilepath("PDB", "1SSU.pdb"), PDBFormat)
        writemmcif(temp_filename, Model[struc[10], struc[5]]; kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 1537
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test modelnumbers(struc_written) == [5, 10]
        @test modelnumber(defaultmodel(struc_written)) == 5
        @test countatoms(struc_written[5]) == 756
        @test countatoms(struc_written[10]) == 756
        @test_throws KeyError struc_written[1]

        # Test disordered residue writing
        struc = read(testfilepath("PDB", "1EN2.pdb"), PDBFormat)
        writemmcif(temp_filename, struc; kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 844
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test countatoms(struc_written) == 754
        @test isa(struc_written['A'][15], Residue)
        @test isa(struc_written['A'][16], DisorderedResidue)
        @test defaultresname(struc_written['A'][16]) == "ARG"
        @test isa(struc_written['A'][16]["N"], DisorderedAtom)
        @test defaultaltlocid(struc_written['A'][16]["N"]) == 'A'
        @test isa(disorderedres(struc_written['A'][16], "TRP")["N"], Atom)
        @test countatoms(struc_written['A'][16]) == 11
        @test countatoms(disorderedres(struc_written['A'][16], "TRP")) == 14
        writemmcif(temp_filename, AbstractResidue[struc['A'][16], struc['A'][10]]; kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 71
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test countresidues(struc_written) == 2
        @test isa(struc_written['A'][10], DisorderedResidue)
        @test isa(struc_written['A'][16], DisorderedResidue)
        @test defaultresname(struc_written['A'][10]) == "SER"
        @test isa(disorderedres(struc_written['A'][10], "GLY")["O"], Atom)
        @test altlocid(disorderedres(struc_written['A'][10], "GLY")["O"]) == 'B'
        @test countatoms(struc_written['A'][10]) == 6
        @test countatoms(struc_written['A'][16]) == 11
        writemmcif(temp_filename, struc, expand_disordered=false; kwargs...)
        @test countlines_gzip(temp_filename; kwargs...) == 779
        struc_written = read(temp_filename, MMCIFFormat; kwargs...)
        @test !any(isdisorderedres.(collectresidues(struc_written)))
        @test !any(isdisorderedatom.(collectatoms(struc_written)))
        @test countresidues(struc_written, expand_disordered=true) == 166
        @test countatoms(struc_written, expand_disordered=true) == 754

        # Test writing multi-character chain IDs
        struc = read(IOBuffer(multichar_str), MMCIFFormat)
        writemmcif(temp_filename, struc; kwargs...)
        struc_back = read(temp_filename, MMCIFFormat; kwargs...)
        @test chainids(struc_back) == ["A1", "A2"]
    end

    # Test readmultimmcif
    # Read from string
    comment_str = """
        data_test
        _test_single foo # Ignore this comment
        loop_
        _test_loop_1
        _test_loop_2
        a b # Ignore this comment
        c d
        """
    underscore_str = """
        data_4Q9R
        loop_
        _pdbx_audit_revision_item.ordinal
        _pdbx_audit_revision_item.revision_ordinal
        _pdbx_audit_revision_item.data_content_type
        _pdbx_audit_revision_item.item
        1  5 'Structure model' '_atom_site.B_iso_or_equiv'
        2  5 'Structure model' '_atom_site.Cartn_x'
        3  5 'Structure model' '_atom_site.Cartn_y'
        4 5 'Structure model' '_atom_site.Cartn_z'
        """
    cifs = readmultimmcif(IOBuffer(comment_str * underscore_str))
    dic = MMCIFDict(IOBuffer(comment_str))
    @test cifs["test"]["_test_single"] == ["foo"]
    @test cifs["test"]["_test_loop_1"] == ["a", "c"]
    @test cifs["test"]["_test_loop_2"] == ["b", "d"]
    @test length(keys(cifs["4Q9R"])) == 5
    @test cifs["4Q9R"]["_pdbx_audit_revision_item.item"] == [
        "_atom_site.B_iso_or_equiv", "_atom_site.Cartn_x",
        "_atom_site.Cartn_y", "_atom_site.Cartn_z"
    ]
    # Read from file - write a small multicif manually and read it back
    for gzip in (false, true)
        transcoder = gzip ? (GzipCompressorStream,) : ()
        open(transcoder..., temp_filename, "w") do f_out
            open(testfilepath("mmCIF", "1AKE.cif")) do f_in
                write(f_out, f_in)
            end
            open(testfilepath("mmCIF", "1EN2.cif")) do f_in
                write(f_out, f_in)
            end
        end
        cifs = readmultimmcif(temp_filename; gzip=gzip)
        # Tests for 1AKE
        dic = cifs["1AKE"]
        @test isa(dic.dict, Dict{String, Vector{String}})
        @test dic["_pdbx_database_status.recvd_initial_deposition_date"] == ["1991-11-08"]
        @test dic["_audit_author.name"] == ["Mueller, C.W.", "Schulz, G.E."]
        @test length(dic["_atom_site.group_PDB"]) == 3816
        dic["_pdbx_database_status.recvd_initial_deposition_date"] = ["changed"]
        @test dic["_pdbx_database_status.recvd_initial_deposition_date"] == ["changed"]
        @test length(keys(dic)) == 610
        @test length(values(dic)) == 610
        @test haskey(dic, "_cell.entry_id")
        @test !haskey(dic, "nokey")
        @test get(dic, "_cell.entry_id", ["default"]) == ["1AKE"]
        @test get(dic, "nokey", ["default"]) == ["default"]
        @test ismissing(get(dic, "nokey", missing))
        # Tests for 1EN2
        struc = MolecularStructure(cifs["1EN2"])
        @test modelnumbers(struc) == [1]
        @test chainids(struc[1]) == ["A"]
        @test serial(struc['A'][48]["CA"]) == 394
        @test isdisorderedres(struc['A'][10])
        @test defaultresname(struc['A'][10]) == "SER"
        @test resname(disorderedres(struc['A'][10], "SER")) == "SER"
        @test resname(disorderedres(struc['A'][10], "GLY")) == "GLY"
        @test !isdisorderedatom(struc['A'][10]["CA"])
        @test altlocid(struc['A'][10]["CA"]) == 'A'
        @test isdisorderedres(struc['A'][16])
        @test defaultresname(struc['A'][16]) == "ARG"
        @test resname(disorderedres(struc['A'][16], "ARG")) == "ARG"
        @test resname(disorderedres(struc['A'][16], "TRP")) == "TRP"
        @test !isdisorderedatom(disorderedres(struc['A'][16], "TRP")["CA"])
        @test altlocid(disorderedres(struc['A'][16], "TRP")["CA"]) == 'C'
        @test isdisorderedatom(struc['A'][16]["CA"])
        @test altlocids(struc['A'][16]["CA"]) == ['A', 'B']
        @test defaultaltlocid(struc['A'][16]["CA"]) == 'A'
        @test occupancy(struc['A'][16]["CA"]) == 0.22
        ats = collectatoms(DisorderedResidue[struc['A'][10], struc['A'][16]])
        @test length(ats) == 17
        @test isa(ats, Vector{AbstractAtom})
        @test serial(ats[10]) == 113
        @test isa(ats[10], DisorderedAtom)
        @test countatoms(DisorderedResidue[struc['A'][10], struc['A'][16]]) == 17
        res = collectresidues(DisorderedResidue[struc['A'][16], struc['A'][10]])
        @test length(res) == 2
        @test isa(res, Vector{DisorderedResidue})
        @test resnumber(res[1]) == 16
        @test countresidues(DisorderedResidue[struc['A'][16], struc['A'][10]]) == 2
    end
    # Test error on duplicate data_ entry somewhere in the middle
    test_multicif = """
        data_test
        _test_single foo # Ignore this comment
        data_test
        _test_single foo # Ignore this comment
        data_foo
        _test_single foo # Ignore this comment
    """
    @test_throws ErrorException readmultimmcif(IOBuffer(test_multicif))
    # Test error on duplicate data_ entry at the end
    test_multicif = """
        data_test
        _test_single foo # Ignore this comment
        data_test
        _test_single foo # Ignore this comment
    """
    @test_throws ErrorException readmultimmcif(IOBuffer(test_multicif))

    # Test writemultimmcif
    test_multicif = Dict{String, MMCIFDict}()
    for pdbid in ("1AKE", "1SSU")
        cif = MMCIFDict(testfilepath("mmCIF", "$pdbid.cif"))
        test_multicif[pdbid] = cif
    end
    # Write to buffer
    for gzip in (false, true)
        open(temp_filename, "w") do io
            writemultimmcif(io, test_multicif; gzip=gzip)
        end
        cifs = readmultimmcif(temp_filename; gzip=gzip)
        for k in keys(cifs)
            @test cifs[k].dict == test_multicif[k].dict
        end
    end
    # Write to filepath
    for gzip in (false, true)
        writemultimmcif(temp_filename, test_multicif; gzip=gzip)
        cifs = readmultimmcif(temp_filename; gzip=gzip)
        for k in keys(cifs)
            @test cifs[k].dict == test_multicif[k].dict
        end
    end
    # test warning
    @test_logs (:warn,
                "writemultimmcif: MMCIFDict for key \"not_1AKE\" has different \"data_\" key (\"1AKE\")"
                ) writemultimmcif(temp_filename, Dict("not_1AKE" => test_multicif["1AKE"]))
end

@testset "MMTF" begin
    # Test MMTF dictionary
    dic = MMTFDict(Dict())
    show(IOBuffer(), MIME("text/plain"), MMTFDict())
    show(IOBuffer(), MIME("text/plain"), MMTFDict(Dict("a" => "b")))

    dic = MMTFDict(testfilepath("MMTF", "1AKE.mmtf"))
    @test isa(dic.dict, Dict{String, Any})
    @test dic["rWork"] == 0.196f0
    @test dic["chainNameList"] == ["A", "B", "A", "B", "A", "B"]
    @test length(dic["groupIdList"]) == 808
    @test length(keys(dic)) == 39
    @test length(values(dic)) == 39
    @test haskey(dic, "chainNameList")
    @test !haskey(dic, "nokey")
    @test get(dic, "chainNameList", ["default"]) == ["A", "B", "A", "B", "A", "B"]
    @test get(dic, "nokey", ["default"]) == ["default"]
    @test ismissing(get(dic, "nokey", missing))
    show(devnull, dic)
    dic_gzip = MMTFDict(testfilepath("MMTF", "1AKE.mmtf.gz"), gzip=true)
    @test dic.dict == dic_gzip.dict
    open(testfilepath("MMTF", "1AKE.mmtf.gz")) do file
        dic_gzip = MMTFDict(file, gzip=true)
        @test dic.dict == dic_gzip.dict
    end
    dic_red = MMTFDict(testfilepath("MMTF", "1AKE_reduced.mmtf"))
    @test dic_red["numAtoms"] == 549
    @test dic_red["entityList"] == dic["entityList"]
    dic["numGroups"] = 100
    @test dic["numGroups"] == 100

    # Test parsing 1AKE
    struc = read(testfilepath("MMTF", "1AKE.mmtf"), MMTFFormat)
    @test structurename(struc) == "1AKE.mmtf"
    @test countmodels(struc) == 1
    @test modelnumbers(struc) == [1]
    @test countchains(struc) == 2
    @test countchains(struc[1]) == 2
    @test chainids(struc) == ["A", "B"]
    @test chainids(struc[1]) == ["A", "B"]
    @test resname(struc['A'][10]) == "GLY"
    @test !isdisorderedres(struc['A'][10])
    @test serial(struc['A'][200]["NZ"]) == 1555
    @test !isdisorderedatom(struc['A'][200]["NZ"])
    @test isdisorderedatom(struc['A'][167]["CD"])
    @test altlocids(struc['A'][167]["CD"]) == ['A', 'B']
    @test isapprox(x(struc['A'][167]["CD"]), 24.502, atol=1e-5)
    @test isapprox(x(struc['A'][167]["CD"]['A']), 24.502, atol=1e-5)
    @test isapprox(x(struc['A'][167]["CD"]['B']), 24.69, atol=1e-5)
    mos = collect(struc)
    @test modelnumber(mos[1]) == 1
    chs = collect(mos[1])
    @test chainid.(chs) == ["A", "B"]
    res = collect(chs[1])
    @test length(res) == 456
    @test resid(res[20]) == "20"
    ats = collect(res[20])
    @test length(ats) == 8
    @test atomname.(ats) == ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"]

    struc = read(testfilepath("MMTF", "1AKE_reduced.mmtf"), MMTFFormat)
    @test countatoms(struc) == 542
    @test countatoms(struc, standardselector, !calphaselector) == 0
    @test chainids(struc) == ["A", "B"]
    @test isapprox(x(struc['A'][167]["CA"]), 23.859, atol=1e-5)
    @test isapprox(tempfactor(struc['A'][167]["CA"]), 15.89, atol=1e-5)

    # Test parsing options
    struc = read(testfilepath("MMTF", "1AKE.mmtf"), MMTFFormat, structure_name="New name")
    @test structurename(struc) == "New name"
    @test countatoms(struc) == 3804

    struc = read(testfilepath("MMTF", "1AKE.mmtf"), MMTFFormat, read_het_atoms=false)
    @test countatoms(struc) == 3312
    # Different to the PDB file due to the lack of TER label serial
    @test serial(collectatoms(struc)[2000]) == 2005
    @test sum(ishetero, collectatoms(struc)) == 0

    struc = read(testfilepath("MMTF", "1AKE.mmtf"), MMTFFormat, read_std_atoms=false)
    @test countatoms(struc) == 492
    @test serial(collectatoms(struc)[400]) == 3724
    @test sum(ishetero, collectatoms(struc)) == 492

    struc = read(testfilepath("MMTF", "1AKE.mmtf"), MMTFFormat, read_het_atoms=false,
                read_std_atoms=false)
    @test countatoms(struc) == 0
    @test countresidues(struc) == 0
    @test countchains(struc) == 0
    @test countmodels(struc) == 0

    struc = read(testfilepath("MMTF", "1AKE.mmtf"), MMTFFormat, remove_disorder=true)
    @test countatoms(struc) == 3804
    @test sum(isdisorderedatom, collectatoms(struc)) == 0
    @test isapprox(tempfactor(struc['A'][167]["NE"]), 23.32, atol=1e-5)

    struc = read(testfilepath("MMTF", "1AKE.mmtf.gz"), MMTFFormat, gzip=true)
    @test chainids(struc) == ["A", "B"]
    @test serial(struc['A'][200]["NZ"]) == 1555
    @test isapprox(x(struc['A'][167]["CD"]['A']), 24.502, atol=1e-5)

    # Test parsing from stream
    open(testfilepath("MMTF", "1AKE.mmtf")) do file
        struc = read(file, MMTFFormat)
        @test countatoms(struc) == 3804
        @test countresidues(struc) == 808
    end
    open(testfilepath("MMTF", "1AKE.mmtf.gz")) do file
        struc = read(file, MMTFFormat, gzip=true)
        @test countatoms(struc) == 3804
        @test countresidues(struc) == 808
    end

    # Test parsing 1EN2
    struc = read(testfilepath("MMTF", "1EN2.mmtf"), MMTFFormat)
    @test modelnumbers(struc) == [1]
    @test chainids(struc[1]) == ["A"]
    @test serial(struc['A'][48]["CA"]) == 394
    @test isdisorderedres(struc['A'][10])
    @test defaultresname(struc['A'][10]) == "SER"
    @test resname(disorderedres(struc['A'][10], "SER")) == "SER"
    @test resname(disorderedres(struc['A'][10], "GLY")) == "GLY"
    @test !isdisorderedatom(struc['A'][10]["CA"])
    @test altlocid(struc['A'][10]["CA"]) == 'A'
    @test isdisorderedres(struc['A'][16])
    @test defaultresname(struc['A'][16]) == "ARG"
    @test resname(disorderedres(struc['A'][16], "ARG")) == "ARG"
    @test resname(disorderedres(struc['A'][16], "TRP")) == "TRP"
    @test !isdisorderedatom(disorderedres(struc['A'][16], "TRP")["CA"])
    @test altlocid(disorderedres(struc['A'][16], "TRP")["CA"]) == 'C'
    @test isdisorderedatom(struc['A'][16]["CA"])
    @test altlocids(struc['A'][16]["CA"]) == ['A', 'B']
    @test defaultaltlocid(struc['A'][16]["CA"]) == 'A'
    @test isapprox(occupancy(struc['A'][16]["CA"]), 0.22, atol=1e-5)
    ats = collectatoms(DisorderedResidue[struc['A'][10], struc['A'][16]])
    @test length(ats) == 17
    @test isa(ats, Vector{AbstractAtom})
    @test serial(ats[10]) == 113
    @test isa(ats[10], DisorderedAtom)
    @test countatoms(DisorderedResidue[struc['A'][10], struc['A'][16]]) == 17
    res = collectresidues(DisorderedResidue[struc['A'][16], struc['A'][10]])
    @test length(res) == 2
    @test isa(res, Vector{DisorderedResidue})
    @test resnumber(res[1]) == 16
    @test countresidues(DisorderedResidue[struc['A'][16], struc['A'][10]]) == 2

    # Test parsing 1SSU
    struc = read(testfilepath("MMTF", "1SSU.mmtf"), MMTFFormat)
    @test countmodels(struc) == 20
    @test modelnumbers(struc) == collect(1:20)
    @test countchains(struc) == 1
    @test countchains(struc[5]) == 1
    @test chainids(struc) == ["A"]
    @test chainids(struc[5]) == ["A"]
    # This is different to the PDB file as in the MMTF file the atom serial
    #   does not reset for new models
    @test serial(struc[10]['A'][40]["HB2"]) == 7378
    @test countatoms(struc) == 756
    @test countatoms.([struc[i] for i in modelnumbers(struc)]) == 756 * ones(Int, 20)
    @test countatoms(struc, hydrogenselector) == 357
    ats = collectatoms(Model[struc[5], struc[10]])
    @test length(ats) == 1512
    # Due to the changes in atom serial this is also different to the PDB file
    @test isapprox(z(ats[20]), -14.782, atol=1e-5)
    @test isapprox(z(ats[1000]), -3.367, atol=1e-5)
    @test countatoms(Model[struc[5], struc[10]]) == 1512
    res = collectresidues(Model[struc[5], struc[10]])
    @test length(res) == 102
    @test isapprox(y(res[10]["O"]), -1.612, atol=1e-5)
    @test isapprox(y(res[100]["O"]), -13.184, atol=1e-5)
    @test countresidues(Model[struc[5], struc[10]]) == 102
    chs = collectchains(Model[struc[5], struc[10]])
    @test length(chs) == 2
    @test chainid.(chs) == ["A", "A"]
    @test isapprox(z(chs[2][5]["CA"]), -5.667, atol=1e-5)
    @test countchains(Model[struc[5], struc[10]]) == 2
    mos = collectmodels(Model[struc[10], struc[5]])
    @test length(mos) == 2
    @test modelnumber.(mos) == [10, 5]
    @test isapprox(z(mos[2]['A'][5]["CA"]), -5.837, atol=1e-5)
    @test countmodels(Model[struc[10], struc[5]]) == 2

    # Test writemmtf
    dic_one = MMTFDict(testfilepath("MMTF", "1AKE.mmtf"))
    writemmtf(temp_filename, dic_one)
    dic_two = MMTFDict(temp_filename)
    @test length(keys(dic_one)) == length(keys(dic_two))
    @test all([haskey(dic_two, k) for k in keys(dic_one)])
    @test all([dic_one[k] == dic_two[k] for k in keys(dic_one)])

    dic = MMTFDict(testfilepath("MMTF", "1SSU.mmtf"))
    struc = read(testfilepath("MMTF", "1SSU.mmtf"), MMTFFormat)
    for gzip in (false, true)
        writemmtf(temp_filename, struc, gzip=gzip)
        dic_written = MMTFDict(temp_filename, gzip=gzip)
        for k in ["chainsPerModel", "atomIdList", "xCoordList", "yCoordList",
                    "zCoordList", "altLocList", "occupancyList", "bFactorList",
                    "insCodeList", "groupIdList", "chainIdList", "chainNameList",
                    "groupsPerChain"]
            @test dic[k] == dic_written[k]
        end
        for k in ["groupTypeList", "groupList"]
            @test length(dic[k]) == length(dic_written[k])
        end
    end

    for (ft, dir_name) in ((PDBFormat, "PDB"), (MMCIFFormat, "mmCIF"), (MMTFFormat, "MMTF"))
        struc = read(testfilepath(dir_name, "1SSU.$(pdbextension[ft])"), ft)
        writemmtf(temp_filename, struc)
        struc_written = read(temp_filename, MMTFFormat)
        @test isa(struc_written, MolecularStructure)
        @test modelnumbers(struc_written) == collect(1:20)
        @test countatoms(struc_written) == 756
        @test isapprox(z(struc_written[4]['A']["30"]["OG"]), -2.177, atol=1e-5)
        @test atomnames(struc_written[15]['A']["39"]) == [
            "N", "CA", "C", "O", "CB", "SG", "H", "HA", "HB2", "HB3"]
    end

    writemmtf(temp_filename, MMTFDict())
    dic_written = MMTFDict(temp_filename)
    @test length(keys(dic_written)) == 39
    struc_written = read(temp_filename, MMTFFormat)
    @test countatoms(struc_written) == 0

    struc_red = read(testfilepath("MMTF", "1AKE_reduced.mmtf"), MMTFFormat)
    writemmtf(temp_filename, struc_red)
    struc_written = read(temp_filename, MMTFFormat)
    @test countatoms(struc_written) == 542
    @test countatoms(struc_written, standardselector, !calphaselector) == 0
    @test chainids(struc_written) == ["A", "B"]
    dic_written = MMTFDict(temp_filename)
    for group in dic_written["groupList"]
        @test length(group["atomNameList"]) == length(group["elementList"]) == length(group["formalChargeList"])
    end

    # Test writing to stream
    for gzip in (false, true)
        open(temp_filename, "w") do file
            writemmtf(file, struc, gzip=gzip)
        end
        struc_written = read(temp_filename, MMTFFormat, gzip=gzip)
        @test countatoms(struc_written) == 756
    end

    # Test selectors
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat)
    writemmtf(temp_filename, struc, heteroselector)
    struc_written = read(temp_filename, MMTFFormat)
    @test modelnumbers(struc_written) == [1]
    @test countatoms(struc_written) == 492
    @test chainids(struc_written) == ["A", "B"]
    @test isapprox(tempfactor(struc_written['B']["H_705"]["O"]), 64.17, atol=1e-5)
    writemmtf(temp_filename, collectatoms(struc, standardselector,
                                                disorderselector))
    struc_written = read(temp_filename, MMTFFormat)
    @test countatoms(struc_written) == 5
    @test sum(isdisorderedatom, collectatoms(struc_written)) == 5
    @test defaultaltlocid(struc_written['A'][167]["NH1"]) == 'A'
    writemmtf(temp_filename, struc, expand_disordered=false)
    struc_written = read(temp_filename, MMTFFormat)
    @test !any(isdisorderedatom.(collectatoms(struc_written)))
    @test countatoms(struc_written) == 3804

    # Test writing different element types
    writemmtf(temp_filename, struc[1])
    struc_written = read(temp_filename, MMTFFormat)
    @test modelnumbers(struc_written) == [1]
    @test countatoms(struc_written) == 3804
    @test countatoms(struc_written, expand_disordered=true) == 3816
    writemmtf(temp_filename, struc['A'])
    struc_written = read(temp_filename, MMTFFormat)
    @test countatoms(struc_written) == 1954
    @test chainids(struc_written) == ["A"]
    writemmtf(temp_filename, struc['A'][50])
    struc_written = read(temp_filename, MMTFFormat)
    @test chainids(struc_written) == ["A"]
    @test countresidues(struc_written) == 1
    @test countatoms(struc_written) == 9
    @test atomnames(struc_written['A'][50]) == [
        "N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"]
    writemmtf(temp_filename, struc['A'][50]["CA"])
    struc_written = read(temp_filename, MMTFFormat)
    @test countatoms(struc_written) == 1
    writemmtf(temp_filename, struc['A'][167]["CZ"])
    struc_written = read(temp_filename, MMTFFormat)
    @test countatoms(struc_written) == 1
    @test countatoms(struc_written, expand_disordered=true) == 2
    @test isdisorderedatom(collectatoms(struc_written)[1])
    @test length(collectatoms(struc_written)[1]) == 2
    @test isapprox(tempfactor(collectatoms(struc_written)[1]), 16.77, atol=1e-5)
    writemmtf(temp_filename, Chain[struc['A'], struc['B']])
    struc_written = read(temp_filename, MMTFFormat)
    @test chainids(struc_written) == ["A", "B"]
    @test countatoms(struc_written['A']) == 1954
    @test countatoms(struc_written['B']) == 1850
    @test altlocids(struc_written['A']["H_215"]["O1G"]) == ['A', 'B']
    writemmtf(temp_filename, AbstractResidue[struc['A'][51], struc['A'][50]])
    struc_written = read(temp_filename, MMTFFormat)
    @test countresidues(struc_written) == 2
    @test resid.(collectresidues(struc_written)) == ["50", "51"]
    @test countatoms(struc_written) == 17
    writemmtf(temp_filename, AbstractAtom[struc['A'][51]["CA"], struc['A'][50]["CA"]])
    struc_written = read(temp_filename, MMTFFormat)
    @test countatoms(struc_written) == 2
    @test !ishetero(struc_written['A'][51]["CA"])
    writemmtf(temp_filename, struc['A'][51]["N"], calphaselector)
    dic_written = MMTFDict(temp_filename)
    @test length(keys(dic_written)) == 39
    struc_written = read(temp_filename, MMTFFormat)
    @test countatoms(struc_written) == 0
    writemmtf(temp_filename, MolecularStructure())
    dic_written = MMTFDict(temp_filename)
    @test length(keys(dic_written)) == 39
    struc_written = read(temp_filename, MMTFFormat)
    @test countatoms(struc_written) == 0

    # Test multiple model writing
    struc = read(testfilepath("PDB", "1SSU.pdb"), PDBFormat)
    writemmtf(temp_filename, Model[struc[10], struc[5]])
    struc_written = read(temp_filename, MMTFFormat)
    # Differs from PDB and mmCIF as model numbers are not recorded by MMTF
    @test modelnumbers(struc_written) == [1, 2]
    @test modelnumber(defaultmodel(struc_written)) == 1
    @test countatoms(struc_written[1]) == 756
    @test countatoms(struc_written[2]) == 756
    @test_throws KeyError struc_written[5]

    # Test disordered residue writing
    struc = read(testfilepath("PDB", "1EN2.pdb"), PDBFormat)
    writemmtf(temp_filename, struc)
    struc_written = read(temp_filename, MMTFFormat)
    @test countatoms(struc_written) == 754
    @test countatoms(struc_written, expand_disordered=true) == 819
    @test isa(struc_written['A'][15], Residue)
    @test isa(struc_written['A'][16], DisorderedResidue)
    @test defaultresname(struc_written['A'][16]) == "ARG"
    @test isa(struc_written['A'][16]["N"], DisorderedAtom)
    @test defaultaltlocid(struc_written['A'][16]["N"]) == 'A'
    @test isa(disorderedres(struc_written['A'][16], "TRP")["N"], Atom)
    @test countatoms(struc_written['A'][16]) == 11
    @test countatoms(disorderedres(struc_written['A'][16], "TRP")) == 14
    writemmtf(temp_filename, AbstractResidue[struc['A'][16], struc['A'][10]])
    struc_written = read(temp_filename, MMTFFormat)
    @test countresidues(struc_written) == 2
    @test countresidues(struc_written, expand_disordered=true) == 4
    @test isa(struc_written['A'][10], DisorderedResidue)
    @test isa(struc_written['A'][16], DisorderedResidue)
    @test defaultresname(struc_written['A'][10]) == "SER"
    @test isa(disorderedres(struc_written['A'][10], "GLY")["O"], Atom)
    @test altlocid(disorderedres(struc_written['A'][10], "GLY")["O"]) == 'B'
    @test countatoms(struc_written['A'][10]) == 6
    @test countatoms(struc_written['A'][16]) == 11
    writemmtf(temp_filename, struc, expand_disordered=false)
    struc_written = read(temp_filename, MMTFFormat)
    @test !any(isdisorderedres.(collectresidues(struc_written)))
    @test !any(isdisorderedatom.(collectatoms(struc_written)))
    @test countresidues(struc_written, expand_disordered=true) == 166
    @test countatoms(struc_written, expand_disordered=true) == 754

    # Test generatechainid
    @test generatechainid.([1, 20, 30, 18283]) == ["A", "T", "DA", "EAAA"]
    @test_throws ArgumentError generatechainid(0)
    @test_throws ArgumentError generatechainid(-10)
end

@testset "Spatial" begin
    # Test coordarray
    res = Residue("ALA", 1, ' ', false, Chain('A'))
    at = Atom(100, "CA", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res)
    cs = coordarray(at)
    @test size(cs) == (3, 1)
    @test cs[1] == 1.0
    @test cs[2] == 2.0
    @test cs[3] == 3.0

    struc_1AKE = read(testfilepath("PDB", "1AKE.pdb"), PDBFormat)
    cs = coordarray(struc_1AKE)
    @test size(cs) == (3, 3804)
    @test cs[1, 3787] == 20.135
    @test cs[2, 3787] == -10.789
    @test cs[3, 3787] == -1.732
    cs = coordarray(struc_1AKE['A'], calphaselector)
    @test size(cs) == (3, 214)
    @test cs[1, 10] == 17.487
    @test cs[2, 10] == 42.426
    @test cs[3, 10] == 19.756
    @test coordarray(cs) == cs
    cs = coordarray(struc_1AKE, expand_disordered=true)
    @test size(cs) == (3, 3816)

    # Test superimposition
    cs_one = [
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 0.0
    ]
    cs_two = [
        0.0 -1.0 0.0
        1.0  0.0 0.0
        1.0  1.0 1.0
    ]
    trans = Transformation(cs_one, cs_two)
    @test isapprox(trans.trans1, [1/3, 1/3, 0])
    @test isapprox(trans.trans2, [-1/3, 1/3, 1])
    rot_real = [
        0.0 -1.0 0.0
        1.0  0.0 0.0
        0.0  0.0 1.0
    ]
    @test isapprox(trans.rot, rot_real)

    # Test rot isn't a reflection
    cs_one = Float64[
        1 -1  0  0  0  0
        0  0  2 -2  0  0
        0  0  0  0  2 -2
    ]
    cs_two = Float64[
        -1  1  0  0  0  0
         0  0  2 -2  0  0
         0  0  0  0  2 -2
    ]
    trans = Transformation(cs_one, cs_two)
    @test isapprox(trans.trans1, [0.0, 0.0, 0.0])
    @test isapprox(trans.trans2, [0.0, 0.0, 0.0])
    @test isapprox(det(trans.rot), 1.0)
    rot_real = [
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
    ]
    @test isapprox(trans.rot, rot_real)

    cs_one = [
         7.0  5.0  3.0
         2.0  0.0 -2.0
        -1.0 -3.0 -5.0
    ]
    cs_two = [
        -1.5 -0.5 0.5
        -1.0  0.0 1.0
         2.0  3.0 4.0
    ]
    trans = Transformation(cs_one, cs_two)
    cs = applytransform(cs_one, trans)
    cs_real = [
        -2.5 -0.5 1.5
        -2.0  0.0 2.0
         1.0  3.0 5.0
    ]
    @test isapprox(cs, cs_real)

    ats = [
        Atom(100, "CA", ' ', [8.0, 3.0, 0.0], 1.0, 10.0, " C", "  ", res),
        Atom(101, "CB", ' ', [4.0, -1.0, -4.0], 1.0, 10.0, " C", "  ", res)
    ]
    applytransform!(ats, trans)
    cs_real = [
        -3.5 0.5
        -3.0 1.0
         0.0 4.0
    ]
    @test isapprox(coordarray(ats), cs_real)

    cs_one = [
        0.0 1.0       2.0
        0.0 sqrt(3.0) 0.0
    ]
    cs_two = [
        2.0           4.0 0.0
        2 * sqrt(3.0) 0.0 0.0
    ]
    trans = Transformation(cs_one, cs_two)
    @test isapprox(trans.trans1, [1.0, sqrt(3.0) / 3])
    @test isapprox(trans.trans2, [2.0, 2 * sqrt(3.0) / 3])
    rot_real = [
         cos(2 * pi / 3) sin(2 * pi / 3)
        -sin(2 * pi / 3) cos(2 * pi / 3)
    ]
    @test isapprox(trans.rot, rot_real)

    struc_1SSU = read(testfilepath("PDB", "1SSU.pdb"), PDBFormat)
    theta = 2 * pi * rand()
    rand_rot = [
        1.0 0.0         0.0
        0.0 cos(theta) -sin(theta)
        0.0 sin(theta)  cos(theta)
    ]
    rand_trans = Transformation(100 * randn(3), 100 * randn(3), rand_rot)
    applytransform!(struc_1SSU[1], rand_trans)
    superimpose!(struc_1SSU[1], struc_1SSU[2], standardselector)
    @test isapprox(coords(struc_1SSU[1]["A"][10]["CA"]),
                    [2.68601, -2.06401, 3.73024], atol=1e-5)
    @test coords(struc_1SSU[2]["A"][10]["CA"]) == [3.983, -3.252, 2.368]
    @test isapprox(rmsd(struc_1SSU[1], struc_1SSU[2], superimpose=false),
                    2.54859, atol=1e-5)
    applytransform!(struc_1SSU[1], rand_trans)
    superimpose!(struc_1SSU[1], struc_1SSU[2], standardselector,
                    alignatoms=cbetaselector)
    @test isapprox(rmsd(struc_1SSU[1], struc_1SSU[2], superimpose=false),
                    2.55015, atol=1e-5)
    trans = Transformation(collectresidues(struc_1SSU[3])[1:40],
                    collectresidues(struc_1SSU[2])[11:end],
                    standardselector)
    rot_real = [
        0.999633     0.00144653 -0.02704
        -0.00108653  0.999911    0.0133232
        0.0270569   -0.013289    0.999546
    ]
    @test isapprox(trans.rot, rot_real, atol=1e-5)
    @test trans.inds1 == collect(11:40)
    @test trans.inds2 == collect(1:30)
    superimpose!(collectresidues(struc_1SSU[3])[1:40],
                    collectresidues(struc_1SSU[2])[11:end],
                    standardselector)
    @test isapprox(rmsd(collectresidues(struc_1SSU[3])[21:30],
                    collectresidues(struc_1SSU[2])[21:30],
                    superimpose=false), 0.365745, atol=1e-5)

    # Test rmsd
    cs_one = [
        0.0 0.0
        1.0 0.0
        1.0 3.0
    ]
    cs_two = [
        0.0 2.0
        2.0 0.0
        1.0 3.0
    ]
    # This line gives a @simd warning when --inline=no
    @test isapprox(rmsd(cs_one, cs_two), sqrt(5/2))
    cs_one = [
        0.0 0.0 1.0
        1.0 0.0 2.0
        1.0 3.0 3.0
    ]
    cs_two = [
        0.0 2.0
        2.0 0.0
        1.0 3.0
    ]
    @test_throws ArgumentError rmsd(cs_one, cs_two)

    struc_1SSU = read(testfilepath("PDB", "1SSU.pdb"), PDBFormat)
    @test isapprox(rmsd(struc_1SSU[1], struc_1SSU[2], superimpose=false), 4.18219, atol=1e-5)
    @test isapprox(rmsd(struc_1SSU[1], struc_1SSU[2]), 2.54859, atol=1e-5)
    @test isapprox(rmsd(struc_1SSU[5], struc_1SSU[6], superimpose=false, rmsdatoms=backboneselector), 5.36997, atol=1e-5)
    @test isapprox(rmsd(struc_1SSU[5], struc_1SSU[6], rmsdatoms=backboneselector), 3.65486, atol=1e-5)
    @test_throws ArgumentError rmsd(struc_1SSU[1]['A'][8], struc_1SSU[1]['A'][9], superimpose=false, rmsdatoms=allselector)

    # Test displacements
    cs_one = [
        0.0 0.0
        1.0 0.0
        1.0 3.0
    ]
    cs_two = [
        0.0 2.0
        2.0 0.0
        1.0 4.0
    ]
    # This line gives a @simd warning when --inline=no
    @test isapprox(displacements(cs_one, cs_two), [1.0, sqrt(5)])
    cs_one = [
        0.0 0.0 1.0
        1.0 0.0 2.0
        1.0 3.0 3.0
    ]
    cs_two = [
        0.0 2.0
        2.0 0.0
        1.0 4.0
    ]
    @test_throws ArgumentError displacements(cs_one, cs_two)

    disps = displacements(struc_1SSU[5], struc_1SSU[10], superimpose=false, dispatoms=allselector)
    @test isa(disps, Vector{Float64})
    @test length(disps) == 756
    @test isapprox(disps[20], sqrt(1.984766))
    disps = displacements(struc_1SSU[5], struc_1SSU[10], dispatoms=allselector)
    @test isapprox(disps[20], 3.62754, atol=1e-5)
    disps = displacements(struc_1SSU[5], struc_1SSU[10], superimpose=false)
    @test length(disps) == 51
    @test isapprox(disps[20], sqrt(0.032822))
    disps = displacements(struc_1SSU[5], struc_1SSU[10])
    @test isapprox(disps[20], 0.717014, atol=1e-5)

    # Test sqdistance and distance
    at_a = Atom(100, "CA", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res)
    at_b = Atom(110, "CA", ' ', [0.0, -1.0, 3.0], 1.0, 10.0, " C", "  ", res)
    @test sqdistance(at_a, at_b) == 10.0
    @test isapprox(distance(at_a, at_b), sqrt(10))

    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B']), sqrt(6.852947))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'][50]), sqrt(530.645746))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'][50]["CA"]), sqrt(574.699125))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'], backboneselector), sqrt(11.752440))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'], standardselector), sqrt(11.252973))
    @test isapprox(distance(struc_1AKE['A'][50]["CA"], struc_1AKE['B'][50]["CA"]), sqrt(2607.154834))

    # Test bondangle
    at_a = Atom(100, "CA", ' ', [1.0, 0.0, 1.0], 1.0, 10.0, " C", "  ", res)
    at_b = Atom(100, "CA", ' ', [0.0, 0.0, 0.0], 1.0, 10.0, " C", "  ", res)
    at_c = Atom(100, "CA", ' ', [3.0, 2.0, 1.0], 1.0, 10.0, " C", "  ", res)
    @test isapprox(bondangle(at_a, at_b, at_c), 0.713724, atol=1e-5)
    vec_a = [2.0, 0.0, 0.0]
    vec_b = [2.0, 1.0, 1.0]
    @test isapprox(bondangle(vec_a, vec_b), 0.615480, atol=1e-5)

    # Test dihedral functions
    at_a = Atom(100, "CA", ' ', [-1.0, -1.0, 0.0], 1.0, 10.0, " C", "  ", res)
    at_b = Atom(100, "CA", ' ', [0.0, 0.0, 0.0], 1.0, 10.0, " C", "  ", res)
    at_c = Atom(100, "CA", ' ', [1.0, 0.0, 0.0], 1.0, 10.0, " C", "  ", res)
    at_d = Atom(100, "CA", ' ', [2.0, 1.0, -1.0], 1.0, 10.0, " C", "  ", res)
    @test isapprox(dihedralangle(at_a, at_b, at_c, at_d), 2.35619, atol=1e-5)
    vec_a = [1.0, 1.0, 0.0]
    vec_b = [1.0, 0.0, 0.0]
    vec_c = [1.0, -1.0, 1.0]
    @test isapprox(dihedralangle(vec_a, vec_b, vec_c), -0.785398, atol=1e-5)

    @test isapprox(omegaangle(struc_1AKE['A'][20], struc_1AKE['A'][19]), -3.09191, atol=1e-5)
    @test isapprox(phiangle(struc_1AKE['A'][7], struc_1AKE['A'][6]), 2.85115, atol=1e-5)
    @test isapprox(psiangle(struc_1AKE['A'][8], struc_1AKE['A'][9]), 2.83827, atol=1e-5)
    @test_throws ArgumentError omegaangle(struc_1AKE['A'][20], Residue("ALA", 19, ' ', false, Chain('A')))
    @test_throws ArgumentError phiangle(struc_1AKE['A'][7], Residue("ALA", 6, ' ', false, Chain('A')))
    @test_throws ArgumentError psiangle(struc_1AKE['A'][8], Residue("ALA", 9, ' ', false, Chain('A')))

    @test isapprox(omegaangle(struc_1AKE['A'], 20), omegaangle(struc_1AKE['A'][20], struc_1AKE['A'][19]), atol=1e-5)
    @test isapprox(phiangle(struc_1AKE['A'], 7), phiangle(struc_1AKE['A'][7], struc_1AKE['A'][6]), atol=1e-5)
    @test isapprox(psiangle(struc_1AKE['A'], 8), psiangle(struc_1AKE['A'][8], struc_1AKE['A'][9]), atol=1e-5)
    @test_throws ArgumentError omegaangle(struc_1AKE['A'], -1 )
    @test_throws ArgumentError phiangle(struc_1AKE['A'], -1 )
    @test_throws ArgumentError psiangle(struc_1AKE['A'], -1 )
    @test_throws ArgumentError omegaangle(struc_1AKE['A'], 1)
    @test_throws ArgumentError phiangle(struc_1AKE['A'], 1)
    @test_throws ArgumentError psiangle(struc_1AKE['A'], 214)

    phis, psis = ramachandranangles(struc_1AKE['A'])
    @test size(phis) == (456,)
    @test size(psis) == (456,)
    @test isapprox(phis[5], -1.76451, atol=1e-5)
    @test isapprox(psis[10], 0.442509, atol=1e-5)
    @test isnan(phis[1])
    @test isnan(psis[214])
    @test sum(isnan, phis) == 243
    @test sum(isnan, psis) == 243
    @test_throws ArgumentError ramachandranangles(struc_1AKE['A'][10]["CA"])

    phis = phiangles(struc_1AKE['A'], standardselector)
    psis = psiangles(struc_1AKE['A'], standardselector)
    omegas = omegaangles(struc_1AKE['A'], standardselector)
    @test length(phis) == countresidues(struc_1AKE['A'], standardselector)
    @test length(phis) == length(psis)
    @test length(phis) == length(omegas)
    @test isapprox(psis[10], psiangle(struc_1AKE['A'], 10), atol=1e-5)
    @test isapprox(phis[10], phiangle(struc_1AKE['A'], 10), atol=1e-5)
    @test isapprox(omegas[10], omegaangle(struc_1AKE['A'], 10), atol=1e-5)

    # Test that the entries in `chitables` are bonded
    sortt((a, b)) = a < b ? (a, b) : (b, a)
    rd = BioStructures.residuedata
    for ct in BioStructures.chitables
        for (rname, alist) in ct
            if rname == "HIS"
                rname = "HID"
            end
            sb = sortt.(rd[rname].bonds)
            for i = 1:3
                @test sortt((alist[i], alist[i+1])) ∈ sb
            end
        end
    end
    chis = chiangles(struc_1AKE['A'], standardselector)
    @test length(chis) == countresidues(struc_1AKE['A'], standardselector)
    @test chis[1] ≈ deg2rad.([-176.231, 172.056, 56.069]) rtol=1e-5   # MET
    @test chis[2] ≈ deg2rad.([-76.551, 171.696, 171.162, -175.969, 0.424]) rtol=1e-5   # ARG
    @test isempty(chis[7])       # GLY
    @test chiangle(struc_1AKE['A'][2], 3) == chis[2][3]
    @test_throws "GLY does not have any χ angles" chiangle(struc_1AKE['A'][7], 1)
    @test_throws "χ angle with index 4 does not exist for residue MET (max is 3)" chiangle(struc_1AKE['A'][1], 4)
    fakeres = Residue("UKN", 10, ' ', false, struc_1AKE['A'])
    @test_throws "no χ angles are defined for residues with name UKN" chiangle(fakeres, 1)
    @test_throws "no χ angles are defined for residues with name UKN" chiangles(fakeres)
    # Test symmetries: two ASPs with OD1 and OD2 swapped
    io = IOBuffer("""
        ATOM    534  N   ASP A   1      -8.068   7.150  -2.008  1.00 95.46           N
        ATOM    535  CA  ASP A   1      -8.464   6.572  -3.297  1.00 95.46           C
        ATOM    536  C   ASP A   1      -8.522   7.636  -4.409  1.00 95.46           C
        ATOM    537  CB  ASP A   1      -9.844   5.896  -3.142  1.00 95.46           C
        ATOM    538  O   ASP A   1      -8.015   7.431  -5.518  1.00 95.46           O
        ATOM    539  CG  ASP A   1      -9.818   4.542  -2.418  1.00 95.46           C
        ATOM    540  OD1 ASP A   1      -8.738   3.930  -2.377  1.00 95.46           O
        ATOM    541  OD2 ASP A   1     -10.899   4.073  -1.981  1.00 95.46           O
        """)
    r = only(only(only(read(io, PDBFormat))))
    χs1 = chiangles(r)
    io = IOBuffer("""
        ATOM    534  N   ASP A   1      -8.068   7.150  -2.008  1.00 95.46           N
        ATOM    535  CA  ASP A   1      -8.464   6.572  -3.297  1.00 95.46           C
        ATOM    536  C   ASP A   1      -8.522   7.636  -4.409  1.00 95.46           C
        ATOM    537  CB  ASP A   1      -9.844   5.896  -3.142  1.00 95.46           C
        ATOM    538  O   ASP A   1      -8.015   7.431  -5.518  1.00 95.46           O
        ATOM    539  CG  ASP A   1      -9.818   4.542  -2.418  1.00 95.46           C
        ATOM    540  OD1 ASP A   1     -10.899   4.073  -1.981  1.00 95.46           O
        ATOM    541  OD2 ASP A   1      -8.738   3.930  -2.377  1.00 95.46           O
        """)
    r = only(only(only(read(io, PDBFormat))))
    χs2 = chiangles(r)
    @test χs1[1] ≈ χs2[1]
    @test χs1[2] ≈ χs2[2] rtol=0.05  # the bonds are not exactly symmetric

    # Test ContactMap
    cas = collectatoms(struc_1AKE, calphaselector)[1:10]
    @test isa(ContactMap(cas, 10).data, BitArray{2})
    @test ContactMap(cas, 10).data == [
        true  true  true  false false false false false false false
        true  true  true  true  false false false false false false
        true  true  true  true  true  true  false false false false
        false true  true  true  true  true  false false false false
        false false true  true  true  true  true  false false false
        false false true  true  true  true  true  true  true  false
        false false false false true  true  true  true  true  true
        false false false false false true  true  true  true  true
        false false false false false true  true  true  true  true
        false false false false false false true  true  true  true
    ]
    @test ContactMap(struc_1AKE[1], 1.0).data == [
        true  false
        false true
    ]
    cmap = ContactMap(struc_1AKE['A'], 5.0)
    @test size(cmap) == (456, 456)
    @test size(cmap, 2) == 456
    @test cmap[196, 110]
    @test !cmap[15, 89]
    cmap[15, 89] = true
    @test cmap[15, 89]
    show(devnull, cmap)

    @test ContactMap(struc_1AKE['A'][10], struc_1AKE['A'][11], 4.0).data == [
        true  false false false false
        true  true  false false false
        true  true  true  false true
        true  true  true  false false
    ]

    # Test the plot recipe
    RecipesBase.apply_recipe(Dict{Symbol, Any}(), cmap)

    showcontactmap(devnull, cmap)

    # Test DistanceMap
    @test isa(DistanceMap(cas).data, Array{Float64, 2})
    @test isapprox(DistanceMap(cas).data, [
        0.0                3.8116564640586357 6.5343511537106735 10.30641673909997  12.934087289020436 16.277628113456824 19.66671792648687  23.24696586653837  24.371038508853083 24.24882335289694
        3.8116564640586357 0.0                3.8115826109373554 7.075483093047433  10.24737590800689  13.368176240609639 16.959531449895664 20.432823079545326 21.590401316325735 21.931261910797566
        6.5343511537106735 3.8115826109373554 0.0                3.844083635926775  6.56909392534465   9.783285439973628  13.273855091871388 16.83042189013692  18.152391219891662 18.44135962991883
        10.30641673909997  7.075483093047433  3.844083635926775  0.0                3.852368622029827  6.369990659333809  10.066898131996766 13.489078471118772 14.9800810411693   15.83387716259034
        12.934087289020436 10.24737590800689  6.56909392534465   3.852368622029827  0.0                3.8467504468057196 6.845814633774421  10.314890304797236 11.724292217443233 12.19677604943208
        16.277628113456824 13.368176240609639 9.783285439973628  6.369990659333809  3.8467504468057196 0.0                3.8515559193655755 7.381047825342957  9.702061069690293  11.058714075334438
        19.66671792648687  16.959531449895664 13.273855091871388 10.066898131996766 6.845814633774421  3.8515559193655755 0.0                3.8486537906130245 6.699431169883006  8.065864739753573
        23.24696586653837  20.432823079545326 16.83042189013692  13.489078471118772 10.314890304797236 7.381047825342957  3.8486537906130245 0.0                3.86929786912303   6.311746667920065
        24.371038508853083 21.590401316325735 18.152391219891662 14.9800810411693   11.724292217443233 9.702061069690293  6.699431169883006  3.86929786912303   0.0                3.8548294385095696
        24.24882335289694  21.931261910797566 18.44135962991883  15.83387716259034  12.19677604943208  11.058714075334438 8.065864739753573  6.311746667920065  3.8548294385095696 0.0
    ], atol=1e-5)
    @test isapprox(DistanceMap(struc_1AKE[1]).data, [
        0.0     2.61781
        2.61781 0.0
    ], atol=1e-5)
    dmap = DistanceMap(struc_1AKE['A'])
    @test size(dmap) == (456, 456)
    @test size(dmap, 2) == 456
    @test isapprox(dmap[196, 110], 3.01887, atol=1e-5)
    dmap[196, 110] = 10.0
    @test dmap[196, 110] == 10.0
    @test isapprox(dmap[begin + 1], 1.3405938982406227, atol=1e-5)
    @test isapprox(dmap[end   - 1], 18.716774962583695, atol=1e-5)
    @test isapprox(dmap[2, begin + 2], 1.329450262326499 , atol=1e-5)
    @test isapprox(dmap[2, end   - 2], 18.595504483611087, atol=1e-5)
    show(devnull, dmap)

    @test isapprox(DistanceMap(struc_1AKE['A'][10], struc_1AKE['A'][11]).data, [
        2.8803883418733673 4.295162278657233  5.211988200293629  6.427722458227334 5.03904008715946
        2.5066699024801813 3.850861981427013  4.557738913101538  5.762204526047301 5.0027798272560435
        1.3632875705440877 2.4469166311911783 3.1850982716393537 4.33450689236965  3.799490623754715
        2.2682515292621335 2.774401016435799  3.235280049702036  4.198412795331114 4.313454647959104
    ], atol=1e-5)

    # Test the plot recipe
    RecipesBase.apply_recipe(Dict{Symbol, Any}(), dmap)

    # Test atom graph constructor
    cbetas = collectatoms(struc_1AKE["A"], cbetaselector)
    mg = MetaGraph(cbetas, 8.0)
    @test nv(mg) == 214
    @test ne(mg) == 1027
    @test get_prop(mg, :contactdist) == 8.0
    @test mg[10, :element] == cbetas[10]

    mg = MetaGraph(struc_1AKE[1], 10.0)
    @test nv(mg) == 2
    @test ne(mg) == 1
    @test get_prop(mg, :contactdist) == 10.0
    @test mg[2, :element] == struc_1AKE["B"]
end

@testset "Secondary structure" begin
    res = Residue("ALA", 1, ' ', false, Chain('A'))
    @test sscode(res) == '-'
    sscode!(res, 'H')
    @test sscode(res) == 'H'

    res = Residue(
        "ALA",
        1,
        ' ',
        false,
        ["CA"],
        Dict(
            "CA" => Atom(
                1,
                "CA",
                ' ',
                [0.0, 0.0, 0.0],
                1.0,
                0.0,
                "  ",
                "  ",
                Residue("ALA", 1, ' ', false, Chain('A')),
            ),
        ),
        Chain('A'),
        'T'
    )
    @test sscode(res) == 'T'
    sscode!(res, 'P')
    @test sscode(res) == 'P'

    pdb_path  = downloadpdb("1BQ0", dir=temp_dir, format=PDBFormat)
    cif_path  = downloadpdb("1BQ0", dir=temp_dir, format=MMCIFFormat)

    struc = read(pdb_path, PDBFormat)
    for res in collectresidues(struc)[1:5]
        sscode!(res, 'T')
        @test sscode(res) == 'T'
    end

    struc = read(pdb_path, PDBFormat)
    rundssp!(struc)

    helix_atoms = collectatoms(struc, helixselector)
    sheet_atoms = collectatoms(struc, sheetselector)
    coil_atoms  = collectatoms(struc, coilselector)

    helixselector2(el) = sscodeselector(el, ['G', 'H', 'I', 'P'])
    helix_atoms2 = collectatoms(struc, helixselector2)
    helix_residues2 = collectresidues(struc, helixselector2)

    @test length(helix_atoms)  == 441
    @test length(sheet_atoms)  == 0
    @test length(coil_atoms)   == 803
    @test length(helix_atoms2) == 441

    helix_residues = collectresidues(struc, helixselector)
    sheet_residues = collectresidues(struc, sheetselector)
    coil_residues  = collectresidues(struc, coilselector)

    @test length(helix_residues)  == 25
    @test length(sheet_residues)  == 0
    @test length(coil_residues)   == 52
    @test length(helix_residues2) == 25

    @test sscode(helix_residues[1])  == 'H'
    @test sscode(helix_residues[11]) == 'H'
    @test sscode(helix_residues[21]) == 'H'

    struc2 = read(pdb_path, PDBFormat; run_dssp=true)
    @test sscode.(collectatoms(struc)) == sscode.(collectatoms(struc2))

    struc3 = read(cif_path, MMCIFFormat; run_dssp=true)
    @test sscode.(collectatoms(struc)) == sscode.(collectatoms(struc3))

    struc = read(pdb_path, PDBFormat)
    runstride!(struc)

    helix_atoms = collectatoms(struc, helixselector)
    sheet_atoms = collectatoms(struc, sheetselector)
    coil_atoms  = collectatoms(struc, coilselector)

    @test length(helix_atoms) == 595
    @test length(sheet_atoms) == 0
    @test length(coil_atoms)  == 649

    helix_residues = collectresidues(struc, helixselector)
    sheet_residues = collectresidues(struc, sheetselector)
    coil_residues  = collectresidues(struc, coilselector)

    @test length(helix_residues) == 34
    @test length(sheet_residues) == 0
    @test length(coil_residues)  == 43

    @test sscode(helix_residues[1])  == 'H'
    @test sscode(helix_residues[11]) == 'G'
    @test sscode(helix_residues[21]) == 'H'

    struc2 = read(pdb_path, PDBFormat; run_stride=true)
    @test sscode.(collectatoms(struc)) == sscode.(collectatoms(struc2))

    struc3 = read(cif_path, MMCIFFormat; run_stride=true)
    @test sscode.(collectatoms(struc)) == sscode.(collectatoms(struc3))

    isfile(temp_filename) && rm(temp_filename)
    rundssp(cif_path, temp_filename)
    @test countlines(temp_filename) > 50
    rm(temp_filename)
    runstride(pdb_path, temp_filename)
    @test countlines(temp_filename) > 100
    rm(temp_filename)
end

# Delete temporary file and temporary directory
rm(temp_filename, force=true)
rm(temp_dir, recursive=true, force=true)

end # TestBioStructures
