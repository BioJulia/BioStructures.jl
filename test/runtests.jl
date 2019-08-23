module TestBioStructures

using Test
using Format
using RecipesBase
using BioCore
using BioStructures
using BioStructures:
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
    requiresquote


# Skip download tests on Linux - required due to timeouts during CI
# Set to false locally to run these tests
skip_linux_download = true

fmtdir = BioCore.Testing.get_bio_fmt_specimens()

# Access files in BioFmtSpecimens to test against
testfilepath(path::AbstractString...) = joinpath(fmtdir, path...)

# All writing is done to one temporary file which is removed at the end
temp_filename, io = mktemp()
close(io)


@testset "PDB handling" begin
    if skip_linux_download && Sys.islinux()
        @info "Skipping download tests on Linux due to timeouts during CI; change skip_linux_download to run these tests"
    else
        @test length(pdbentrylist()) > 100000

        # This may be empty on a given date so we just check it has the correct type
        @test isa(pdbstatuslist("ftp://ftp.wwpdb.org/pub/pdb/data/status/latest/added.pdb"), Vector{String})
        # Invalid URL
        # The error type changes from ErrorException to ProcessFailedException in Julia v1.2
        #   therefore we check for the more general Exception type
        # This also applies to two examples below
        @test_throws Exception pdbstatuslist("ftp://ftp.wwpdb.org/pub/pdb/data/status/latest/dummy.pdb")

        addedlist, modifiedlist, obsoletelist = pdbrecentchanges()

        @test length(pdbobsoletelist()) > 3600

        pdb_dir = joinpath(tempdir(), "PDB")
        # Invalid PDB ID format
        @test_throws ArgumentError downloadpdb("1a df")
        # Valid PDB ID format but PDB does not exist
        @test_throws Exception downloadpdb("no1e", pdb_dir=pdb_dir)
        # Invalid PDB file_format
        @test_throws ArgumentError downloadpdb("1alw", pdb_dir=pdb_dir, file_format=String)
        # Biological assembly not available in PDBXML and MMTF
        @test_throws ArgumentError downloadpdb("1alw", pdb_dir=pdb_dir, file_format=PDBXML, ba_number=1)
        # Invalid ba_number for this PDB entry
        @test_throws Exception downloadpdb("1alw", pdb_dir=pdb_dir, file_format=MMCIF, ba_number=10)
        # Tests if downloadpdb returns the path to the downloaded file
        @test isfile(downloadpdb("1crn", pdb_dir=pdb_dir))

        # PDB format
        downloadpdb("1alw", pdb_dir=pdb_dir, file_format=PDB)
        pdbpath = joinpath(pdb_dir, "1ALW$(pdbextension[PDB])")
        @test isfile(pdbpath) && filesize(pdbpath) > 0
        # PDBXML format
        downloadpdb("1alw", pdb_dir=pdb_dir, file_format=PDBXML)
        pdbpath = joinpath(pdb_dir, "1ALW$(pdbextension[PDBXML])")
        @test isfile(pdbpath) && filesize(pdbpath) > 0
        # MMCIF format
        downloadpdb("1alw", pdb_dir=pdb_dir, file_format=MMCIF)
        pdbpath = joinpath(pdb_dir, "1ALW$(pdbextension[MMCIF])")
        @test isfile(pdbpath) && filesize(pdbpath) > 0
        # MMTF format
        downloadpdb("1alw", pdb_dir=pdb_dir, file_format=MMTF)
        pdbpath = joinpath(pdb_dir, "1ALW$(pdbextension[MMTF])")
        @test isfile(pdbpath) && filesize(pdbpath) > 0
        # Obsolete PDB
        downloadpdb("116l", pdb_dir=pdb_dir, file_format=PDB, obsolete=true)
        pdbpath = joinpath(pdb_dir, "obsolete", "116L$(pdbextension[PDB])")
        @test isfile(pdbpath) && filesize(pdbpath) > 0
        # Biological Assembly - PDB format
        downloadpdb("1alw", pdb_dir=pdb_dir, file_format=PDB, ba_number=1)
        pdbpath = joinpath(pdb_dir, "1ALW_ba1$(pdbextension[PDB])")
        @test isfile(pdbpath) && filesize(pdbpath) > 0
        # Biological Assembly - MMCIF format
        downloadpdb("5a9z", pdb_dir=pdb_dir, file_format=MMCIF, ba_number=1)
        pdbpath = joinpath(pdb_dir, "5A9Z_ba1$(pdbextension[MMCIF])")
        @test isfile(pdbpath) && filesize(pdbpath) > 0
        # Download multiple PDB files
        pdbidlist = ["1ent", "1en2"]
        downloadpdb(pdbidlist, pdb_dir=pdb_dir, file_format=PDB)
        for pdbid in pdbidlist
            pdbpath = joinpath(pdb_dir, "$(uppercase(pdbid))$(pdbextension[PDB])")
            @test isfile(pdbpath) && filesize(pdbpath) > 0
        end

        # Test function as first argument to downloadpdb
        @test downloadpdb(countlines, "1alw") == 3600
        @test downloadpdb("1alw") do fp
            s = read(fp, PDB)
            return countresidues(s, standardselector)
        end == 346
        @test downloadpdb(["169l", "2lzm"]) do fps
            n_chains = Int[]
            for fp in fps
                s = read(fp, PDB)
                push!(n_chains, countchains(s))
            end
            return n_chains
        end == [5, 1]

        # Test retrievepdb and readpdb
        struc = retrievepdb("1AKE", pdb_dir=pdb_dir, structure_name="New name")
        @test structurename(struc) == "New name"
        @test countatoms(struc) == 3804

        struc = retrievepdb("1AKE", pdb_dir=pdb_dir, obsolete=true, read_het_atoms=false)
        @test countatoms(struc) == 3312
        @test serial(collectatoms(struc)[2000]) == 2006
        @test sum(ishetero, collectatoms(struc)) == 0

        struc = retrievepdb("1AKE", pdb_dir=pdb_dir, ba_number=1, read_het_atoms=false, read_std_atoms=false)
        @test countatoms(struc) == 0
        @test countresidues(struc) == 0
        @test countchains(struc) == 0
        @test countmodels(struc) == 0

        struc = readpdb("1AKE", pdb_dir=pdb_dir, read_std_atoms=false)
        @test countatoms(struc) == 492
        @test serial(collectatoms(struc)[400]) == 3726
        @test sum(ishetero, collectatoms(struc)) == 492

        struc = readpdb("1AKE", pdb_dir=pdb_dir, ba_number=1, remove_disorder=true)
        @test countatoms(struc) == 1954
        @test sum(isdisorderedatom, collectatoms(struc)) == 0
        @test tempfactor(struc['A'][167]["NE"]) == 23.32
    end
end


@testset "Model" begin
    # Test constructors and indexing
    struc = ProteinStructure("Test structure")
    struc[1] = Model(1, struc)
    mod = struc[1]
    @test isa(mod, Model)
    struc[3] = Model(3, struc)
    struc[1]['A'] = Chain('A', mod)
    ch = struc[1]['A']
    @test isa(ch, Chain)
    struc['B'] = Chain('B', mod)
    struc['A'][10] = Residue("ALA", 10, ' ', false, ch)
    res = struc['A'][10]
    @test isa(res, Residue)
    struc['A']["H_20A"] = DisorderedResidue(Dict(
        " VA"=> Residue(" VA", 20, 'A', true, ch),
        "ILE"=> Residue("ILE", 20, 'A', true, ch)
    ), " VA")
    dis_res = struc['A']["H_20A"]
    @test isa(dis_res, DisorderedResidue)
    struc['A'][10][" CA "] = Atom(
        100, " CA ", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res)
    at = struc['A'][10]["CA"]
    @test isa(at, Atom)
    struc['A'][10][" CB "] = DisorderedAtom(Dict(
        'A'=> Atom(200, " CB ", 'A', [10.0, 20.0, 30.0], 0.6, 20.0, " C", "1+", res),
        'B'=> Atom(201, " CB ", 'B', [11.0, 21.0, 31.0], 0.4, 30.0, " C", "1+", res)
    ), 'A')
    dis_at = struc['A'][10]["CB"]
    @test isa(dis_at, DisorderedAtom)
    struc['A']["H_20A"][" CG "] = Atom(
        300, " CG ", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", defaultresidue(dis_res))
    disorderedres(dis_res, "ILE")[" O  "] = Atom(
        400, " O  ", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " O", "  ", disorderedres(dis_res, "ILE"))
    fixlists!(struc)


    # Test alternate constructors
    ProteinStructure("struc", Dict(1=> Model()))
    ProteinStructure()

    Model(1, Dict("A"=> Chain('A')), ProteinStructure())
    Model(1, ProteinStructure())
    Model(1)
    Model()

    Chain("A", ["1"], Dict("1"=> Residue("ALA", 1, ' ', false, Chain('A'))), Model())
    Chain("A", Model())
    Chain('A', Model())
    Chain("A")
    Chain('A')

    Residue("ALA", 1, ' ', false, ["CA"], Dict("CA"=>
        Atom(1, "CA", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, "  ", "  ",
        Residue("ALA", 1, ' ', false, Chain('A')))), Chain('A'))
    Residue("ALA", 1, ' ', false, Chain('A'))


    # Test show
    show(devnull, at)
    show(devnull, dis_at)
    show(devnull, res)
    show(devnull, dis_res)
    show(devnull, ch)
    show(devnull, mod)
    show(devnull, struc)
    show(devnull, Model())
    show(devnull, ProteinStructure())


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

    @test resids(ch) == ["10", "H_20A"]

    @test isa(residues(ch), Dict{String, AbstractResidue})
    @test length(residues(ch)) == 2
    @test serial(residues(ch)["10"]["CA"]) == 100

    @test model(at) == mod
    @test model(dis_at) == mod
    @test model(res) == mod
    @test model(dis_res) == mod
    @test model(ch) == mod
    @test model(mod) == mod

    @test modelnumber(at) == 1
    @test modelnumber(dis_at) == 1
    @test modelnumber(res) == 1
    @test modelnumber(dis_res) == 1
    @test modelnumber(ch) == 1
    @test modelnumber(mod) == 1

    @test chainids(mod) == ["A", "B"]
    @test chainids(struc) == ["A", "B"]
    @test chainids(ProteinStructure()) == String[]
    @test chainids(Model()) == String[]

    @test isa(chains(mod), Dict{String, Chain})
    @test length(chains(mod)) == 2
    @test resname(chains(mod)["A"]["H_20A"]) == "VA"
    @test isa(chains(Model()), Dict{String, Chain})
    @test isa(chains(struc), Dict{String, Chain})
    @test length(chains(struc)) == 2
    @test isa(chains(ProteinStructure()), Dict{String, Chain})

    @test structure(at) == struc
    @test structure(dis_at) == struc
    @test structure(res) == struc
    @test structure(dis_res) == struc
    @test structure(ch) == struc
    @test structure(mod) == struc
    @test structure(struc) == struc

    @test structurename(at) == "Test structure"
    @test structurename(dis_at) == "Test structure"
    @test structurename(res) == "Test structure"
    @test structurename(dis_res) == "Test structure"
    @test structurename(ch) == "Test structure"
    @test structurename(mod) == "Test structure"
    @test structurename(struc) == "Test structure"

    @test modelnumbers(struc) == [1, 3]

    @test isa(models(struc), Dict{Int, Model})
    @test length(models(struc)) == 2
    @test resids(models(struc)[1]['A']) == ["10", "H_20A"]

    @test isa(defaultmodel(struc), Model)
    @test modelnumber(defaultmodel(struc)) == 1


    # Test iteration over elements
    at_col = collect(at)
    @test isa(at_col, Vector{Atom})
    @test length(at_col) == 1
    @test serial(at_col[1]) == 100

    dis_at_col = collect(dis_at)
    @test isa(dis_at_col, Vector{Atom})
    @test length(dis_at_col) == 2
    @test serial(dis_at_col[2]) == 201

    res_col = collect(res)
    @test isa(res_col, Vector{AbstractAtom})
    @test length(res_col) == 2
    @test serial(res_col[2]) == 200

    dis_res_col = collect(dis_res)
    @test isa(dis_res_col, Vector{AbstractAtom})
    @test length(dis_res_col) == 1
    @test serial(dis_res_col[1]) == 300

    ch_col = collect(ch)
    @test isa(ch_col, Vector{AbstractResidue})
    @test length(ch_col) == 2
    @test resname(ch_col[2]) == "VA"

    mod_col = collect(mod)
    @test isa(mod_col, Vector{Chain})
    @test length(mod_col) == 2
    @test chainid(mod_col[2]) == "B"

    struc_col = collect(struc)
    @test isa(struc_col, Vector{Model})
    @test length(struc_col) == 2
    @test modelnumber(struc_col[2]) == 3


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

    @test isa(mod["A"], Chain)
    @test isa(mod['A'], Chain)
    @test resname(mod['A']["H_20A"]) == "VA"
    @test_throws KeyError mod['C']

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

    # Order when sorting a residue list is chain ID, then stdres/hetres,
    # then residue number, then insertion code
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
    @test map(res -> resid(res, full=true), res_ord) == [
        "201A:A", "203:A", "H_200:A", "201B:A", "202:A", "300:B", "H_201:A",
        "201:A", "H_201A:A", "100:B", "H_203:A", "200:A"]
    @test map(res -> resid(res, full=true), sort(res_ord)) == [
        "200:A", "201:A", "201A:A", "201B:A", "202:A", "203:A", "H_200:A",
        "H_201:A", "H_201A:A", "H_203:A", "100:B", "300:B"]
    sort!(res_ord)
    @test map(res -> resid(res, full=true), res_ord) == [
        "200:A", "201:A", "201A:A", "201B:A", "202:A", "203:A", "H_200:A",
        "H_201:A", "H_201A:A", "H_203:A", "100:B", "300:B"]

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
    mod_ord = Model(1, Dict(
        "AB" => Chain("AB"),
        "A" => Chain("A"),
        " " => Chain(" "),
        "1" => Chain("1"),
        "AAA" => Chain("AAA"),
        "BC" => Chain("BC"),
        "a" => Chain("a"),
        "X" => Chain("X"),
        "AC" => Chain("AC"),
    ), ProteinStructure())
    @test chainids(mod_ord) == ["1", "A", "X", "a", "AB", "AC", "BC", "AAA", " "]


    # Test sequence extraction
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDB)
    seq = AminoAcidSequence(struc['B'])
    @test seq == AminoAcidSequence(
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
    seq = AminoAcidSequence(struc['B'], standardselector)
    @test seq == AminoAcidSequence(
        "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNG" *
        "FLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRIVGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQ" *
        "EETVRKRLVEYHQMTAPLIGYYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG"
    )
    seq = AminoAcidSequence(AbstractResidue[
        Residue("VAL", 20, 'B', true, Chain('B')),
        Residue("ALA", 10, 'A', false, Chain('A')),
    ])
    @test seq == AminoAcidSequence("VA")
end


@testset "Parsing" begin
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
    @test_throws PDBParseError    parseserial("ATOM         C   GLY A  80      29.876  54.131  35.806  1.00 40.97           C1+")
    @test_throws PDBParseError  parseatomname(line_short)
    @test_throws PDBParseError    parsealtloc(line_short)
    @test_throws PDBParseError   parseresname(line_short)
    @test_throws PDBParseError   parsechainid(line_short)
    @test_throws PDBParseError parseresnumber("ATOM    591  C   GLY A          29.876  54.131  35.806  1.00 40.97           C1+")
    @test_throws PDBParseError   parseinscode(line_short)
    @test_throws PDBParseError    parsecoordx("ATOM    591  C   GLY A  80      xxxxxx  54.131  35.806  1.00 40.97           C1+")
    @test_throws PDBParseError    parsecoordy("ATOM    591  C   GLY A  80      29.876  xxxxxx  35.806  1.00 40.97           C1+")
    @test_throws PDBParseError    parsecoordz("ATOM    591  C   GLY A  80      29.876  54.131  xxxxxx  1.00 40.97           C1+")
    @test parseoccupancy(line_short) == 1.0
    @test parsetempfac(line_short) == 0.0
    @test parseelement(line_short) == "  "
    @test parsecharge(line_short) == "  "


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
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDB)
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
    mods = collect(struc)
    @test modelnumber(mods[1]) == 1
    chs = collect(mods[1])
    @test chainid.(chs) == ["A", "B"]
    res = collect(chs[1])
    @test length(res) == 456
    @test resid(res[20]) == "20"
    ats = collect(res[20])
    @test length(ats) == 8
    @test atomname.(ats) == ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"]


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
    @test map(res -> resid(res, full=true), res_min) == map(res -> resid(res, full=true), res)
    applyselectors!(res_min)
    @test length(res_min) == length(res)
    @test map(res -> resid(res, full=true), res_min) == map(res -> resid(res, full=true), res)
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
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDB, structure_name="New name")
    @test structurename(struc) == "New name"
    @test countatoms(struc) == 3804

    struc = read(testfilepath("PDB", "1AKE.pdb"), PDB, read_het_atoms=false)
    @test countatoms(struc) == 3312
    @test serial(collectatoms(struc)[2000]) == 2006
    @test sum(ishetero, collectatoms(struc)) == 0

    struc = read(testfilepath("PDB", "1AKE.pdb"), PDB, read_std_atoms=false)
    @test countatoms(struc) == 492
    @test serial(collectatoms(struc)[400]) == 3726
    @test sum(ishetero, collectatoms(struc)) == 492

    struc = read(testfilepath("PDB", "1AKE.pdb"), PDB, read_het_atoms=false, read_std_atoms=false)
    @test countatoms(struc) == 0
    @test countresidues(struc) == 0
    @test countchains(struc) == 0
    @test countmodels(struc) == 0

    struc = read(testfilepath("PDB", "1AKE.pdb"), PDB, remove_disorder=true)
    @test countatoms(struc) == 3804
    @test sum(isdisorderedatom, collectatoms(struc)) == 0
    @test tempfactor(struc['A'][167]["NE"]) == 23.32

    # Test parsing from stream
    open(testfilepath("PDB", "1AKE.pdb")) do file
        struc = read(file, PDB)
        @test countatoms(struc) == 3804
        @test countresidues(struc) == 808
    end


    # Test parsing 1EN2 (disordered residue)
    struc = read(testfilepath("PDB", "1EN2.pdb"), PDB)
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
    @test resnumber(res[1]) == 10
    @test countresidues(DisorderedResidue[struc['A'][16], struc['A'][10]]) == 2


    # Test parsing 1SSU (multiple models)
    struc = read(testfilepath("PDB", "1SSU.pdb"), PDB)
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
    @test z(ats[20]) == -16.522
    @test z(ats[1000]) == -0.394
    @test countatoms(Model[struc[5], struc[10]]) == 1512
    res = collectresidues(Model[struc[5], struc[10]])
    @test length(res) == 102
    @test y(res[10]["O"]) == -6.421
    @test y(res[100]["O"]) == -15.66
    @test countresidues(Model[struc[5], struc[10]]) == 102
    chs = collectchains(Model[struc[5], struc[10]])
    @test length(chs) == 2
    @test chainid.(chs) == ["A", "A"]
    @test z(chs[2][5]["CA"]) == -5.667
    @test countchains(Model[struc[5], struc[10]]) == 2
    mods = collectmodels(Model[struc[10], struc[5]])
    @test length(mods) == 2
    @test modelnumber.(mods) == [5, 10]
    @test z(mods[2]['A'][5]["CA"]) == -5.667
    @test countmodels(Model[struc[10], struc[5]]) == 2


    # Test collectatoms
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDB)
    ats = collectatoms(struc)
    @test length(ats) == 3804
    @test isa(ats, Vector{AbstractAtom})
    @test isa(ats[70], Atom)
    @test isa(ats[1290], DisorderedAtom)
    @test serial(ats[1660]) == 1666
    ats = collectatoms(struc, heteroselector)
    @test length(ats) == 492
    @test serial(ats[80]) == 3406
    ats = collectatoms(struc, disorderselector)
    @test length(ats) == 12
    @test serial(ats[10]) == 3338
    ats = collectatoms(struc, standardselector, disorderselector)
    @test length(ats) == 5
    @test serial(ats[4]) == 1294
    ats = collectatoms(struc[1])
    @test length(ats) == 3804
    @test serial(ats[1660]) == 1666
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
    @test serial(ats[5]) == 5
    ats = collectatoms(Residue[struc['A'][51], struc['A'][50]])
    @test length(ats) == 17
    @test serial(ats[5]) == 359
    ats = collectatoms(AbstractResidue[struc['A'][51], struc['A'][50]])
    @test length(ats) == 17
    @test serial(ats[5]) == 359
    ats = collectatoms(Atom[struc['A'][51]["CA"], struc['A'][50]["CA"]])
    @test length(ats) == 2
    @test isa(ats, Vector{Atom})
    @test serial(ats[2]) == 365
    ats = collectatoms(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]])
    @test length(ats) == 2
    @test isa(ats, Vector{DisorderedAtom})
    @test serial(ats[2]) == 1292
    ats = collectatoms(AbstractAtom[struc['A'][50]["CA"], struc['A'][167]["CZ"]])
    @test length(ats) == 2
    @test isa(ats, Vector{AbstractAtom})
    @test serial(ats[2]) == 1292


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
    @test countatoms(struc['A'], heteroselector) == 298
    @test countatoms(struc['A'], standardselector, disorderselector) == 5

    @test countatoms(ProteinStructure()) == 0
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
    @test resid(res[5], full=true) == "5:A"
    res = collectresidues(Residue[struc['A'][51], struc['A'][50]])
    @test length(res) == 2
    @test isa(res, Vector{Residue})
    @test resnumber(res[1]) == 50
    res = collectresidues(AbstractResidue[struc['A'][51], struc['A'][50]])
    @test length(res) == 2
    @test isa(res, Vector{AbstractResidue})
    @test resnumber(res[1]) == 50
    res = collectresidues(Atom[struc['A'][51]["CA"], struc['A'][50]["CA"]])
    @test length(res) == 2
    @test isa(res, Vector{AbstractResidue})
    @test resnumber(res[2]) == 51
    res = collectresidues(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]])
    @test length(res) == 1
    @test isa(res, Vector{AbstractResidue})
    @test atomnames(res[1]) == ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"]
    res = collectresidues(AbstractAtom[struc['A'][50]["CA"], struc['A'][167]["CZ"]])
    @test length(res) == 2
    @test isa(res, Vector{AbstractResidue})
    @test atomnames(res[1]) == ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"]


    # Test countresidues
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

    @test countresidues(ProteinStructure()) == 0
    @test countresidues(Model()) == 0
    @test countresidues(Chain('X')) == 0
    @test countresidues(Residue("ALA", 100, ' ', false, Chain('A'))) == 1


    # Test collectchains
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
    @test chainid(chs[2]) == "B"
    chs = collectchains(Residue[struc['A'][51], struc['B'][50]])
    @test length(chs) == 2
    @test chainid(chs[2]) == "B"
    chs = collectchains(AbstractResidue[struc['A'][51], struc['B'][50]])
    @test length(chs) == 2
    @test chainid(chs[2]) == "B"
    chs = collectchains(Atom[struc['B'][51]["CA"], struc['A'][50]["CA"]])
    @test length(chs) == 2
    @test chainid(chs[2]) == "B"
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

    @test countchains(ProteinStructure()) == 0
    @test countchains(Model()) == 0
    @test countchains(Chain('X')) == 1
    @test countchains(Residue("ALA", 100, ' ', false, Chain('A'))) == 1


    # Test collectmodels
    struc_1SSU = read(testfilepath("PDB", "1SSU.pdb"), PDB)
    mods = collectmodels(struc_1SSU)
    @test length(mods) == 20
    @test isa(mods, Vector{Model})
    @test modelnumber.(mods) == collect(1:20)
    mods = collectmodels(struc_1SSU, mod -> modelnumber(mod) < 4)
    @test length(mods) == 3
    @test modelnumber(mods[2]) == 2
    mods = collectmodels(struc_1SSU[10])
    @test length(mods) == 1
    @test modelnumber(mods[1]) == 10
    mods = collectmodels(struc_1SSU['A'])
    @test length(mods) == 1
    @test modelnumber(mods[1]) == 1
    mods = collectmodels(struc_1SSU[7]['A'][50])
    @test length(mods) == 1
    @test modelnumber(mods[1]) == 7
    mods = collectmodels(struc_1SSU['A'][50]["CA"])
    @test length(mods) == 1
    @test modelnumber(mods[1]) == 1
    mods = collectmodels(Chain[struc['B'], struc['A']])
    @test length(mods) == 1
    @test modelnumber(mods[1]) == 1
    mods = collectmodels(Residue[struc['A'][51], struc['B'][50]])
    @test length(mods) == 1
    @test modelnumber(mods[1]) == 1
    mods = collectmodels(AbstractResidue[struc_1SSU[12]['A'][51], struc_1SSU['A'][50]])
    @test length(mods) == 2
    @test modelnumber(mods[2]) == 12
    mods = collectmodels(Atom[struc_1SSU['A'][51]["CA"], struc_1SSU['A'][50]["CA"]])
    @test length(mods) == 1
    @test modelnumber(mods[1]) == 1
    mods = collectmodels(DisorderedAtom[struc['A'][167]["CZ"], struc['A'][167]["CD"]])
    @test length(mods) == 1
    @test modelnumber(mods[1]) == 1
    mods = collectmodels(AbstractAtom[struc['A'][50]["CA"], struc['A'][167]["CZ"]])
    @test length(mods) == 1
    @test modelnumber(mods[1]) == 1


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

    @test countmodels(struc_1SSU, mod -> modelnumber(mod) < 4) == 3

    @test countmodels(ProteinStructure()) == 0
    @test countmodels(Model()) == 1
    @test countmodels(Chain('X')) == 1
    @test countmodels(Residue("ALA", 100, ' ', false, Chain('A'))) == 1


    # Test parser error handling
    error = PDBParseError("message", 10, "line")
    showerror(devnull, error)
    # Missing coordinate (blank string)
    @test_throws PDBParseError read(testfilepath("PDB", "1AKE_err_a.pdb"), PDB)
    # Missing chain ID (line ends early)
    @test_throws PDBParseError read(testfilepath("PDB", "1AKE_err_b.pdb"), PDB)
    # Bad MODEL record
    @test_throws PDBParseError read(testfilepath("PDB", "1SSU_err.pdb"), PDB)
    # Truncated MODEL record from ASTRAL/SCOP95
    struc = read(testfilepath("PDB", "d9pcya_.ent"), PDB)
    @test isa(struc, ProteinStructure)
    @test countmodels(struc) == 16
    # Duplicate atom names in same residue
    @test_throws ErrorException read(testfilepath("PDB", "1AKE_err_c.pdb"), PDB)
    # Non-existent file
    @test_throws SystemError read(testfilepath("PDB", "non_existent_file.pdb"), PDB)


    # Test parsing empty file
    struc = read(IOBuffer(""), PDB)
    @test isa(struc, ProteinStructure)
    @test countmodels(struc) == 0
end


@testset "Writing" begin
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

    line_a = "ATOM    669  CA  ILE A  90      31.743  33.110  31.221  1.00 25.76           C  "
    line_b = "HETATM 3474  O  B XX A 334A      8.802  62.000   8.672  1.00 39.15           O1-"
    at_rec = AtomRecord(line_a)
    @test pdbline(at_rec) == "ATOM    669  CA  ILE A  90      31.743  33.110  31.221  1.00 25.76           C  "
    at_rec = AtomRecord(line_b)
    @test pdbline(at_rec) == "HETATM 3474  O  B XX A 334A      8.802  62.000   8.672  1.00 39.15           O1-"


    # Test writepdb
    struc = read(testfilepath("PDB", "1SSU.pdb"), PDB)
    writepdb(temp_filename, struc)
    @test countlines(temp_filename) == 15160
    struc_written = read(temp_filename, PDB)
    @test isa(struc_written, ProteinStructure)
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
    struc_written = read(temp_filename, PDB)
    @test modelnumbers(struc_written) == collect(1:20)
    @test countatoms(struc_written) == 756
    @test z(struc_written[4]['A']["30"]["OG"]) == -2.177
    @test atomnames(struc_written[15]['A']["39"]) == [
        "N", "CA", "C", "O", "CB", "SG", "H", "HA", "HB2", "HB3"]


    # Test selectors
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDB)
    writepdb(temp_filename, struc, heteroselector)
    @test countlines(temp_filename) == 499
    struc_written = read(temp_filename, PDB)
    @test modelnumbers(struc_written) == [1]
    @test countatoms(struc_written) == 492
    @test chainids(struc_written) == ["A", "B"]
    @test tempfactor(struc_written['B']["H_705"]["O"]) == 64.17
    writepdb(temp_filename, struc, standardselector, disorderselector)
    @test countlines(temp_filename) == 10
    struc_written = read(temp_filename, PDB)
    @test countatoms(struc_written) == 5
    @test sum(isdisorderedatom, collectatoms(struc_written)) == 5
    @test defaultaltlocid(struc_written['A'][167]["NH1"]) == 'A'


    # Test writing different element types
    writepdb(temp_filename, struc[1])
    @test countlines(temp_filename) == 3816
    struc_written = read(temp_filename, PDB)
    @test modelnumbers(struc_written) == [1]
    @test countatoms(struc_written) == 3804
    writepdb(temp_filename, struc['A'])
    @test countlines(temp_filename) == 1966
    struc_written = read(temp_filename, PDB)
    @test chainids(struc_written) == ["A"]
    writepdb(temp_filename, struc['A'][50])
    @test countlines(temp_filename) == 9
    struc_written = read(temp_filename, PDB)
    @test chainids(struc_written) == ["A"]
    @test countresidues(struc_written) == 1
    @test atomnames(struc_written['A'][50]) == [
        "N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"]
    writepdb(temp_filename, struc['A'][50]["CA"])
    @test countlines(temp_filename) == 1
    struc_written = read(temp_filename, PDB)
    @test countatoms(struc_written) == 1
    @test !isdisorderedatom(collectatoms(struc_written)[1])
    @test serial(collectatoms(struc_written)[1]) == 356
    writepdb(temp_filename, struc['A'][167]["CZ"])
    @test countlines(temp_filename) == 2
    struc_written = read(temp_filename, PDB)
    @test countatoms(struc_written) == 1
    @test isdisorderedatom(collectatoms(struc_written)[1])
    @test length(collectatoms(struc_written)[1]) == 2
    @test tempfactor(collectatoms(struc_written)[1]) == 16.77
    writepdb(temp_filename, Chain[struc['A'], struc['B']])
    @test countlines(temp_filename) == 3816
    struc_written = read(temp_filename, PDB)
    @test chainids(struc_written) == ["A", "B"]
    @test countatoms(struc_written['A']) == 1954
    @test countatoms(struc_written['B']) == 1850
    @test altlocids(struc_written['A']["H_215"]["O1G"]) == ['A', 'B']
    writepdb(temp_filename, AbstractResidue[struc['A'][51], struc['A'][50]])
    @test countlines(temp_filename) == 17
    struc_written = read(temp_filename, PDB)
    @test countresidues(struc_written) == 2
    @test resid.(collectresidues(struc_written)) == ["50", "51"]
    @test countatoms(struc_written) == 17
    writepdb(temp_filename, AbstractAtom[struc['A'][51]["CA"], struc['A'][50]["CA"]])
    @test countlines(temp_filename) == 2
    struc_written = read(temp_filename, PDB)
    @test countatoms(struc_written) == 2
    @test !ishetero(struc_written['A'][51]["CA"])


    # Test multiple model writing
    struc = read(testfilepath("PDB", "1SSU.pdb"), PDB)
    writepdb(temp_filename, Model[struc[10], struc[5]])
    @test countlines(temp_filename) == 1516
    struc_written = read(temp_filename, PDB)
    @test modelnumbers(struc_written) == [5, 10]
    @test modelnumber(defaultmodel(struc_written)) == 5
    @test countatoms(struc_written[5]) == 756
    @test countatoms(struc_written[10]) == 756
    @test_throws KeyError struc_written[1]


    # Test disordered residue writing
    struc = read(testfilepath("PDB", "1EN2.pdb"), PDB)
    writepdb(temp_filename, struc)
    @test countlines(temp_filename) == 819
    struc_written = read(temp_filename, PDB)
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
    struc_written = read(temp_filename, PDB)
    @test countresidues(struc_written) == 2
    @test isa(struc_written['A'][10], DisorderedResidue)
    @test isa(struc_written['A'][16], DisorderedResidue)
    @test defaultresname(struc_written['A'][10]) == "SER"
    @test isa(disorderedres(struc_written['A'][10], "GLY")["O"], Atom)
    @test altlocid(disorderedres(struc_written['A'][10], "GLY")["O"]) == 'B'
    @test countatoms(struc_written['A'][10]) == 6
    @test countatoms(struc_written['A'][16]) == 11

    @test_throws ArgumentError writepdb(temp_filename, Atom(
        1, "11H11", ' ', [0.0, 0.0, 0.0], 1.0, 0.0, " H", "  ", res))

    checkchainerror(Chain("A"))
    @test_throws ArgumentError checkchainerror(Chain("AA"))
    @test_throws ArgumentError writepdb(temp_filename,
        Residue("ALA", 10, ' ', false, [], Dict(), Chain("AA")))
end


@testset "mmCIF" begin
    # Test mmCIF dictionary
    dic = MMCIFDict()
    dic = MMCIFDict(Dict())
    dic = MMCIFDict(testfilepath("mmCIF", "1AKE.cif"))
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
    show(devnull, dic)

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
    struc = read(testfilepath("mmCIF", "1AKE.cif"), MMCIF)
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
    mods = collect(struc)
    @test modelnumber(mods[1]) == 1
    chs = collect(mods[1])
    @test chainid.(chs) == ["A", "B"]
    res = collect(chs[1])
    @test length(res) == 456
    @test resid(res[20]) == "20"
    ats = collect(res[20])
    @test length(ats) == 8
    @test atomname.(ats) == ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"]


    # Test parsing options
    struc = read(testfilepath("mmCIF", "1AKE.cif"), MMCIF, structure_name="New name")
    @test structurename(struc) == "New name"
    @test countatoms(struc) == 3804

    struc = read(testfilepath("mmCIF", "1AKE.cif"), MMCIF, read_het_atoms=false)
    @test countatoms(struc) == 3312
    # Different to the PDB file due to the lack of TER label serial
    @test serial(collectatoms(struc)[2000]) == 2005
    @test sum(ishetero, collectatoms(struc)) == 0

    struc = read(testfilepath("mmCIF", "1AKE.cif"), MMCIF, read_std_atoms=false)
    @test countatoms(struc) == 492
    @test serial(collectatoms(struc)[400]) == 3724
    @test sum(ishetero, collectatoms(struc)) == 492

    struc = read(testfilepath("mmCIF", "1AKE.cif"), MMCIF, read_het_atoms=false,
                read_std_atoms=false)
    @test countatoms(struc) == 0
    @test countresidues(struc) == 0
    @test countchains(struc) == 0
    @test countmodels(struc) == 0

    struc = read(testfilepath("mmCIF", "1AKE.cif"), MMCIF, remove_disorder=true)
    @test countatoms(struc) == 3804
    @test sum(isdisorderedatom, collectatoms(struc)) == 0
    @test tempfactor(struc['A'][167]["NE"]) == 23.32


    # Test parsing from stream
    open(testfilepath("mmCIF", "1AKE.cif")) do file
        struc = read(file, MMCIF)
        @test countatoms(struc) == 3804
        @test countresidues(struc) == 808
    end


    # Test parsing 1EN2
    struc = read(testfilepath("mmCIF", "1EN2.cif"), MMCIF)
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
    @test resnumber(res[1]) == 10
    @test countresidues(DisorderedResidue[struc['A'][16], struc['A'][10]]) == 2


    # Test parsing 1SSU
    struc = read(testfilepath("mmCIF", "1SSU.cif"), MMCIF)
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
    @test y(res[10]["O"]) == -6.421
    @test y(res[100]["O"]) == -15.66
    @test countresidues(Model[struc[5], struc[10]]) == 102
    chs = collectchains(Model[struc[5], struc[10]])
    @test length(chs) == 2
    @test chainid.(chs) == ["A", "A"]
    @test z(chs[2][5]["CA"]) == -5.667
    @test countchains(Model[struc[5], struc[10]]) == 2
    mods = collectmodels(Model[struc[10], struc[5]])
    @test length(mods) == 2
    @test modelnumber.(mods) == [5, 10]
    @test z(mods[2]['A'][5]["CA"]) == -5.667
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
    struc = read(IOBuffer(multichar_str), MMCIF)
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
    struc = read(IOBuffer(multlinestruc_str), MMCIF)
    @test coords(struc['A'][1]["CA"]) == [26.091, 52.849, 39.889]
    @test serial(struc['A'][1]["O"]) == 4


    # Test files that should not parse
    @test_throws Exception read(testfilepath("mmCIF", "1AKE_err.cif"), MMCIF)
    @test_throws ErrorException read(testfilepath("mmCIF", "1EN2_err.cif"), MMCIF)


    # Test parsing empty file
    struc = read(IOBuffer(""), MMCIF)
    @test isa(struc, ProteinStructure)
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
    dic_one = MMCIFDict(testfilepath("mmCIF", "1AKE.cif"))
    writemmcif(temp_filename, dic_one)
    dic_two = MMCIFDict(temp_filename)
    @test length(keys(dic_one)) == length(keys(dic_two))
    for k in keys(dic_one)
        @test haskey(dic_two, k)
        @test dic_one[k] == dic_two[k]
    end

    writemmcif(temp_filename, MMCIFDict(Dict("key.key"=> ["value"])))
    @test_throws ArgumentError writemmcif(temp_filename, MMCIFDict(Dict("key"=> ["value"])))
    @test_throws ArgumentError writemmcif(temp_filename, MMCIFDict(Dict("key.key.key"=> ["value"])))
    @test_throws ArgumentError writemmcif(temp_filename, MMCIFDict(Dict(
        "key.one"=> ["value"],
        "key.two"=> ["value1", "value2"],
    )))

    struc = read(testfilepath("PDB", "1SSU.pdb"), PDB)
    writemmcif(temp_filename, struc)
    @test countlines(temp_filename) == 15145
    dic_written = MMCIFDict(temp_filename)
    @test length(keys(dic_written)) == 22
    @test length(dic_written["_atom_site.group_PDB"]) == 15120
    struc_written = read(temp_filename, MMCIF)
    @test isa(struc_written, ProteinStructure)
    @test modelnumbers(struc_written) == collect(1:20)
    @test countatoms(struc_written) == 756
    @test z(struc_written[4]['A']["30"]["OG"]) == -2.177
    @test atomnames(struc_written[15]['A']["39"]) == [
        "N", "CA", "C", "O", "CB", "SG", "H", "HA", "HB2", "HB3"]


    # Test writing to stream
    open(temp_filename, "w") do file
        writemmcif(file, struc)
    end
    @test countlines(temp_filename) == 15145


    # Test selectors
    struc = read(testfilepath("PDB", "1AKE.pdb"), PDB)
    writemmcif(temp_filename, struc, heteroselector)
    @test countlines(temp_filename) == 524
    struc_written = read(temp_filename, MMCIF)
    @test modelnumbers(struc_written) == [1]
    @test countatoms(struc_written) == 492
    @test chainids(struc_written) == ["A", "B"]
    @test tempfactor(struc_written['B']["H_705"]["O"]) == 64.17
    writemmcif(temp_filename, struc, standardselector, disorderselector)
    @test countlines(temp_filename) == 35
    struc_written = read(temp_filename, MMCIF)
    @test countatoms(struc_written) == 5
    @test sum(isdisorderedatom, collectatoms(struc_written)) == 5
    @test defaultaltlocid(struc_written['A'][167]["NH1"]) == 'A'


    # Test writing different element types
    writemmcif(temp_filename, struc[1])
    @test countlines(temp_filename) == 3841
    struc_written = read(temp_filename, MMCIF)
    @test modelnumbers(struc_written) == [1]
    @test countatoms(struc_written) == 3804
    writemmcif(temp_filename, struc['A'])
    @test countlines(temp_filename) == 1991
    struc_written = read(temp_filename, MMCIF)
    @test chainids(struc_written) == ["A"]
    writemmcif(temp_filename, struc['A'][50])
    @test countlines(temp_filename) == 34
    struc_written = read(temp_filename, MMCIF)
    @test chainids(struc_written) == ["A"]
    @test countresidues(struc_written) == 1
    @test atomnames(struc_written['A'][50]) == [
        "N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"]
    writemmcif(temp_filename, struc['A'][50]["CA"])
    @test countlines(temp_filename) == 24
    writemmcif(temp_filename, struc['A'][167]["CZ"])
    @test countlines(temp_filename) == 27
    struc_written = read(temp_filename, MMCIF)
    @test countatoms(struc_written) == 1
    @test isdisorderedatom(collectatoms(struc_written)[1])
    @test length(collectatoms(struc_written)[1]) == 2
    @test tempfactor(collectatoms(struc_written)[1]) == 16.77
    writemmcif(temp_filename, Chain[struc['A'], struc['B']])
    @test countlines(temp_filename) == 3841
    struc_written = read(temp_filename, MMCIF)
    @test chainids(struc_written) == ["A", "B"]
    @test countatoms(struc_written['A']) == 1954
    @test countatoms(struc_written['B']) == 1850
    @test altlocids(struc_written['A']["H_215"]["O1G"]) == ['A', 'B']
    writemmcif(temp_filename, AbstractResidue[struc['A'][51], struc['A'][50]])
    @test countlines(temp_filename) == 42
    struc_written = read(temp_filename, MMCIF)
    @test countresidues(struc_written) == 2
    @test resid.(collectresidues(struc_written)) == ["50", "51"]
    @test countatoms(struc_written) == 17
    writemmcif(temp_filename, AbstractAtom[struc['A'][51]["CA"], struc['A'][50]["CA"]])
    @test countlines(temp_filename) == 27
    struc_written = read(temp_filename, MMCIF)
    @test countatoms(struc_written) == 2
    @test !ishetero(struc_written['A'][51]["CA"])


    # Test multiple model writing
    struc = read(testfilepath("PDB", "1SSU.pdb"), PDB)
    writemmcif(temp_filename, Model[struc[10], struc[5]])
    @test countlines(temp_filename) == 1537
    struc_written = read(temp_filename, MMCIF)
    @test modelnumbers(struc_written) == [5, 10]
    @test modelnumber(defaultmodel(struc_written)) == 5
    @test countatoms(struc_written[5]) == 756
    @test countatoms(struc_written[10]) == 756
    @test_throws KeyError struc_written[1]


    # Test disordered residue writing
    struc = read(testfilepath("PDB", "1EN2.pdb"), PDB)
    writemmcif(temp_filename, struc)
    @test countlines(temp_filename) == 844
    struc_written = read(temp_filename, MMCIF)
    @test countatoms(struc_written) == 754
    @test isa(struc_written['A'][15], Residue)
    @test isa(struc_written['A'][16], DisorderedResidue)
    @test defaultresname(struc_written['A'][16]) == "ARG"
    @test isa(struc_written['A'][16]["N"], DisorderedAtom)
    @test defaultaltlocid(struc_written['A'][16]["N"]) == 'A'
    @test isa(disorderedres(struc_written['A'][16], "TRP")["N"], Atom)
    @test countatoms(struc_written['A'][16]) == 11
    @test countatoms(disorderedres(struc_written['A'][16], "TRP")) == 14
    writemmcif(temp_filename, AbstractResidue[struc['A'][16], struc['A'][10]])
    @test countlines(temp_filename) == 71
    struc_written = read(temp_filename, MMCIF)
    @test countresidues(struc_written) == 2
    @test isa(struc_written['A'][10], DisorderedResidue)
    @test isa(struc_written['A'][16], DisorderedResidue)
    @test defaultresname(struc_written['A'][10]) == "SER"
    @test isa(disorderedres(struc_written['A'][10], "GLY")["O"], Atom)
    @test altlocid(disorderedres(struc_written['A'][10], "GLY")["O"]) == 'B'
    @test countatoms(struc_written['A'][10]) == 6
    @test countatoms(struc_written['A'][16]) == 11


    # Test writing multi-character chain IDs
    struc = read(IOBuffer(multichar_str), MMCIF)
    writemmcif(temp_filename, struc)
    struc_back = read(temp_filename, MMCIF)
    @test chainids(struc_back) == ["A1", "A2"]
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

    struc_1AKE = read(testfilepath("PDB", "1AKE.pdb"), PDB)
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

    struc_1SSU = read(testfilepath("PDB", "1SSU.pdb"), PDB)
    @test isapprox(rmsd(struc_1SSU[1], struc_1SSU[2], calphaselector), 4.1821925809691889)
    @test isapprox(rmsd(struc_1SSU[5], struc_1SSU[6], backboneselector), 5.2878196391279939)
    @test_throws ArgumentError rmsd(struc_1SSU[1]['A'][8], struc_1SSU[1]['A'][9])


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

    disps = displacements(struc_1SSU[5], struc_1SSU[10])
    @test isa(disps, Vector{Float64})
    @test length(disps) == 756
    @test isapprox(disps[20], sqrt(1.984766))
    disps = displacements(struc_1SSU[5], struc_1SSU[10], calphaselector)
    @test length(disps) == 51
    @test isapprox(disps[20], sqrt(0.032822))


    # Test sqdistance and distance
    at_a = Atom(100, "CA", ' ', [1.0, 2.0, 3.0], 1.0, 10.0, " C", "  ", res)
    at_b = Atom(110, "CA", ' ', [0.0, -1.0, 3.0], 1.0, 10.0, " C", "  ", res)
    @test sqdistance(at_a, at_b) == 10.0
    @test isapprox(distance(at_a, at_b), sqrt(10))

    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B']), sqrt(6.852947))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'][50]), sqrt(530.645746))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'][50]["CA"]), sqrt(574.699125))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'], backboneselector), sqrt(17.350083))
    @test isapprox(distance(struc_1AKE['A'], struc_1AKE['B'], standardselector), sqrt(11.252973))
    @test isapprox(distance(struc_1AKE['A'][50]["CA"], struc_1AKE['B'][50]["CA"]), sqrt(2607.154834))


    # Test bondangle
    at_a = Atom(100, "CA", ' ', [1.0, 0.0, 1.0], 1.0, 10.0, " C", "  ", res)
    at_b = Atom(100, "CA", ' ', [0.0, 0.0, 0.0], 1.0, 10.0, " C", "  ", res)
    at_c = Atom(100, "CA", ' ', [3.0, 2.0, 1.0], 1.0, 10.0, " C", "  ", res)
    @test isapprox(bondangle(at_a, at_b, at_c), 0.713724378944765)
    vec_a = [2.0, 0.0, 0.0]
    vec_b = [2.0, 1.0, 1.0]
    @test isapprox(bondangle(vec_a, vec_b), 0.615479708670387)


    # Test dihedral functions
    at_a = Atom(100, "CA", ' ', [-1.0, -1.0, 0.0], 1.0, 10.0, " C", "  ", res)
    at_b = Atom(100, "CA", ' ', [0.0, 0.0, 0.0], 1.0, 10.0, " C", "  ", res)
    at_c = Atom(100, "CA", ' ', [1.0, 0.0, 0.0], 1.0, 10.0, " C", "  ", res)
    at_d = Atom(100, "CA", ' ', [2.0, 1.0, -1.0], 1.0, 10.0, " C", "  ", res)
    @test isapprox(dihedralangle(at_a, at_b, at_c, at_d), 2.356194490192345)
    vec_a = [1.0, 1.0, 0.0]
    vec_b = [1.0, 0.0, 0.0]
    vec_c = [1.0, -1.0, 1.0]
    @test isapprox(dihedralangle(vec_a, vec_b, vec_c), -0.785398163397448)

    @test isapprox(omegaangle(struc_1AKE['A'][20], struc_1AKE['A'][19]), -3.091913621551854, atol=1e-5)
    @test isapprox(phiangle(struc_1AKE['A'][7], struc_1AKE['A'][6]), 2.851151641716221, atol=1e-5)
    @test isapprox(psiangle(struc_1AKE['A'][8], struc_1AKE['A'][9]), 2.838265381719911, atol=1e-5)
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
    @test isapprox(phis[5], -1.764512005880236, atol=1e-5)
    @test isapprox(psis[10], 0.4425094841355222, atol=1e-5)
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
    ])
    @test isapprox(DistanceMap(struc_1AKE[1]).data, [
        0.0                2.6178134005310594
        2.6178134005310594 0.0
    ])
    dmap = DistanceMap(struc_1AKE['A'])
    @test size(dmap) == (456, 456)
    @test size(dmap, 2) == 456
    @test isapprox(dmap[196, 110], 3.0188665091388214)
    dmap[196, 110] = 10.0
    @test dmap[196, 110] == 10.0
    show(devnull, dmap)

    @test isapprox(DistanceMap(struc_1AKE['A'][10], struc_1AKE['A'][11]).data, [
        2.8803883418733673 4.295162278657233  5.211988200293629  6.427722458227334 5.03904008715946
        2.5066699024801813 3.850861981427013  4.557738913101538  5.762204526047301 5.0027798272560435
        1.3632875705440877 2.4469166311911783 3.1850982716393537 4.33450689236965  3.799490623754715
        2.2682515292621335 2.774401016435799  3.235280049702036  4.198412795331114 4.313454647959104
    ])

    # Test the plot recipe
    RecipesBase.apply_recipe(Dict{Symbol, Any}(), dmap)
end

# Delete temporary file
rm(temp_filename)

end # TestBioStructures
