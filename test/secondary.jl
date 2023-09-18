using BioStructures
using Test

@testset "Secondary Structure Information" begin
    @testset "Residue Manipulation 1" begin
        res = Residue("ALA", 1, ' ', false, Chain('A'))
        ss_code!(res, "H")
        @test ss_code(res) == "H"

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
            "T"
        )
        @test ss_code(res) == "T"

        for at in atoms(res)
            @test ss_code(at) == "T"
        end
    end

    @testset "DSSP" begin
        pdb_path = downloadpdb("1BQ0", dir = tempname() * ".cif", format = MMCIF)
        struc = read(pdb_path, MMCIF)
        struc = run_dssp(struc)

        helix_atoms = collectatoms(struc, helixselector)
        sheet_atoms = collectatoms(struc, sheetselector)
        coil_atoms = collectatoms(struc, coilselector)

        @test length(helix_atoms) > 0
        @test length(sheet_atoms) == 0
        @test length(coil_atoms) > 0

        helix_residues = collectresidues(struc, helixselector)
        sheet_residues = collectresidues(struc, sheetselector)
        coil_residues = collectresidues(struc, coilselector)

        @test length(helix_residues) > 0
        @test length(sheet_residues) == 0
        @test length(coil_residues) > 0

        rm(pdb_path)
    end

    @testset "STRIDE" begin
        pdb_path = downloadpdb("1BQ0", dir = tempname() * ".cif", format = MMCIF)
        struc = read(pdb_path, MMCIF)
        struc = run_stride(struc)

        helix_atoms = collectatoms(struc, helixselector)
        sheet_atoms = collectatoms(struc, sheetselector)
        coil_atoms = collectatoms(struc, coilselector)

        @test length(helix_atoms) > 0
        @test length(sheet_atoms) == 0
        @test length(coil_atoms) > 0

        helix_residues = collectresidues(struc, helixselector)
        sheet_residues = collectresidues(struc, sheetselector)
        coil_residues = collectresidues(struc, coilselector)

        @test length(helix_residues) > 0
        @test length(sheet_residues) == 0
        @test length(coil_residues) > 0

        rm(pdb_path)
    end
end
