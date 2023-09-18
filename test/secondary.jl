using BioStructures
using Test

@testset "Secondary Structure Prediction" begin
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
