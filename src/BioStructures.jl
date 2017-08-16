# BioStructures.jl
# ================
#
# A julia package to read, write and manipulate macromolecular structures.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioStructures.jl/blob/master/LICENSE

__precompile__()

module BioStructures

using Libz
import BioCore
import BioCore.distance
import BioSymbols
import BioSequences.AminoAcidSequence

include("model.jl")
include("pdb.jl")
include("spatial.jl")

end # BioStructures
