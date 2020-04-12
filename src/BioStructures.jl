# BioStructures.jl
# ================
#
# A Julia package to read, write and manipulate macromolecular structures.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioStructures.jl/blob/master/LICENSE.md

"Read, write and manipulate macromolecular structures."
module BioStructures

using LinearAlgebra
using Statistics

using CodecZlib
using Format
using RecipesBase
using LightGraphs
using MetaGraphs
using DataFrames
using BioCore
using BioSymbols
using BioSequences
using BioAlignments
import MMTF: parsemmtf, writemmtf # Imported to avoid clash with MMTF name

include("model.jl")
include("pdb.jl")
include("mmcif.jl")
include("mmtf.jl")
include("spatial.jl")

end # BioStructures
