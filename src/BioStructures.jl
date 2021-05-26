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

using BioAlignments
import BioCore # Imported to avoid clash with BioGenerics distance
using BioGenerics
using BioSequences
using BioSymbols
using CodecZlib
using DataFrames
using Downloads
using Format
using LightGraphs
using MetaGraphs
import MMTF: parsemmtf, writemmtf # Imported to avoid clash with MMTF name
using RecipesBase

include("model.jl")
include("pdb.jl")
include("mmcif.jl")
include("mmtf.jl")
include("spatial.jl")

end # BioStructures
