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
using BioGenerics
using BioSequences
using BioSymbols
using CodecZlib
using Downloads
using Format
import MMTF: parsemmtf, writemmtf # Imported to avoid clash with MMTF name
using RecipesBase
using STRIDE_jll
using DSSP_jll

include("model.jl")
include("pdb.jl")
include("mmcif.jl")
include("mmtf.jl")
include("spatial.jl")
include("secondary.jl")

if !isdefined(Base, :get_extension)
    include("../ext/BioStructuresDataFramesExt.jl")
    include("../ext/BioStructuresGraphsExt.jl")
end

end # BioStructures
