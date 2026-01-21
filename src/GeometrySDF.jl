module GeometrySDF

using LinearAlgebra
using StaticArrays
using GeometryBasics

using FileIO
import WriteVTK
import ProgressMeter

public SignedDistanceField, sdf

include("grid.jl")
include("bvh.jl")
include("sdf.jl")
include("fileio.jl")

end # module GeometrySDF
