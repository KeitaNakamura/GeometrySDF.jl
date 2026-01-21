# GeometrySDF.jl

*Fast signed distance field (SDF) generation from 3D triangle meshes and 2D polylines in Julia*

```julia
julia> using GeometrySDF

julia> sdf = GeometrySDF.sdf("test/stanford_dragon.stl");

julia> GeometrySDF.write_vtk("dragon", sdf)
1-element Vector{String}:
 "dragon.vti"
```
