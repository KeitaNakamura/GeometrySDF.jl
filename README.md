# GeometrySDF.jl

*Fast signed distance field (SDF) generation from 3D triangle meshes and 2D polylines in Julia*

[![CI](https://github.com/KeitaNakamura/GeometrySDF.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/KeitaNakamura/GeometrySDF.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/KeitaNakamura/GeometrySDF.jl/graph/badge.svg?token=hKUB58zHkF)](https://codecov.io/gh/KeitaNakamura/GeometrySDF.jl)


```julia
julia> using GeometrySDF

julia> sdf = GeometrySDF.sdf("test/stanford_dragon.stl");

julia> GeometrySDF.write_vtk("dragon", sdf)
1-element Vector{String}:
 "dragon.vti"
```

<img width="400" src="https://github.com/user-attachments/assets/cd523eaa-909e-4d03-ac74-1e6a52fcb776" />
