# GeometrySDF.jl

*Fast signed distance field (SDF) generation from 3D triangle meshes and 2D polylines in Julia*

[![CI](https://github.com/KeitaNakamura/GeometrySDF.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/KeitaNakamura/GeometrySDF.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/KeitaNakamura/GeometrySDF.jl/graph/badge.svg?token=hKUB58zHkF)](https://codecov.io/gh/KeitaNakamura/GeometrySDF.jl)
[![Stable](https://img.shields.io/badge/docs-latest%20release-blue.svg)](https://KeitaNakamura.github.io/GeometrySDF.jl/stable)


```julia
julia> using GeometrySDF

julia> sdf3d = GeometrySDF.sdf("test/stanford_bunny.stl");

julia> GeometrySDF.writevtk("bunny", sdf3d)
1-element Vector{String}:
 "bunny.vti"

julia> sdf2d = GeometrySDF.sdf([0.0 1.0; 0.2886751346 0.5; 0.8660254038 0.5; 0.5773502692 0.0; 0.8660254038 -0.5; 0.2886751346 -0.5; 0.0 -1.0; -0.2886751346 -0.5; -0.8660254038 -0.5; -0.5773502692 0.0; -0.8660254038 0.5; -0.2886751346 0.5]');

julia> GeometrySDF.writevtk("hexagram", sdf2d)
1-element Vector{String}:
 "hexagram.vti"
```

<h1 align="center">
    <img width="180" src="https://github.com/user-attachments/assets/706348b6-c149-4688-9894-84d949d82f4b" />
    <img width="210" src="https://github.com/user-attachments/assets/9408b4e9-c812-486d-bece-8c8f8c33af8e" />
    <img width="170" src="https://github.com/user-attachments/assets/806915fd-a750-4eac-8946-3ba778258e5e" />
    <img width="150" src="https://github.com/user-attachments/assets/71a89afe-54d7-4732-ada0-d744b496a5d0" />
</h1>
