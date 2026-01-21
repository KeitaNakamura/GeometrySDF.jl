# 3D version
function sdf(path::String; datatype::Type{T}=Float64, kwargs...) where {T}
    mesh = load(path; pointtype=Point{3, T})
    _sdf(mesh; datatype, kwargs...)
end

# 2D version
function sdf(points::AbstractMatrix{T}; kwargs...) where {T <: Real}
    @assert size(points, 1) == 2
    sdf(vec(reinterpret(SVector{2, T}, points)); kwargs...)
end
function sdf(points::AbstractVector{<: SVector{2}}; datatype::Type{T}=Float64, kwargs...) where {T}
    _sdf(convert(Vector{SVector{2, T}}, points); datatype, kwargs...)
end

function _sdf(
        mesh; spacing=:auto, minres=128,
        padding::Int=10, datatype::Type{T}=Float64, progress::Bool=false
    ) where {T}

    prims = mesh_primitives(mesh)
    normals = mesh_normals(mesh)
    lims = mesh_bbox(mesh)

    # grid
    h = spacing === :auto ? T(auto_spacing_bbox(lims; minres)) : T(spacing)
    axes = map(lims) do (xmin, xmax)
        range(T(xmin-padding*h), T(xmax+padding*h); step=h)
    end
    grid = Grid(h, axes)

    sdf(prims, normals, grid; progress)
end

function auto_spacing_bbox(lims; minres::Int)
    Lmax = maximum((xmax - xmin) for (xmin, xmax) in lims)
    return Lmax / minres
end

function mesh_bbox(mesh::Mesh{dim, T, <: NgonFace{3}}) where {dim, T}
    position = coordinates(mesh)
    NTuple{dim, Tuple{T,T}}(extrema(reinterpret(reshape, T, position), dims=2))
end
function mesh_bbox(position::AbstractVector{<: SVector{dim, T}}) where {dim, T}
    NTuple{dim, Tuple{T,T}}(extrema(reinterpret(reshape, T, position), dims=2))
end

function mesh_primitives(mesh::Mesh{dim, T, <: NgonFace{3}}) where {dim, T}
    # create a vector of `Triangle` because `mesh[i]` is slow
    collect(Triangle{dim, T}, mesh)
end
function mesh_primitives(points::AbstractVector{<: SVector})
    segs = map(p -> Segment(points[p], points[p+1]), 1:length(points)-1)
    push!(segs, Segment(points[end], points[1]))
    return segs
end

function mesh_normals(mesh::Mesh{3, <: Any, <: NgonFace{3}})
    # computation is duplicated if normals are already stored in the mesh
    map(mesh) do tri
        a, b, c = tri
        n = (b-a) × (c-a)
        iszero(n) ? zero(n) : normalize(n)
    end
end
function mesh_normals(points::AbstractVector{<: SVector})
    normals = map(p -> outward_normal(Segment(points[p], points[p+1])), 1:length(points)-1)
    push!(normals, outward_normal(Segment(points[end], points[1])))
    normals
end

function write_vtk(path::String, sdf::SignedDistanceField)
    grid, ϕ = sdf.grid, sdf.ϕ
    vtk = WriteVTK.vtk_grid(path, grid.axes...)
    vtk["Signed distance"] = ϕ
    WriteVTK.vtk_save(vtk)
end
