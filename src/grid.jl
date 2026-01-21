struct Grid{dim, T, V <: AbstractVector{T}} <: AbstractArray{SVector{dim, T}, dim}
    h::T
    axes::NTuple{dim, V}
end

Base.size(grid::Grid) = map(length, grid.axes)
spacing(grid::Grid) = grid.h

@generated function Base.getindex(grid::Grid{dim}, I::Vararg{Int, dim}) where {dim}
    quote
        @inline
        @boundscheck checkbounds(grid, I...)
        @inbounds SVector(Base.Cartesian.@ntuple $dim d -> grid.axes[d][I[d]])
    end
end
