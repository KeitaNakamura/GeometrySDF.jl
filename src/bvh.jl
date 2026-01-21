struct Segment{dim, T}
    a::SVector{dim, T}
    b::SVector{dim, T}
end

@inline centroid(seg::Segment) = (seg.a + seg.b) / 2

@inline function closest_point(p::SVector{dim, T}, seg::Segment{dim, T}) where {dim, T}
    a, b = seg.a, seg.b
    ab = b - a
    denom = ab ⋅ ab
    iszero(denom) && return a
    t = ((p - a) ⋅ ab) / denom
    t = clamp(t, zero(T), one(T))
    a + t*ab
end

@inline function outward_normal(seg::Segment{2, T}) where {T}
    t = seg.b - seg.a
    n = SVector(t[2], -t[1])
    nn = n ⋅ n
    iszero(nn) && return zero(SVector{2, T})
    n * inv(sqrt(nn))
end

@inline function centroid(tri::Triangle)
    a, b, c = tri
    (a + b + c) / 3
end

@inline function closest_point(p::SVector, tri::Triangle)
    closest_point(p, tri...)
end

struct AABB{dim, T}
    bmin::SVector{dim, T}
    bmax::SVector{dim, T}
end

@inline function Base.merge(a::AABB{dim, T}, b::AABB{dim, T}) where {dim, T}
    AABB{dim, T}(min.(a.bmin, b.bmin), max.(a.bmax, b.bmax))
end

@inline function aabb(seg::Segment{dim, T}) where {dim, T}
    a, b = seg.a, seg.b
    AABB{dim, T}(min.(a,b), max.(a,b))
end

@inline function aabb(tri::Triangle{dim, T}) where {dim, T}
    a, b, c = tri
    bmin = min.(min.(a, b), c)
    bmax = max.(max.(a, b), c)
    AABB{dim, T}(bmin, bmax)
end

@inline function dist2_point_aabb(x::SVector{dim, T}, box::AABB{dim, T}) where {dim, T}
    s = zero(T)
    @inbounds for k in 1:dim
        if x[k] < box.bmin[k]
            d = box.bmin[k] - x[k]
            s += d*d
        elseif x[k] > box.bmax[k]
            d = x[k] - box.bmax[k]
            s += d*d
        end
    end
    return s
end

struct BVHNode{dim, T}
    aabb::AABB{dim, T}
    left::Int
    right::Int
    start::Int # range into `indices`
    count::Int # 0 => internal, >0 => leaf
end

isleaf(node::BVHNode) = node.count > 0
leaf_range(node::BVHNode) = node.start : node.start + node.count - 1
@inline dist2_point_aabb(x::SVector, node::BVHNode) = dist2_point_aabb(x, node.aabb)

struct BVH{dim,T}
    nodes::Vector{BVHNode{dim,T}}
    indices::Vector{Int} # primitive indices, permuted by build
    aabbs::Vector{AABB{dim, T}}
    centroids::Vector{SVector{dim, T}}
end

getnode(bvh::BVH, i::Int) = bvh.nodes[i]
leaf_prims(bvh::BVH, node::BVHNode) = view(bvh.indices, leaf_range(node))
left_child(bvh::BVH, node::BVHNode) = bvh.nodes[node.left]
right_child(bvh::BVH, node::BVHNode) = bvh.nodes[node.right]

@inline centroid_axis(centroids, i, axis) = centroids[i][axis]

@inline function median3!(indices::Vector{Int}, lo::Int, mid::Int, hi::Int,
                          axis::Int, centroids)
    i_lo  = indices[lo]
    i_mid = indices[mid]
    i_hi  = indices[hi]

    a_lo  = centroid_axis(centroids, i_lo,  axis)
    a_mid = centroid_axis(centroids, i_mid, axis)
    a_hi  = centroid_axis(centroids, i_hi,  axis)

    # return index position (lo/mid/hi) whose key is median
    if a_lo < a_mid
        if     a_mid < a_hi; return mid
        elseif a_lo  < a_hi; return hi
        else               ; return lo
        end
    else
        if     a_lo  < a_hi; return lo
        elseif a_mid < a_hi; return hi
        else               ; return mid
        end
    end
end

@inline function partition!(indices::Vector{Int}, lo::Int, hi::Int,
                            axis::Int, pivot_key, centroids)
    i = lo
    j = hi
    while true
        while centroid_axis(centroids, indices[i], axis) < pivot_key
            i += 1
        end
        while centroid_axis(centroids, indices[j], axis) > pivot_key
            j -= 1
        end
        if i >= j
            return j
        end
        indices[i], indices[j] = indices[j], indices[i] # swap
        i += 1
        j -= 1
    end
end

function select_kth!(indices::Vector{Int}, lo::Int, hi::Int, k::Int,
                     axis::Int, centroids)
    while lo < hi
        mid = (lo + hi) ÷ 2
        ppos = median3!(indices, lo, mid, hi, axis, centroids)
        pivot_key = centroid_axis(centroids, indices[ppos], axis)

        cut = partition!(indices, lo, hi, axis, pivot_key, centroids)

        if k <= cut
            hi = cut
        else
            lo = cut + 1
        end
    end
end

function build_bvh(prims::AbstractVector{<: P}; leaf_size::Int=8) where {dim, T, P <: Union{Segment{dim, T}, Triangle{dim, T}}}
    n = length(prims)

    aabbs = Vector{AABB{dim, T}}(undef, n)
    centroids = Vector{SVector{dim, T}}(undef, n)
    @inbounds for i in 1:n
        prim = prims[i]
        aabbs[i] = aabb(prim)
        centroids[i] = centroid(prim)
    end

    indices = collect(1:n)
    nodes = BVHNode{dim, T}[]
    sizehint!(nodes, 2n)

    if n == 0
        return BVH{dim, T}(nodes, indices, aabbs, centroids)
    end

    function node_aabb(lo::Int, hi::Int)
        box = aabbs[indices[lo]]
        @inbounds for j in (lo+1):hi
            box = merge(box, aabbs[indices[j]])
        end
        return box
    end

    function centroid_bounds(lo::Int, hi::Int)
        c0 = centroids[indices[lo]]
        cmin = c0; cmax = c0
        @inbounds for j in (lo+1):hi
            c = centroids[indices[j]]
            cmin = min.(cmin, c)
            cmax = max.(cmax, c)
        end
        return cmin, cmax
    end

    function build_node(lo::Int, hi::Int)
        box = node_aabb(lo, hi)
        count = hi - lo + 1

        if count <= leaf_size
            push!(nodes, BVHNode{dim, T}(box, 0, 0, lo, count))
            return length(nodes)
        end

        cmin, cmax = centroid_bounds(lo, hi)
        ext = cmax - cmin
        _, axis = findmax(ext) # axis::Int

        mid = (lo + hi) ÷ 2
        select_kth!(indices, lo, hi, mid, axis, centroids)

        push!(nodes, BVHNode{dim, T}(box, 0, 0, 0, 0)) # placeholder
        idx = length(nodes)

        left  = build_node(lo, mid)
        right = build_node(mid+1, hi)
        nodes[idx] = BVHNode{dim, T}(box, left, right, 0, 0)

        return idx
    end

    build_node(1, n)
    return BVH{dim, T}(nodes, indices, aabbs, centroids)
end
