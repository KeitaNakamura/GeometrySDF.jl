struct SignedDistanceField{dim, T, V, A <: AbstractArray{T}}
    grid::Grid{dim, T, V}
    ϕ::A
end

function sdf(
        prims::AbstractVector{<: Union{Triangle{dim, T}, Segment{dim, T}}}, normals::AbstractVector{<: StaticVector{dim, T}}, grid::Grid{dim, T};
        progress::Bool = false,
    ) where {dim, T}
    @assert length(prims) == length(normals)

    ϕ = Array{T}(undef, size(grid))
    p = ProgressMeter.Progress(length(grid); desc = "Generating SDF...", enabled=progress)

    bvh = build_bvh(prims; leaf_size=8)
    atol = sqrt(eps(T)) * spacing(grid)^2
    rtol = sqrt(eps(T))

    Threads.@threads for I in eachindex(grid)
        x = grid[I]
        d²_min = T(Inf)
        dₙ_max = T(0)

        stack = [1]
        while !isempty(stack)
            node = getnode(bvh, pop!(stack))

            dist2_point_aabb(x, node) > d²_min && continue

            if isleaf(node)
                @inbounds for i in leaf_prims(bvh, node)
                    prim = prims[i]
                    n = normals[i]

                    v = x - closest_point(x, prim)
                    d² = v ⋅ v
                    dₙ = v ⋅ n

                    if abs(d² - d²_min) < max(atol, rtol * d²_min)
                        dₙ_max = ifelse(abs(dₙ) > abs(dₙ_max), dₙ, dₙ_max)
                    elseif d² < d²_min
                        d²_min = d²
                        dₙ_max = dₙ
                    end
                end
            else
                dl = dist2_point_aabb(x, left_child(bvh, node))
                dr = dist2_point_aabb(x, right_child(bvh, node))
                if dl < dr
                    dr ≤ d²_min && push!(stack, node.right)
                    dl ≤ d²_min && push!(stack, node.left)
                else
                    dl ≤ d²_min && push!(stack, node.left)
                    dr ≤ d²_min && push!(stack, node.right)
                end
            end
        end
        ϕ[I] = sign(dₙ_max) * sqrt(d²_min)

        ProgressMeter.next!(p)
    end
    ProgressMeter.finish!(p)

    SignedDistanceField(grid, ϕ)
end

# Book: Real-Time Collision Detection
@inline function closest_point(p, a, b, c)
    # check if P in vertex region outside A
    ab = b - a
    ac = c - a
    ap = p - a
    d1 = ab ⋅ ap
    d2 = ac ⋅ ap
    (d1 ≤ 0 && d2 ≤ 0) && return a

    # check if P in vertex region outside B
    bp = p - b
    d3 = ab ⋅ bp
    d4 = ac ⋅ bp
    (d3 ≥ 0 && d4 ≤ d3) && return b

    # check if P in edge region of AB, if so return projection of P onto AB
    vc = d1*d4 - d3*d2
    if vc ≤ 0 && d1 ≥ 0 && d3 ≤ 0
        v = d1 / (d1 - d3)
        return a + v * ab
    end

    # check if P in vertex region outside C
    cp = p - c
    d5 = ab ⋅ cp
    d6 = ac ⋅ cp
    (d6 ≥ 0 && d5 ≤ d6) && return c

    # check if P in edge region of AC, if so return projection of P onto AC
    vb = d5*d2 - d1*d6
    if vb ≤ 0 && d2 ≥ 0 && d6 ≤ 0
        w = d2 / (d2 - d6)
        return a + w * ac
    end

    # check if P in edge region of BC, if so return projection of P onto BC
    va = d3*d6 - d5*d4
    if va ≤ 0 && (d4 - d3) ≥ 0 && (d5 - d6) ≥ 0
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return b + w * (c - b)
    end

    # P inside face region. compute Q through its barycentric coordinates
    denom = inv(va + vb + vc)
    v = vb * denom
    w = vc * denom
    return a + ab*v + ac*w
end

function measure(sdf::SignedDistanceField{dim}) where {dim}
    grid, ϕ = sdf.grid, sdf.ϕ
    l = spacing(grid)
    H = h -> heaviside_function(h, l)
    l^dim * mapreduce((ϕ,x)->H(-ϕ), +, ϕ, grid)
end
area(sdf::SignedDistanceField{2}) = measure(sdf)
volume(sdf::SignedDistanceField{3}) = measure(sdf)

function centroid(sdf::SignedDistanceField{dim}) where {dim}
    grid, ϕ = sdf.grid, sdf.ϕ
    l = spacing(grid)
    H = h -> heaviside_function(h, l)
    l^dim * mapreduce((ϕ,x)->H(-ϕ)*SVector(x), +, ϕ, grid) / measure(sdf)
end

# per density
function inertia(sdf::SignedDistanceField{dim}; origin = centroid(sdf)) where {dim}
    grid, ϕ = sdf.grid, sdf.ϕ
    l = spacing(grid)
    H = h -> heaviside_function(h, l)
    l^dim * mapreduce((ϕ,x) -> H(-ϕ) * inertia_kernel(x, origin), +, ϕ, grid)
end

function inertia_kernel(x::SVector{3}, c::SVector{3})
    I₁₁ = (x[2]-c[2])^2 + (x[3]-c[3])^2
    I₂₂ = (x[1]-c[1])^2 + (x[3]-c[3])^2
    I₃₃ = (x[1]-c[1])^2 + (x[2]-c[2])^2
    I₂₃ = I₃₂ = -(x[2]-c[2]) * (x[3]-c[3])
    I₁₃ = I₃₁ = -(x[1]-c[1]) * (x[3]-c[3])
    I₁₂ = I₂₁ = -(x[1]-c[1]) * (x[2]-c[2])
    @SMatrix [I₁₁ I₁₂ I₁₃
              I₂₁ I₂₂ I₂₃
              I₃₁ I₃₂ I₃₃]
end

@inline function inertia_kernel(x::SVector{2}, c::SVector{2})
    dx = x[1] - c[1]
    dy = x[2] - c[2]
    dx*dx + dy*dy
end

function heaviside_function(ϕ::T, Δ::T, δ::T=T(1.5)) where {T <: Real}
    ϵ = δ * Δ
    ξ = ϕ / ϵ
    ξ < -1 && return T(0)
    ξ >  1 && return T(1)
    (1 + ξ + sin(π*ξ)/π) / 2
end
