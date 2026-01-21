using GeometrySDF
using Test

using StaticArrays

@testset "GeometrySDF 3D" begin
    @testset "Sphere" begin
        sdf = GeometrySDF.sdf("sphere.stl")
        r = 0.5
        c = [0.1,0.2,0.3]
        ρ = 1
        V = (4/3)*π*r^3
        m = ρ*V
        I = (2/5)*m*r^2 * [1 0 0; 0 1 0; 0 0 1]
        @test GeometrySDF.volume(sdf) ≈ V rtol=1e-2
        @test GeometrySDF.centroid(sdf) ≈ c rtol=1e-2
        @test GeometrySDF.inertia(sdf) ≈ I rtol=1e-2
    end
    @testset "Cone" begin
        sdf = GeometrySDF.sdf("cone.stl")
        r = 1000.0
        h = 2000.0
        c = [1000,1000,3h/4]
        ρ = 1
        V = π*r^2*h/3
        m = ρ*V
        I = [(3/80)*m*(4r^2+h^2) 0 0; 0 (3/80)*m*(4r^2+h^2) 0; 0 0 (3/10)*m*r^2]
        @test GeometrySDF.volume(sdf) ≈ V rtol=1e-2
        @test GeometrySDF.centroid(sdf) ≈ c rtol=1e-2
        @test GeometrySDF.inertia(sdf) ≈ I rtol=1e-2
        # check inertia tensor for other axes
        @test GeometrySDF.inertia(sdf; origin=SVector{3,Float32}(1000,1000,0)) ≈ [(3/5)*m*h^2+(3/20)*m*r^2 0 0; 0 (3/5)*m*h^2+(3/20)*m*r^2 0; 0 0 (3/10)*m*r^2] rtol=1e-2
        @test GeometrySDF.inertia(sdf; origin=SVector{3,Float32}(1000,1000,h)) ≈ [(1/10)*m*h^2+(3/20)*m*r^2 0 0; 0 (1/10)*m*h^2+(3/20)*m*r^2 0; 0 0 (3/10)*m*r^2] rtol=1e-2
    end
    @testset "Tube" begin
        sdf = GeometrySDF.sdf("tube.stl")
        r₁ = 500.0
        r₂ = 1000.0
        h = 2000.0
        c = [1000,1000,h/2]
        ρ = 1
        V = π*(r₂^2-r₁^2)*h
        m = ρ*V
        I = (1/12)*m * [3(r₂^2+r₁^2)+h^2 0 0; 0 3(r₂^2+r₁^2)+h^2 0; 0 0 6(r₂^2+r₁^2)]
        @test GeometrySDF.volume(sdf) ≈ V rtol=1e-2
        @test GeometrySDF.centroid(sdf) ≈ c rtol=1e-2
        @test GeometrySDF.inertia(sdf) ≈ I rtol=1e-2
    end
end

@testset "GeometrySDF 2D" begin
    sdf = GeometrySDF.sdf([0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]')
    @test GeometrySDF.area(sdf) ≈ 1.0 rtol=1e-2
    @test GeometrySDF.centroid(sdf) ≈ SVector(0.5, 0.5) rtol=1e-2
    @test GeometrySDF.inertia(sdf) ≈ 1/6 rtol=1e-2  # per density, about centroid
end
