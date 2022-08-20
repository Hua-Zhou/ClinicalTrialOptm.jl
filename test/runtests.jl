using ClinicalTrialOptm
using Distributions, ForwardDiff, Random, Test

@testset "Country" begin
    m = 1.5
    s² = 2.0
    l, u = 0, 10
    c₀ = 20_000.0
    c = 5_000.0
    q = 1_500.0
    d = 0.05
    T₀ = Uniform(0.0, 6.0)
    Td = 24.0
    ctry = Country(m, s², l, u, c₀, c, q, d, T₀, Td)
    @test ctry.m == m
    @info "mean/var of a center"
    @show mean(ctry)
    @show var(ctry)
    @info "pgf of a center"
    @test pgf(ctry, 1.0) == 1.0
    # derivative of pgf at 1.0 == mean
    @test ForwardDiff.derivative(z -> pgf(ctry, z), 1.0) ≈ mean(ctry)
    # 2nd derivative of pgf at 1.0 == E[X(X-1)]
    @test ForwardDiff.derivative(
        t -> ForwardDiff.derivative(z -> pgf(ctry, z), t), 
        1.0
        ) ≈ var(ctry) + mean(ctry)^2 - mean(ctry)
end

@testset "ClinicalTrial" begin
    rng = MersenneTwister(123)
    nc = 5 # number of countries
    m = rand(rng, Uniform(1.0, 2.0), nc)
    s² = rand(rng, Uniform(1.0, 2.0), nc)
    l = rand(rng, 0:2, nc)
    u = l .+ rand(rng, 10:30, nc)
    c₀ = rand(rng, 15000.0:1000.0:25000.0, nc)
    c = rand(rng, 4000.0:200.0:5000.0, nc)
    q = rand(rng, 1000.0:200.0:2000.0, nc)
    d = rand(rng, 0.01:0.01:0.15, nc)
    T₀ = fill(Uniform(0.0, 6.0), nc)
    Td = 24.0
    centers = rand(rng, 10:20, nc)
    ct = ClinicalTrial(m, s², l, u, c₀, c, q, d, T₀, Td, centers)
    @info "mean/var of a clinical trial"
    @show mean(ct)
    @show var(ct)
    @info "pgf of a clinical trial"
    @test pgf(ct, 1.0) == 1.0
    # derivative of pgf at 1.0 == mean
    @test ForwardDiff.derivative(z -> pgf(ct, z), 1.0) ≈ mean(ct)
    # 2nd derivative of pgf at 1.0 == E[X(X-1)]
    @test ForwardDiff.derivative(
        t -> ForwardDiff.derivative(z -> pgf(ct, z), t), 
        1.0
        ) ≈ var(ct) + mean(ct)^2 - mean(ct)
end
