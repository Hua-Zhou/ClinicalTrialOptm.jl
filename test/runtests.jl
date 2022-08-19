using ClinicalTrialOptm
using Distributions
using Test

@testset "ClinicalTrialOptm.jl" begin
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
    @info "pgf of a center"
    @test pgf(ctry, 1.0) == 1.0
end
