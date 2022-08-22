using ClinicalTrialOptm
using Distributions, ForwardDiff, LinearAlgebra, Random, Test, UnicodePlots

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
    @info "pmf of a center"
    ctry_pmf = pmf(ctry)
    pltrange = 0:200
    display(
        lineplot(
            pltrange, 
            ctry_pmf[pltrange .+ 1], 
            xlabel = "Number of patients", 
            ylabel = "Probability"
        )
    )
    println()        
    @info "empirical mean/var"
    # @show ctry_pmf
    @show length(ctry_pmf)
    @show sum(ctry_pmf)
    μ̂ = dot(0:(length(ctry_pmf)-1), ctry_pmf)
    σ̂² = dot(abs2.(0:(length(ctry_pmf)-1)), ctry_pmf) - abs2(μ̂)
    println("μ̂ = $μ̂")
    println("σ̂² = $σ̂²")
end

@testset "ClinicalTrial" begin
    rng = MersenneTwister(123)
    nc = 5 # number of countries
    m = rand(rng, Uniform(1.0, 2.0), nc)
    s² = rand(rng, Uniform(1.0, 2.0), nc)
    l = rand(rng, 1:3, nc)
    u = l .+ rand(rng, 10:30, nc)
    c₀ = rand(rng, 15000.0:1000.0:25000.0, nc)
    c = rand(rng, 4000.0:200.0:5000.0, nc)
    q = rand(rng, 1000.0:200.0:2000.0, nc)
    d = rand(rng, 0.1:0.02:0.2, nc)
    T₀ = fill(Uniform(0.0, 6.0), nc)
    Td = 15.0
    centers = rand(rng, 5:15, nc)
    ct = ClinicalTrial(m, s², l, u, c₀, c, q, d, T₀, Td, centers)
    show(ct)
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
    @info "pmf of a clinical trial"
    ct_pmf = pmf(ct)
    pltrange = 450:1300
    display(
        lineplot(
            pltrange, 
            ct_pmf[pltrange .+ 1], 
        xlabel = "Number of patients", 
        ylabel = "Probability"
        )
    )
    println()
    @info "empirical mean/var"
    # @show ct_pmf
    @show length(ct_pmf)
    @show sum(ct_pmf)
    μ̂ = dot(0:(length(ct_pmf)-1), ct_pmf)
    σ̂² = dot(abs2.(0:(length(ct_pmf)-1)), ct_pmf) - abs2(μ̂)
    println("μ̂ = $μ̂")
    println("σ̂² = $σ̂²")
    @info "optdes!"
    optdes!(ct, 500, ps = 0.49)
    show(ct)
end
