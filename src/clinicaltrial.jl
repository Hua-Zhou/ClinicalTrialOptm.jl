struct ClinicalTrial{T <: Real, TI <: Integer}
    countries :: Vector{Country{T}}
    centers   :: Vector{TI} # number of centors in each country
    ntarget   :: Vector{TI} # target enrollment
    isoptm    :: Vector{Bool} 
end

# constructor
function ClinicalTrial(
    countries :: Vector{Country{T}}, 
    centers = zeros(Int, length(countries)),
    ntarget = 0
    ) where T
    ClinicalTrial{T, Int}(countries, centers, [ntarget], [false])
end

function ClinicalTrial(
    m  :: AbstractVector{T}, # mean of Gamma-distributed enrollment rate 
    s² :: AbstractVector{T}, # var of Gamma-distributed enrollment rate 
    l  :: AbstractVector{<:Integer}, # lower bound of number of centers
    u  :: AbstractVector{<:Integer}, # upper bound of number of centers
    c₀ :: AbstractVector{T}, # cost of initializing one center
    c  :: AbstractVector{T}, # cost of running one center per unit of time
    q  :: AbstractVector{T}, # cost per one enrolled patient
    d  :: AbstractVector{T}, # probability of an enrolled patient dropping out
    T₀ :: AbstractVector{TD}, # distribution of center initialization time    
    Td :: T,
    centers = zeros(Int, length(m)),
    ntarget = 0
    ) where {T, TD <: Distribution{Univariate, Continuous}}
    n = length(m)
    countries = Vector{Country{T}}(undef, n)
    for i in eachindex(countries)
        countries[i] = Country{T}(
            m[i], s²[i], l[i], u[i], c₀[i], c[i], q[i], d[i], T₀[i], Td
            )
    end
    ClinicalTrial{T, Int}(countries, centers, [ntarget], [false])
end

""" 
    mean(ct :: ClinicalTrial)

Mean number of patients enrolled in a clinical trial `ct`.
"""
mean(ct :: ClinicalTrial) = 
    mapreduce((c, x) -> mean(c) * x, +, ct.countries, ct.centers)

""" 
    var(ct :: ClinicalTrial)

Variance of the number of patients enrolled a clinical trial `ct`.
"""
var(ct :: ClinicalTrial) = 
    mapreduce((c, x) -> var(c) * x, +, ct.countries, ct.centers)

"""
    pgf(ct :: ClinicalTrial, z)

Probability generating function of the number of patients enrolled in a
clinical trial `ct`.
"""
pgf(ct :: ClinicalTrial, z) = 
    mapreduce((c, x) -> pgf(c, z)^x, *, ct.countries, ct.centers)

"""
    pmf(ct :: ClinicalTrial)

Probability mass function of the number of patients enrolled in a
clinical trial `ct`, by FFT of its pgf.
"""
function pmf(ct :: ClinicalTrial)
    # support of pmf is [0, n - 1]
    μ, σ² = mean(ct), var(ct)
    n = 2^ceil(log2(μ + 10sqrt(σ²)))
    # fft
    pmf = max.(real.(fft!([pgf(ct, exp(2π * im * j / n)) for j in 0:(n-1)])), 0)
    pmf ./= sum(pmf)
end

function Base.show(io::IO, ct::ClinicalTrial)
    println(io)
    println(io, "Global Clinical Trial:")
    println(io)
    if ct.isoptm[1]
        println(io, "Optimal center assignment successfully found")
        println(io, "Total duration (months): $(ct.countries[1].Td)")
        println(io, "Target enrollment: $(ct.ntarget[1])")
        μ, σ² = mean(ct), var(ct)
        println(io, "Probability of success: $(ccdf(Normal(μ, sqrt(σ²)), ct.ntarget[1]))")
        μcost = mapreduce((c, x) -> mean_cost(c) * x, +, ct.countries, ct.centers)
        println(io, "Expected cost (\$): $(μcost)")
    else
        println(io, "Optimal center assignment not computed")
    end
    J = length(ct.countries) # number of countries
    tbl = hcat(
        1:J,
        [ct.countries[j].m for j in 1:J], 
        [ct.countries[j].s² for j in 1:J],
        [ct.countries[j].c₀ for j in 1:J],
        [ct.countries[j].c for j in 1:J],
        [ct.countries[j].q for j in 1:J],
        [ct.countries[j].d for j in 1:J],
        [ct.countries[j].l for j in 1:J],
        [ct.countries[j].u for j in 1:J],
        ct.isoptm[1] ? [ct.centers[j] for j in 1:J] : fill("NA", J)
    )
    pretty_table(
        io, 
        tbl,
        header = (
            ["Country", "mean(λ)", "var(λ)", 
            "init. cost", "maint. cost", "enroll. cost",
            "drop out rate", "min. centers", "max. centers", "centers"],
            ["", "", "", "\$/center", "\$/center/month", "\$/patient", "", "", "", ""]
        ),
        formatters = (
            ft_printf("%5.2f", [2, 3, 7]), 
            ft_printf("%8.0f", [4, 5, 6]),
            ft_printf("%i", [1, 8, 9, 10])
        )
    )
    println(io)
    println(io)
end
