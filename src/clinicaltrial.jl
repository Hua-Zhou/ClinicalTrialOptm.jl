struct ClinicalTrial{T <: Real, TI <: Integer}
    countries :: Vector{Country{T}}
    centers   :: Vector{TI} # number of centers in each country
    ntarget   :: Vector{TI} # target enrollment
    isoptm   :: Vector{Bool} # solved status
    solution_status :: Vector{String} # solution status with solver

end

# constructor
function ClinicalTrial(
    countries :: Vector{Country{T}}, 
    centers = zeros(Int, length(countries)),
    ntarget = 0
    ) where T
    ClinicalTrial{T, Int}(countries, centers, [ntarget], [false],["Unsolved"])
end

""" 
    ClinicalTrial(
        m  :: AbstractVector{T},
        s² :: AbstractVector{T}, 
        l  :: AbstractVector{<:Integer}, 
        u  :: AbstractVector{<:Integer}, 
        c₀ :: AbstractVector{T}, 
        c  :: AbstractVector{T}, 
        q  :: AbstractVector{T}, 
        d  :: AbstractVector{T}, 
        T₀ :: AbstractVector{TD},    
        Td :: T,
        centers = zeros(Int, length(m)),
        ntarget = 0
    )

Store parameters for a clinical trial. 
    
# Arguments
- `m`: expected mean enrollment rates in country
- `s²`: variances of enrollment rates in country
- `l`: lower bound of centers in country
- `u`: upper bound of centers in country
- `c₀`: cost of initializing one center in country
- `c`: cost of running one center in country per unit of time
- `q`: cost of one enrolled patient in country
- `d`: probability of an enrolled patient dropping out in country
- `T₀`: distribution of center initialization time in country
- `Td`: duration of clinical trial 
- `centers`: number of centers in each country for the optimal trial design
- `ntarget`: target enrollment of patients

When inputting values, each parameter's index must corresponding with the other country values.
Be sure to verify your values with the displayed table before applying the solver.

See also: [`mean(ct :: ClinicalTrial)`](@ref), [`var(ct :: ClinicalTrial)`](@ref),
[`pgf(ct :: ClinicalTrial, z)`](@ref), [`pmf(ct :: ClinicalTrial)`](@ref),
[`cdf`](@ref), [`ccdf`](@ref),
[`optdes!`](@ref)
"""
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
    ClinicalTrial{T, Int}(countries, centers, [ntarget], [false], ["Unsolved"])
end

""" 
    mean(ct :: ClinicalTrial)

Compute the mean number of patients enrolled in a clinical trial `ct`.

See also [`var`](@ref), [`pgf`](@ref)
"""
mean(ct :: ClinicalTrial) = 
    mapreduce((c, x) -> mean(c) * x, +, ct.countries, ct.centers)

""" 
    var(ct :: ClinicalTrial)

Compute the variance of the number of patients enrolled a clinical trial `ct`.

See also [`mean`](@ref), [`pgf`](@ref)
"""
var(ct :: ClinicalTrial) = 
    mapreduce((c, x) -> var(c) * x, +, ct.countries, ct.centers)

"""
    pgf(ct :: ClinicalTrial, z)

Compute the probability generating function of the number of patients enrolled in a
clinical trial `ct`.

See also [`pmf`](@ref)
"""
pgf(ct :: ClinicalTrial, z) = 
    mapreduce((c, x) -> pgf(c, z)^x, *, ct.countries, ct.centers)

"""
    pmf(ct :: ClinicalTrial)

Compute the probability mass function of the number of patients enrolled in a
clinical trial `ct`, by a Fast Fourier Transform of its pgf.

See also [`pgf`](@ref)
"""
function pmf(ct :: ClinicalTrial)
    # support of pmf is [0, n - 1]
    μ, σ² = mean(ct), var(ct)
    n = 2^ceil(log2(μ + 10sqrt(σ²)))
    # fft
    pmf = max.(real.(fft!([pgf(ct, exp(2π * im * j / n)) for j in 0:(n-1)])), 0)
    pmf ./= sum(pmf)
end

"""
    cdf(ct :: ClinicalTrial, x :: Real)

Compute the cumulative distribution function of the number of patients enrolled in a
clinical trial `ct` evaluated at x, P(X ≤ x), by applying the Gil-Pelaez inversion formula 
(<https://en.wikipedia.org/wiki/Characteristic_function_(probability_theory)#Inversion_formula>)
to the characteristic function. 

See also [`ccdf`](@ref)
"""
function cdf(ct :: ClinicalTrial, x :: Real)
    z = round(Int, x) == x ? x + 0.5 : Float64(x)
    I, _ = quadgk(t -> imag(exp(- im * t * z) * pgf(ct, exp(im * t))) / t, 0, Inf)
    0.5 - inv(π) * I
end

"""
    ccdf(ct :: ClinicalTrial, x :: Real) 

Compute the complementary distribution function of the number of patients enrolled in a 
clinical trial `ct` evaluated by 1 - `cdf(ct, x)`.

See also [`cdf`](@ref)
"""
ccdf(ct :: ClinicalTrial, x :: Real) = 1 - cdf(ct, x)

function Base.show(io::IO, ct::ClinicalTrial)
    println(io)
    println(io, "Global Clinical Trial:")
    println(io)
    if ct.isoptm[1]
        println(io, "Optimal center assignment calculated.")
        println(io, ct.solution_status[1])
        println(io, "Total duration (months): $(ct.countries[1].Td)")
        println(io, "Target enrollment: $(ct.ntarget[1])")
        μ, σ² = mean(ct), var(ct)
        println(io, "Probability of success (based on normal approximation): $(ccdf(Normal(μ, sqrt(σ²)), ct.ntarget[1]))")
        println(io, "Probability of success (based on Poisson-Gamma model): $(ccdf(ct, ct.ntarget[1]))")
        if ccdf(ct, ct.ntarget[1]) < ccdf(Normal(μ, sqrt(σ²)), ct.ntarget[1])
            println(io, "WARNING: Probability of success used in optimization is lower than actual. Consider adjusting it.")
        end
        μcost = mapreduce((c, x) -> mean_cost(c) * x, +, ct.countries, ct.centers)
        println(io, "Expected cost (\$): $(μcost)")
    else
        println(io, "Optimal center assignment not computed.")
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
        ),
        crop = :none
    )
    println(io)
    println(io)
end
