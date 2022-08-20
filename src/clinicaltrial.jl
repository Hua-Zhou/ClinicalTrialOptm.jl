struct ClinicalTrial{T <: Real, TI <: Integer}
    countries :: Vector{Country{T}}
    centers   :: Vector{TI} # number of centors in each country
end

# constructor
function ClinicalTrial(
    countries :: Vector{Country{T}}, 
    centers = zeros(Int, length(countries))
    ) where T
    ClinicalTrial{T, Int}(countries, centers)
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
    centers = zeros(Int, length(countries))
    ) where {T, TD <: Distribution{Univariate, Continuous}}
    n = length(m)
    countries = Vector{Country{T}}(undef, n)
    for i in eachindex(countries)
        countries[i] = Country{T}(
            m[i], s²[i], l[i], u[i], c₀[i], c[i], q[i], d[i], T₀[i], Td
            )
    end
    ClinicalTrial{T, Int}(countries, centers)
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
