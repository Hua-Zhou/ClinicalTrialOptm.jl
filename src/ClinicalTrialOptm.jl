module ClinicalTrialOptm

using Distributions, QuadGK
export Country, pgf

struct Country{T <: Real}
m  :: T # mean of Gamma-distributed enrollment rate 
s² :: T # var of Gamma-distributed enrollment rate 
l  :: Integer # lower bound of number of centers
u  :: Integer # upper bound of number of centers
c₀ :: T # cost of initializing one center
c  :: T # cost of running one center per unit of time
q  :: T # cost per one enrolled patient
d  :: T # probability of an enrolled patient dropping out
T₀ :: Distribution{Univariate, Continuous} # distribution of center initialization time
Td :: T # total duration of clinical trial
end

# constructor
function Country(
    m  :: T, # mean of Gamma-distributed enrollment rate 
    s² :: T, # var of Gamma-distributed enrollment rate 
    l  :: Integer, # lower bound of number of centers
    u  :: Integer, # upper bound of number of centers
    c₀ :: T, # cost of initializing one center
    c  :: T, # cost of running one center per unit of time
    q  :: T, # cost per one enrolled patient
    d  :: T, # probability of an enrolled patient dropping out
    T₀ :: Distribution{Univariate, Continuous}, # distribution of center initialization time    
    Td :: T
) where T <: Real
    Country{T}(m, s², l, u, c₀, c, q, d, T₀, Td)
end

function pgf(
    ctry :: Country{<:Real},
    z :: Union{Real, Complex}
    )
    θ = ctry.s² / ctry.m
    d = ctry.d
    α = ctry.m^2 / ctry.s²
    T₀ = ctry.T₀
    Td = ctry.Td
    (I, E) = quadgk(t -> (1 + θ * d * (Td - t) * (1 - z))^(-α) * pdf(T₀, t), minimum(ctry.T₀), maximum(ctry.T₀))
    I
end

end
