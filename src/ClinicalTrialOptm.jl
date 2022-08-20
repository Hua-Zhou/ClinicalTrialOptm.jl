module ClinicalTrialOptm

using Distributions, QuadGK
import Statistics: mean, var
export ClinicalTrial, Country, mean, pgf, var

include("country.jl")
include("clinicaltrial.jl")

end
