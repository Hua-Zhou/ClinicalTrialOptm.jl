module ClinicalTrialOptm

using Distributions, QuadGK
import Statistics: mean, var
export Country, mean, pgf, var

include("country.jl")

end
