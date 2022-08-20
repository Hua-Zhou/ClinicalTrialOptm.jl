module ClinicalTrialOptm

using Distributions, FFTW, QuadGK
import Statistics: mean, var
export ClinicalTrial, Country, mean, pgf, pmf, var

include("country.jl")
include("clinicaltrial.jl")

end
