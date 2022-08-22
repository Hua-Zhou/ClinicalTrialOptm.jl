module ClinicalTrialOptm

using Distributions, FFTW, 
    HiGHS, Hypatia, QuadGK, JuMP, Pajarito, PrettyTables
import Statistics: mean, var
export ClinicalTrial, Country, mean, mean_cost, pgf, pmf, var, optdes!

include("country.jl")
include("clinicaltrial.jl")
include("optim.jl")

end
