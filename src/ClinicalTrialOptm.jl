module ClinicalTrialOptm

using Distributions, FFTW, HiGHS, QuadGK, Ipopt, JuMP, Juniper, 
    Pajarito, PrettyTables, SCS
import Statistics: mean, var
export ClinicalTrial, Country, mean, mean_cost, pgf, pmf, var, optdes!

include("country.jl")
include("clinicaltrial.jl")
include("optim.jl")

end
