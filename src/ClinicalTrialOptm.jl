module ClinicalTrialOptm

using Cbc, Distributions, FFTW, HiGHS, QuadGK, Ipopt, JuMP, Juniper,Parajito
    KNITRO, SCIP, PrettyTables, SCS
import Statistics: mean, var
import Distributions: ccdf, cdf
export ccdf, cdf, ClinicalTrial, Country, mean, mean_cost, pgf, pmf, var, optdes!

include("country.jl")
include("clinicaltrial.jl")
include("optim.jl")

end
