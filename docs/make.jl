using Distributions, Documenter, HiGHS, JuMP, MathOptInterface,
Pajarito, Plots, SCIP, StatsPlots

makedocs(
    format = Documenter.HTML(),
    sitename = "ClinicalTrialOptm.jl",
    authors = "Arin Vansomphone, Hua Zhou",
    clean = true,
    debug = true,
    pages = [
        "index.md"
    ]
)

deploydocs(
    repo = "github.com/Hua-Zhou/ClinicalTrialOptm.jl.git"
    target = "build",
    deps = nothing,
    make = nothing
)