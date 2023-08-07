# Inside make.jl
using Documenter, ClinicalTrialOptm

makedocs(
    format = Documenter.HTML(),
    sitename = "ClinicalTrialOptm.jl",
    authors = "Hua Zhou, Arin Vansomphone",
    clean = true,
    debug = true,
    pages = [
        "index.md"
    ]    
)

deploydocs(
    repo   = "github.com/Hua-Zhou/ClinicalTrialOptm.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
