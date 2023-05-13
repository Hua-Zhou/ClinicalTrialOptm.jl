# Inside make.jl
push!(LOAD_PATH,"../src/")
using ClinicalTrialOptm
using Documenter
makedocs(
         sitename = "ClinicalTrialOptm.jl",
         modules  = [ClinicalTrialOptm],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/Hua-Zhou/ClinicalTrialOptm.jl",
)
