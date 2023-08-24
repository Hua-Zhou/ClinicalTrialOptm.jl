var documenterSearchIndex = {"docs":
[{"location":"#ClinicalTrialOptm.jl","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"","category":"section"},{"location":"#The-Problem","page":"ClinicalTrialOptm.jl","title":"The Problem","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"(Image: clinicaltrialimage.png)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"In the last five years, clinical trials have become increasingly more difficult to conduct due to staffing, budget, and protocol complications. According to the 2020 Tufts University Impact Report, more than 20% of clinical trials fail to recruit enough patients in time. Biomedical and phameceutical companies rely on trial data to progress in treatment development, making effective clinical trial design necessary to address.","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"ClinicalTrialOptm.jl solves these multi-center, multi-state recruitment problems using mixed-integer algorithms. It seeks to optimize the number of clinics for each country, minimizing cost while maintaining high probabilities of successful recruitment. Details on the calculations are described in the paper: (Insert paper here)","category":"page"},{"location":"#Installation","page":"ClinicalTrialOptm.jl","title":"Installation","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"ClinicalTrialOptm.jl requires Julia v1.7 or later. The package has not been registered yet and must be isntalled using the repository location. To do so, start Julia and use the ] key to switch to the package manager REPL:","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"(@v1.8) Pkg> add https://github.com/Hua-Zhou/ClinicalTrialOptm.jl.git","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Use the backspace key to return to the Julia REPL.","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"versioninfo()","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Julia Version 1.8.1\nCommit afb6c60d69a (2022-09-06 15:09 UTC)\nPlatform Info:\n  OS: macOS (arm64-apple-darwin21.5.0)\n  CPU: 8 × Apple M1\n  WORD_SIZE: 64\n  LIBM: libopenlibm\n  LLVM: libLLVM-13.0.1 (ORCJIT, apple-m1)\n  Threads: 4 on 4 virtual cores\nEnvironment:\n  JULIA_NUM_THREADS = 4","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"# for use in this tutorial\nusing ClinicalTrialOptm\nusing Distributions, HiGHS, SCS, Pajarito, SCIP, MathOptInterface, JuMP, Plots, StatsPlots\nconst MOI = MathOptInterface;","category":"page"},{"location":"#Inputting-Parameters","page":"ClinicalTrialOptm.jl","title":"Inputting Parameters","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"There are two types of data structs within this package: Country and ClinicalTrial:","category":"page"},{"location":"#Country","page":"ClinicalTrialOptm.jl","title":"Country","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Country","category":"page"},{"location":"#ClinicalTrialOptm.Country","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.Country","text":"Country(\n    m  :: T, \n    s² :: T, \n    l  :: Integer, \n    u  :: Integer,\n    c₀ :: T, \n    c  :: T, \n    q  :: T, \n    d  :: T, \n    T₀ :: Distribution{Univariate, Continuous},\n    Td :: T \n)\n\nStores parameters for a country's centers.\n\nArguments\n\nm: expected mean enrollment rates in country\ns²: variances of enrollment rates in country\nl: lower bound of centers in country\nu: upper bound of centers in country\nc₀: cost of initializing one center in country\nc: cost of running one center in country per unit of time\nq: cost of one enrolled patient in country\nd: probability of an enrolled patient dropping out in country\nT₀: distribution of center initialization time in country\nTd: duration of clinical trial \n\nSee also mean, var, mean_cost, pgf, pmf\n\n\n\n\n\n","category":"type"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Example use:","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"m = 1.5\ns² = 2.0\nl, u = 0, 10\nc₀ = 20_000.0\nc = 5_000.0\nq = 1_500.0\nd = 0.05\nT₀ = Uniform(0.0, 6.0)\nTd = 24.0\nctry = Country(m, s², l, u, c₀, c, q, d, T₀, Td)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Country{Float64}(1.5, 2.0, 0, 10, 20000.0, 5000.0, 1500.0, 0.05, Uniform{Float64}(a=0.0, b=6.0), 24.0)","category":"page"},{"location":"#ClinicalTrial","page":"ClinicalTrialOptm.jl","title":"ClinicalTrial","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"ClinicalTrial","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Example use:","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"m = [1, 1.2, 1.4, 1.2]\ns² = [0.2, 0.4, 0.8, 0.6]\nl = [0, 3, 2, 1]\nu = [8, 5, 7, 4]\nc₀ = [15000.0, 13000.0, 16000.0, 17000.0]\nc = [3000.0, 2000.0, 5000.0, 8000.0]\nq = [1000.0, 1300.0, 900.0, 800.0]\nd = [0.01, 0.05, 0.09, 0.15]\nT₀ = fill(Uniform(0.0, 6.0), 4)\nTd = 24.0\nct = ClinicalTrial(m, s², l, u, c₀, c, q, d, T₀, Td)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Global Clinical Trial:\n\nOptimal center assignment not computed.\n┌─────────┬─────────┬────────┬────────────┬────────────────┬──────────────┬───────────────┬──────────────┬──────────────┬─────────┐\n│\u001b[1m Country \u001b[0m│\u001b[1m mean(λ) \u001b[0m│\u001b[1m var(λ) \u001b[0m│\u001b[1m init. cost \u001b[0m│\u001b[1m    maint. cost \u001b[0m│\u001b[1m enroll. cost \u001b[0m│\u001b[1m drop out rate \u001b[0m│\u001b[1m min. centers \u001b[0m│\u001b[1m max. centers \u001b[0m│\u001b[1m centers \u001b[0m│\n│\u001b[90m         \u001b[0m│\u001b[90m         \u001b[0m│\u001b[90m        \u001b[0m│\u001b[90m   $/center \u001b[0m│\u001b[90m $/center/month \u001b[0m│\u001b[90m    $/patient \u001b[0m│\u001b[90m               \u001b[0m│\u001b[90m              \u001b[0m│\u001b[90m              \u001b[0m│\u001b[90m         \u001b[0m│\n├─────────┼─────────┼────────┼────────────┼────────────────┼──────────────┼───────────────┼──────────────┼──────────────┼─────────┤\n│       1 │    1.00 │   0.20 │      15000 │           3000 │         1000 │          0.01 │            0 │            8 │      NA │\n│       2 │    1.20 │   0.40 │      13000 │           2000 │         1300 │          0.05 │            3 │            5 │      NA │\n│       3 │    1.40 │   0.80 │      16000 │           5000 │          900 │          0.09 │            2 │            7 │      NA │\n│       4 │    1.20 │   0.60 │      17000 │           8000 │          800 │          0.15 │            1 │            4 │      NA │\n└─────────┴─────────┴────────┴────────────┴────────────────┴──────────────┴───────────────┴──────────────┴──────────────┴─────────┘","category":"page"},{"location":"#Using-optdes!","page":"ClinicalTrialOptm.jl","title":"Using optdes!","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"optdes!","category":"page"},{"location":"#ClinicalTrialOptm.optdes!","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.optdes!","text":"optdes!(ct :: ClinicalTrial, ntarget :: Integer, ps :: Real, solver)\n\nFind the optimal clinical trial design, with ntarget enrollment and guaranteed  probability of success ps.\n\nThe solver can be changed to the user's preferance by setting solver equal to your solver choice (see documentation page for examples).  The function overwrites ct.centers by the optimal solution.\n\nSee also lbtest, ubtest, ClinicalTrial\n\n\n\n\n\n","category":"function"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Example use:","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"# input parameters\nm = [1, 1.2, 1.4, 1.2]\ns² = [0.2, 0.4, 0.8, 0.6]\nl = [0, 3, 2, 1]\nu = [8, 5, 7, 4]\nc₀ = [15000.0, 13000.0, 16000.0, 17000.0]\nc = [3000.0, 2000.0, 5000.0, 8000.0]\nq = [1000.0, 1300.0, 900.0, 800.0]\nd = [0.01, 0.05, 0.09, 0.15]\nT₀ = fill(Uniform(0.0, 6.0), 4)\nTd = 24.0\nct = ClinicalTrial(m, s², l, u, c₀, c, q, d, T₀, Td)\n\noptdes!(ct, 400, ps = 0.85)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Lower bound probability of success: 1.1138624195141726e-11\nThe optimal solution is not the lower bound of the centers.\nUpper bound probability of success: 0.9878987929122786\nThe optimal solution is feasible.\n  2.427250 seconds (6.64 M allocations: 356.902 MiB, 1.82% gc time, 99.28% compilation time: 25% of which was recompilation)\nsolution_summary(model) = * Solver : SCIP\n\n* Status\n  Result count       : 4\n  Termination status : OPTIMAL\n  Message from the solver:\n  \"SCIP_STATUS_OPTIMAL\"\n\n* Candidate solution (result #1)\n  Primal status      : FEASIBLE_POINT\n  Dual status        : NO_SOLUTION\n  Objective value    : 2.29354e+06\n  Objective bound    : 2.29354e+06\n  Relative gap       : 0.00000e+00\n\n* Work counters\n  Solve time (sec)   : 8.66800e-03\n  Simplex iterations : 8\n  Node count         : 1\n\ntermination_status(model) = MathOptInterface.OPTIMAL\nprimal_status(model) = MathOptInterface.FEASIBLE_POINT\nobjective_value(model) = 2.293537599999999e6\n\n\n\n\n\n\nGlobal Clinical Trial:\n\nOptimal center assignment calculated.\nAn optimal solution has been found.\nTotal duration (months): 24.0\nTarget enrollment: 400\nProbability of success (based on normal approximation): 0.8587747637053408\nProbability of success (based on Poisson-Gamma model): 0.8598777881463103\nExpected cost ($): 2.2935376e6\n┌─────────┬─────────┬────────┬────────────┬────────────────┬──────────────┬───────────────┬──────────────┬──────────────┬─────────┐\n│\u001b[1m Country \u001b[0m│\u001b[1m mean(λ) \u001b[0m│\u001b[1m var(λ) \u001b[0m│\u001b[1m init. cost \u001b[0m│\u001b[1m    maint. cost \u001b[0m│\u001b[1m enroll. cost \u001b[0m│\u001b[1m drop out rate \u001b[0m│\u001b[1m min. centers \u001b[0m│\u001b[1m max. centers \u001b[0m│\u001b[1m centers \u001b[0m│\n│\u001b[90m         \u001b[0m│\u001b[90m         \u001b[0m│\u001b[90m        \u001b[0m│\u001b[90m   $/center \u001b[0m│\u001b[90m $/center/month \u001b[0m│\u001b[90m    $/patient \u001b[0m│\u001b[90m               \u001b[0m│\u001b[90m              \u001b[0m│\u001b[90m              \u001b[0m│\u001b[90m         \u001b[0m│\n├─────────┼─────────┼────────┼────────────┼────────────────┼──────────────┼───────────────┼──────────────┼──────────────┼─────────┤\n│       1 │    1.00 │   0.20 │      15000 │           3000 │         1000 │          0.01 │            0 │            8 │       8 │\n│       2 │    1.20 │   0.40 │      13000 │           2000 │         1300 │          0.05 │            3 │            5 │       5 │\n│       3 │    1.40 │   0.80 │      16000 │           5000 │          900 │          0.09 │            2 │            7 │       6 │\n│       4 │    1.20 │   0.60 │      17000 │           8000 │          800 │          0.15 │            1 │            4 │       1 │\n└─────────┴─────────┴────────┴────────────┴────────────────┴──────────────┴───────────────┴──────────────┴──────────────┴─────────┘","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"optdes! uses the normal approximation of the probability of success (PoS) to compute solutions, however the Poisson-Gamma model probability is displayed because it is the more accurate estimation of the PoS. If the Poisson-Gamma model PoS is less than the normal approximation PoS, then a warning will be outputted to let users know.","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Example:","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"optdes!(ct, 300, ps = 0.73)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Lower bound probability of success: 2.57712065043568e-5\nThe optimal solution is not the lower bound of the centers.\nUpper bound probability of success: 0.999879502628558\nThe optimal solution is feasible.\n  0.013808 seconds (1.48 k allocations: 51.547 KiB)\nsolution_summary(model) = * Solver : SCIP\n\n* Status\n  Result count       : 3\n  Termination status : OPTIMAL\n  Message from the solver:\n  \"SCIP_STATUS_GAPLIMIT\"\n\n* Candidate solution (result #1)\n  Primal status      : FEASIBLE_POINT\n  Dual status        : NO_SOLUTION\n  Objective value    : 1.61443e+06\n  Objective bound    : 1.61289e+06\n  Relative gap       : 9.59704e-04\n\n* Work counters\n  Solve time (sec)   : 7.47200e-03\n  Simplex iterations : 18\n  Node count         : 1\n\ntermination_status(model) = MathOptInterface.OPTIMAL\nprimal_status(model) = MathOptInterface.FEASIBLE_POINT\nobjective_value(model) = 1.6144332000000002e6\n\n\n\n\n\n\nGlobal Clinical Trial:\n\nOptimal center assignment calculated.\nAn optimal solution has been found.\nTotal duration (months): 24.0\nTarget enrollment: 300\nProbability of success (based on normal approximation): 0.7852131181816588\nProbability of success (based on Poisson-Gamma model): 0.7772523002943241\nWARNING: Probability of success used in optimization is lower than actual. Consider adjusting it.\nExpected cost ($): 1.6144332e6\n┌─────────┬─────────┬────────┬────────────┬────────────────┬──────────────┬───────────────┬──────────────┬──────────────┬─────────┐\n│\u001b[1m Country \u001b[0m│\u001b[1m mean(λ) \u001b[0m│\u001b[1m var(λ) \u001b[0m│\u001b[1m init. cost \u001b[0m│\u001b[1m    maint. cost \u001b[0m│\u001b[1m enroll. cost \u001b[0m│\u001b[1m drop out rate \u001b[0m│\u001b[1m min. centers \u001b[0m│\u001b[1m max. centers \u001b[0m│\u001b[1m centers \u001b[0m│\n│\u001b[90m         \u001b[0m│\u001b[90m         \u001b[0m│\u001b[90m        \u001b[0m│\u001b[90m   $/center \u001b[0m│\u001b[90m $/center/month \u001b[0m│\u001b[90m    $/patient \u001b[0m│\u001b[90m               \u001b[0m│\u001b[90m              \u001b[0m│\u001b[90m              \u001b[0m│\u001b[90m         \u001b[0m│\n├─────────┼─────────┼────────┼────────────┼────────────────┼──────────────┼───────────────┼──────────────┼──────────────┼─────────┤\n│       1 │    1.00 │   0.20 │      15000 │           3000 │         1000 │          0.01 │            0 │            8 │       7 │\n│       2 │    1.20 │   0.40 │      13000 │           2000 │         1300 │          0.05 │            3 │            5 │       5 │\n│       3 │    1.40 │   0.80 │      16000 │           5000 │          900 │          0.09 │            2 │            7 │       2 │\n│       4 │    1.20 │   0.60 │      17000 │           8000 │          800 │          0.15 │            1 │            4 │       1 │\n└─────────┴─────────┴────────┴────────────┴────────────────┴──────────────┴───────────────┴──────────────┴──────────────┴─────────┘","category":"page"},{"location":"#Solver-Options","page":"ClinicalTrialOptm.jl","title":"Solver Options","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"optdes! uses Pajarito.jl to solve convex cases (when ps < 0.5) and switches to SCIP.jl for non-convex cases (when ps ≥ 0.5). The exact specifications are listed below.","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"# Pajarito solver\noa_solver = optimizer_with_attributes(HiGHS.Optimizer,\nMOI.Silent() => true,\n\"mip_feasibility_tolerance\" => 1e-8,\n\"mip_rel_gap\" => 1e-6,\n)\nconic_solver = optimizer_with_attributes(SCS.Optimizer, \nMOI.Silent() => true,\n)\nconvex_default_solver = optimizer_with_attributes(Pajarito.Optimizer,\n    \"time_limit\" => 600, \n    \"oa_solver\" => oa_solver, \n    \"conic_solver\" => conic_solver,\n    \"tol_rel_gap\" => 0.001\n)\n\n# SCIP solver\nnonconvex_default_solver =  optimizer_with_attributes(\n    SCIP.Optimizer, \"display/verblevel\"=>0, \"limits/gap\"=>0.001 \n);","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"As stated in the documentation, users can change the solver to their preference by specifying it in the solver argument. Here is an example using the KNITRO solver.","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"knitro_solver = optimizer_with_attributes(\n        KNITRO.Optimizer,\n        \"mip_opt_gap_rel\" => 0.001\n    )\n    \noptdes!(ct, 300, ps = 0.8, solver = knitro_solver)","category":"page"},{"location":"#Other-Functions","page":"ClinicalTrialOptm.jl","title":"Other Functions","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"ClinicalTrialOptm.jl contains numerous other functions for patient amounts, cost, and probabilities of countries and clinical trials. ","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"# for example uses\nm = 1.5\ns² = 2.0\nl, u = 0, 10\nc₀ = 20_000.0\nc = 5_000.0\nq = 1_500.0\nd = 0.05\nT₀ = Uniform(0.0, 6.0)\nTd = 24.0\nctry = Country(m, s², l, u, c₀, c, q, d, T₀, Td)\n\nm = [1, 1.2, 1.4, 1.2]\ns² = [0.2, 0.4, 0.8, 0.6]\nl = [0, 3, 2, 1]\nu = [8, 5, 7, 4]\nc₀ = [15000.0, 13000.0, 16000.0, 17000.0]\nc = [3000.0, 2000.0, 5000.0, 8000.0]\nq = [1000.0, 1300.0, 900.0, 800.0]\nd = [0.01, 0.05, 0.09, 0.15]\nT₀ = fill(Uniform(0.0, 6.0), 4)\nTd = 24.0\ncenters = [2, 4, 7, 3]\nct = ClinicalTrial(m, s², l, u, c₀, c, q, d, T₀, Td, centers);","category":"page"},{"location":"#mean","page":"ClinicalTrialOptm.jl","title":"mean","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"ClinicalTrialOptm.mean","category":"page"},{"location":"#Statistics.mean","page":"ClinicalTrialOptm.jl","title":"Statistics.mean","text":"mean(itr)\n\nCompute the mean of all elements in a collection.\n\nnote: Note\nIf itr contains NaN or missing values, the result is also NaN or missing (missing takes precedence if array contains both). Use the skipmissing function to omit missing entries and compute the mean of non-missing values.\n\nExamples\n\njulia> using Statistics\n\njulia> mean(1:20)\n10.5\n\njulia> mean([1, missing, 3])\nmissing\n\njulia> mean(skipmissing([1, missing, 3]))\n2.0\n\n\n\n\n\nmean(f, itr)\n\nApply the function f to each element of collection itr and take the mean.\n\njulia> using Statistics\n\njulia> mean(√, [1, 2, 3])\n1.3820881233139908\n\njulia> mean([√1, √2, √3])\n1.3820881233139908\n\n\n\n\n\nmean(f, A::AbstractArray; dims)\n\nApply the function f to each element of array A and take the mean over dimensions dims.\n\ncompat: Julia 1.3\nThis method requires at least Julia 1.3.\n\njulia> using Statistics\n\njulia> mean(√, [1, 2, 3])\n1.3820881233139908\n\njulia> mean([√1, √2, √3])\n1.3820881233139908\n\njulia> mean(√, [1 2 3; 4 5 6], dims=2)\n2×1 Matrix{Float64}:\n 1.3820881233139908\n 2.2285192400943226\n\n\n\n\n\nmean(A::AbstractArray; dims)\n\nCompute the mean of an array over the given dimensions.\n\ncompat: Julia 1.1\nmean for empty arrays requires at least Julia 1.1.\n\nExamples\n\njulia> using Statistics\n\njulia> A = [1 2; 3 4]\n2×2 Matrix{Int64}:\n 1  2\n 3  4\n\njulia> mean(A, dims=1)\n1×2 Matrix{Float64}:\n 2.0  3.0\n\njulia> mean(A, dims=2)\n2×1 Matrix{Float64}:\n 1.5\n 3.5\n\n\n\n\n\nmean(A::AbstractArray, w::AbstractWeights[, dims::Int])\n\nCompute the weighted mean of array A with weight vector w (of type AbstractWeights). If dim is provided, compute the weighted mean along dimension dims.\n\nExamples\n\nn = 20\nx = rand(n)\nw = rand(n)\nmean(x, weights(w))\n\n\n\n\n\nmean(d::UnivariateDistribution)\n\nCompute the expectation.\n\n\n\n\n\nmean(d::MultivariateDistribution)\n\nCompute the mean vector of distribution d.\n\n\n\n\n\nmean(d::MatrixDistribution)\n\nReturn the mean matrix of d.\n\n\n\n\n\nmean(d::Union{UnivariateMixture, MultivariateMixture})\n\nCompute the overall mean (expectation).\n\n\n\n\n\nmean(ctry :: Country)\n\nCalculate the mean number of patients enrolled by a center in country ctry.\n\nSee also var, pgf\n\n\n\n\n\nmean(ct :: ClinicalTrial)\n\nCompute the mean number of patients enrolled in a clinical trial ct.\n\nSee also var, pgf\n\n\n\n\n\n","category":"function"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Example uses:","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"mean(ctry)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"29.924999999999997","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"mean(ct)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"388.878","category":"page"},{"location":"#var","page":"ClinicalTrialOptm.jl","title":"var","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"ClinicalTrialOptm.var","category":"page"},{"location":"#Statistics.var","page":"ClinicalTrialOptm.jl","title":"Statistics.var","text":"var(itr; corrected::Bool=true, mean=nothing[, dims])\n\nCompute the sample variance of collection itr.\n\nThe algorithm returns an estimator of the generative distribution's variance under the assumption that each entry of itr is a sample drawn from the same unknown distribution, with the samples uncorrelated. For arrays, this computation is equivalent to calculating sum((itr .- mean(itr)).^2) / (length(itr) - 1). If corrected is true, then the sum is scaled with n-1, whereas the sum is scaled with n if corrected is false where n is the number of elements in itr.\n\nIf itr is an AbstractArray, dims can be provided to compute the variance over dimensions.\n\nA pre-computed mean may be provided. When dims is specified, mean must be an array with the same shape as mean(itr, dims=dims) (additional trailing singleton dimensions are allowed).\n\nnote: Note\nIf array contains NaN or missing values, the result is also NaN or missing (missing takes precedence if array contains both). Use the skipmissing function to omit missing entries and compute the variance of non-missing values.\n\n\n\n\n\nvar(x::AbstractArray, w::AbstractWeights, [dim]; mean=nothing, corrected=false)\n\nCompute the variance of a real-valued array x, optionally over a dimension dim. Observations in x are weighted using weight vector w. The uncorrected (when corrected=false) sample variance is defined as:\n\nfrac1sumw sum_i=1^n w_ileft(x_i - μright)^2 \n\nwhere n is the length of the input and μ is the mean. The unbiased estimate (when corrected=true) of the population variance is computed by replacing frac1sumw with a factor dependent on the type of weights used:\n\nAnalyticWeights: frac1sum w - sum w^2  sum w\nFrequencyWeights: frac1sumw - 1\nProbabilityWeights: fracn(n - 1) sum w where n equals count(!iszero, w)\nWeights: ArgumentError (bias correction not supported)\n\n\n\n\n\nvar(ce::CovarianceEstimator, x::AbstractVector; mean=nothing)\n\nCompute the variance of the vector x using the estimator ce.\n\n\n\n\n\nvar(d::UnivariateDistribution)\n\nCompute the variance. (A generic std is provided as std(d) = sqrt(var(d)))\n\n\n\n\n\nvar(d::MultivariateDistribution)\n\nCompute the vector of element-wise variances for distribution d.\n\n\n\n\n\nvar(d::MatrixDistribution)\n\nCompute the matrix of element-wise variances for distribution d.\n\n\n\n\n\nvar(d::UnivariateMixture)\n\nCompute the overall variance (only for UnivariateMixture).\n\n\n\n\n\nvar(ctry :: Country)\n\nCalculate the variance of the number of patients enrolled by a center in country ctry.\n\nSee also mean, pgf\n\n\n\n\n\nvar(ct :: ClinicalTrial)\n\nCompute the variance of the number of patients enrolled a clinical trial ct.\n\nSee also mean, pgf\n\n\n\n\n\n","category":"function"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Example uses:","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"var(ctry)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"837.436875","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"var(ct)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"3905.413596","category":"page"},{"location":"#mean_cost","page":"ClinicalTrialOptm.jl","title":"mean_cost","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"mean_cost","category":"page"},{"location":"#ClinicalTrialOptm.mean_cost","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.mean_cost","text":"mean_cost(ctry :: Country)\n\nCalculate the mean cost of a center in country ctry.\n\n\n\n\n\n","category":"function"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Example use:","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"mean_cost(ctry)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"169887.5","category":"page"},{"location":"#pgf","page":"ClinicalTrialOptm.jl","title":"pgf","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"pgf","category":"page"},{"location":"#ClinicalTrialOptm.pgf","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.pgf","text":"pgf(ctry :: Country, z)\n\nCalculate the probability generating function of the number of patients enrolled by one  center in country ctry.\n\nSee also pmf\n\n\n\n\n\npgf(ct :: ClinicalTrial, z)\n\nCompute the probability generating function of the number of patients enrolled in a clinical trial ct.\n\nSee also pmf\n\n\n\n\n\n","category":"function"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Example uses:","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"pgf(ctry, 1.0)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"1.0","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"pgf(ct, 1.0)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"1.0","category":"page"},{"location":"#pmf","page":"ClinicalTrialOptm.jl","title":"pmf","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"pmf","category":"page"},{"location":"#ClinicalTrialOptm.pmf","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.pmf","text":"pmf(ctry :: Country)\n\nCalculate the probability mass function of the number of patients enrolled in a center in country ctry, by Fast Fourier transform of its pgf.\n\nSee also pgf\n\n\n\n\n\npmf(ct :: ClinicalTrial)\n\nCompute the probability mass function of the number of patients enrolled in a clinical trial ct, by a Fast Fourier Transform of its pgf.\n\nSee also pgf\n\n\n\n\n\n","category":"function"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Example uses:","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"ctry_pmf = pmf(ctry)\npltrange = 0:200\n\nx = pltrange\ny = ctry_pmf[pltrange .+ 1]\n\nplot(x, y, xlabel = \"Number of patients\", ylabel = \"Probability\", legend = false)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"(Image: svg)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"ct_pmf = pmf(ct)\npltrange = 0:1000\n\nx = pltrange\ny = ct_pmf[pltrange .+ 1]\n\nplot(x, y, xlabel = \"Number of patients\", ylabel = \"Probability\", legend = false)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"(Image: svg)","category":"page"},{"location":"#cdf","page":"ClinicalTrialOptm.jl","title":"cdf","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"cdf","category":"page"},{"location":"#Distributions.cdf","page":"ClinicalTrialOptm.jl","title":"Distributions.cdf","text":"cdf(d::UnivariateDistribution, x::Real)\n\nEvaluate the cumulative probability at x.\n\nSee also ccdf, logcdf, and logccdf.\n\n\n\n\n\ncdf(d::Skellam, t::Real)\n\nImplementation based on SciPy: https://github.com/scipy/scipy/blob/v0.15.1/scipy/stats/discretedistns.py\n\nRefer to Eqn (5) in On an Extension of the Connexion Between Poisson and χ2 Distributions, N.L Johnson(1959) Vol 46, No 3/4, doi:10.2307/2333532 It relates the Skellam and Non-central chisquare PDFs, which is very similar to their CDFs computation as well.\n\nComputing cdf of the Skellam distribution.\n\n\n\n\n\ncdf(ct :: ClinicalTrial, x :: Real)\n\nCompute the cumulative distribution function of the number of patients enrolled in a clinical trial ct evaluated at x, P(X ≤ x), by applying the Gil-Pelaez inversion formula  (https://en.wikipedia.org/wiki/Characteristic_function_(probability_theory)#Inversion_formula) to the characteristic function. \n\nSee also ccdf\n\n\n\n\n\n","category":"function"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Example use:","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"x = [0:1:1000;]\ny = Vector{Float64}(undef, 1001)\n\nfor n in 1:1001\n    y[n] = cdf(ct, n - 1)\nend\n\nplot(x, y, xlabel = \"Number of Patients\", ylabel = \"Cumulative Probability\", legend = false)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"(Image: svg)","category":"page"},{"location":"#ccdf","page":"ClinicalTrialOptm.jl","title":"ccdf","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"ccdf","category":"page"},{"location":"#Distributions.ccdf","page":"ClinicalTrialOptm.jl","title":"Distributions.ccdf","text":"ccdf(d::UnivariateDistribution, x::Real)\n\nThe complementary cumulative function evaluated at x, i.e. 1 - cdf(d, x).\n\n\n\n\n\nccdf(ct :: ClinicalTrial, x :: Real)\n\nCompute the complementary distribution function of the number of patients enrolled in a  clinical trial ct evaluated by 1 - cdf(ct, x).\n\nSee also cdf\n\n\n\n\n\n","category":"function"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Example use:","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"x = [0:1:1000;]\ny = Vector{Float64}(undef, 1001)\n\nfor n in 1:1001\n    y[n] = ccdf(ct, n - 1)\nend\n\nplot(x, y, xlabel = \"Number of Patients\", ylabel = \"Complementary Cumulative Probability\", \n    legend = false)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"(Image: svg)","category":"page"},{"location":"#Example-Case","page":"ClinicalTrialOptm.jl","title":"Example Case","text":"","category":"section"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"With individual examples shown for all functions in the package, here is an example workflow of ClinicalTrialOptm.jl using estimated parameters.","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"m = [1, 1.2, 1.4, 1.9]\ns² = [1.5, 1.7, 1.3, 1.1]\nl = [0, 4, 2, 1]\nu = [10, 24, 20, 15]\nc₀ = [19000.0, 15000.0, 14000.0, 16000.0]\nc = [7000.0, 5000.0, 5000.0, 6000.0]\nq = [1000.0, 2000.0, 1500.0, 1600.0]\nd = [0.10, 0.09, 0.04, 0.07]\nT₀ = fill(Uniform(0.0, 6.0), 4)\nTd = 12.0\nct = ClinicalTrial(m, s², l, u, c₀, c, q, d, T₀, Td)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Global Clinical Trial:\n\nOptimal center assignment not computed.\n┌─────────┬─────────┬────────┬────────────┬────────────────┬──────────────┬───────────────┬──────────────┬──────────────┬─────────┐\n│\u001b[1m Country \u001b[0m│\u001b[1m mean(λ) \u001b[0m│\u001b[1m var(λ) \u001b[0m│\u001b[1m init. cost \u001b[0m│\u001b[1m    maint. cost \u001b[0m│\u001b[1m enroll. cost \u001b[0m│\u001b[1m drop out rate \u001b[0m│\u001b[1m min. centers \u001b[0m│\u001b[1m max. centers \u001b[0m│\u001b[1m centers \u001b[0m│\n│\u001b[90m         \u001b[0m│\u001b[90m         \u001b[0m│\u001b[90m        \u001b[0m│\u001b[90m   $/center \u001b[0m│\u001b[90m $/center/month \u001b[0m│\u001b[90m    $/patient \u001b[0m│\u001b[90m               \u001b[0m│\u001b[90m              \u001b[0m│\u001b[90m              \u001b[0m│\u001b[90m         \u001b[0m│\n├─────────┼─────────┼────────┼────────────┼────────────────┼──────────────┼───────────────┼──────────────┼──────────────┼─────────┤\n│       1 │    1.00 │   1.50 │      19000 │           7000 │         1000 │          0.10 │            0 │           10 │      NA │\n│       2 │    1.20 │   1.70 │      15000 │           5000 │         2000 │          0.09 │            4 │           24 │      NA │\n│       3 │    1.40 │   1.30 │      14000 │           5000 │         1500 │          0.04 │            2 │           20 │      NA │\n│       4 │    1.90 │   1.10 │      16000 │           6000 │         1600 │          0.07 │            1 │           15 │      NA │\n└─────────┴─────────┴────────┴────────────┴────────────────┴──────────────┴───────────────┴──────────────┴──────────────┴─────────┘","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Solving for the optimal trial design:","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"optdes!(ct, 400, ps = 0.95)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"Lower bound probability of success: 7.112034326831168e-28\nThe optimal solution is not the lower bound of the centers.\nUpper bound probability of success: 0.9999940867534107\nThe optimal solution is feasible.\n  0.026162 seconds (2.16 k allocations: 57.289 KiB)\nsolution_summary(model) = * Solver : SCIP\n\n* Status\n  Result count       : 7\n  Termination status : OPTIMAL\n  Message from the solver:\n  \"SCIP_STATUS_OPTIMAL\"\n\n* Candidate solution (result #1)\n  Primal status      : FEASIBLE_POINT\n  Dual status        : NO_SOLUTION\n  Objective value    : 3.27739e+06\n  Objective bound    : 3.27739e+06\n  Relative gap       : 0.00000e+00\n\n* Work counters\n  Solve time (sec)   : 2.14900e-02\n  Simplex iterations : 89\n  Node count         : 28\n\ntermination_status(model) = MathOptInterface.OPTIMAL\nprimal_status(model) = MathOptInterface.FEASIBLE_POINT\nobjective_value(model) = 3.277387200000001e6\n\n\n\n\n\n\nGlobal Clinical Trial:\n\nOptimal center assignment calculated.\nAn optimal solution has been found.\nTotal duration (months): 12.0\nTarget enrollment: 400\nProbability of success (based on normal approximation): 0.9550669693996681\nProbability of success (based on Poisson-Gamma model): 0.9636054293917478\nExpected cost ($): 3.2773872e6\n┌─────────┬─────────┬────────┬────────────┬────────────────┬──────────────┬───────────────┬──────────────┬──────────────┬─────────┐\n│\u001b[1m Country \u001b[0m│\u001b[1m mean(λ) \u001b[0m│\u001b[1m var(λ) \u001b[0m│\u001b[1m init. cost \u001b[0m│\u001b[1m    maint. cost \u001b[0m│\u001b[1m enroll. cost \u001b[0m│\u001b[1m drop out rate \u001b[0m│\u001b[1m min. centers \u001b[0m│\u001b[1m max. centers \u001b[0m│\u001b[1m centers \u001b[0m│\n│\u001b[90m         \u001b[0m│\u001b[90m         \u001b[0m│\u001b[90m        \u001b[0m│\u001b[90m   $/center \u001b[0m│\u001b[90m $/center/month \u001b[0m│\u001b[90m    $/patient \u001b[0m│\u001b[90m               \u001b[0m│\u001b[90m              \u001b[0m│\u001b[90m              \u001b[0m│\u001b[90m         \u001b[0m│\n├─────────┼─────────┼────────┼────────────┼────────────────┼──────────────┼───────────────┼──────────────┼──────────────┼─────────┤\n│       1 │    1.00 │   1.50 │      19000 │           7000 │         1000 │          0.10 │            0 │           10 │       0 │\n│       2 │    1.20 │   1.70 │      15000 │           5000 │         2000 │          0.09 │            4 │           24 │       5 │\n│       3 │    1.40 │   1.30 │      14000 │           5000 │         1500 │          0.04 │            2 │           20 │      20 │\n│       4 │    1.90 │   1.10 │      16000 │           6000 │         1600 │          0.07 │            1 │           15 │      14 │\n└─────────┴─────────┴────────┴────────────┴────────────────┴──────────────┴───────────────┴──────────────┴──────────────┴─────────┘","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"From optdes!, we determine that 0 centers should be placed in country 1, 5 centers should be placed in country 2, 20 in country 3, and 14 in country 4 to minimize the trial cost and ensure a high success rate with recruitment.","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"We can also visualize how the normal approximation of total patient recruitment compares to the more accurate Poisson-Gamma model by plotting the normal distribution, with mean(ct) and var(ct), against pmf(ct). ","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"# normal approximation\nμ = mean(ct)\nσ = sqrt(var(ct))\nplot(Normal(μ, σ), label = \"Normal Approximation\")\n\n# adding the poisson-gamma model\nct_pmf = pmf(ct)\npltrange = 200:800\n\nx = pltrange\ny = ct_pmf[pltrange .+ 1]\n\nplot!(x, y, xlabel = \"Number of patients\", ylabel = \"Probability\", \n    label = \"Possion-Gamma Model\")","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"(Image: svg)","category":"page"},{"location":"","page":"ClinicalTrialOptm.jl","title":"ClinicalTrialOptm.jl","text":"From the graph, the functions are very similar, demonstrating that the normal approximation works well in calculating the patient enrollment.","category":"page"}]
}
