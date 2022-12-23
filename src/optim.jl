"""
    optdes!(ct, ntarget, ps, solver = SCIP)

Find the optimal clinical trial design, with `ntarget` enrollment and guaranteed 
probability of success `ps`. Default `solver` is SCIP, can be changed to KNITRO by using 
"KNITRO". Overwrite `ct.centers` by the optimal solution.
"""
function optdes!(
    ct :: ClinicalTrial,
    ntarget :: Integer;
    ps :: Real = 0.5,
    solver = optimizer_with_attributes(
        SCIP.Optimizer, "display/verblevel"=>0, "limits/gap"=>0.01 
    )
    )
    # test whether the lb is optimal
    if (lbtest(ct, ntarget, ps = ps) == :lb_optimal)
        return
    else
        # test whether the problem is infeasible by checking upper bound
        if (ubtest(ct, ntarget, ps = ps) == :ub_infeasible)
            return
        else
            # MIP
            ct.ntarget[1] = ntarget
            # number of countries
            J = length(ct.countries)
            # coefs for linear objective function
            c_f = [mean_cost(ct.countries[j]) for j in 1:J]
            # coefs for mean of total enrollment
            c_μ = [mean(ct.countries[j]) for j in 1:J]
            # coefs for variance of total enrollment
            c_σ² = [var(ct.countries[j]) for j in 1:J]
            # Φ⁻¹(ϵ)
            c_ϵ = quantile(Normal(), 1 - ps)
            if c_ϵ ≥ 0
                # # set up Pajarito solvers    
                # oa_solver = optimizer_with_attributes(HiGHS.Optimizer,
                #     MOI.Silent() => true,
                #     "mip_feasibility_tolerance" => 1e-8,
                #     "mip_rel_gap" => 1e-6,
                # )
                # # conic_solver = optimizer_with_attributes(Hypatia.Optimizer, 
                # #     MOI.Silent() => true,
                # # )
                # conic_solver = optimizer_with_attributes(SCS.Optimizer, 
                #     MOI.Silent() => true,
                # )
                # opt = optimizer_with_attributes(Pajarito.Optimizer,
                #     "time_limit" => 600, 
                #     "oa_solver" => oa_solver, 
                #     "conic_solver" => conic_solver,
                # )
                # set up KNITRO solver
                opt = optimizer_with_attributes(
                    KNITRO.Optimizer
                )
                # set up model
                model = Model(opt)
                @variable(
                    model, 
                    ct.countries[j].l ≤ x[j = 1:J] ≤ ct.countries[j].u,
                    integer = true
                )
                @variable(model, μ)
                @variable(model, σ²)
                @objective(model, Min, sum(c_f[j] * x[j] for j in 1:J))
                @constraint(model, μ == sum(c_μ[j] * x[j] for j in 1:J))
                @constraint(model, σ² == sum(c_σ²[j] * x[j] for j in 1:J))
                @constraint(model, [σ²; 0.5abs2(c_ϵ); ntarget - μ] in RotatedSecondOrderCone())
            else # c_ϵ < 0
                if solver == "KNITRO"
                    # set up KNITRO solver
                    minlp = optimizer_with_attributes(
                        KNITRO.Optimizer
                    )
                    # set up model
                    model = Model(minlp)
                    @variable(
                        model, 
                        ct.countries[j].l ≤ x[j = 1:J] ≤ ct.countries[j].u,
                        integer = true
                    )
                    @variable(model, μ ≥ ntarget)
                    @variable(model, σ²)
                    @objective(model, Min, sum(c_f[j] * x[j] for j in 1:J))
                    @constraint(model, μ == sum(c_μ[j] * x[j] for j in 1:J))
                    @constraint(model, σ² == sum(c_σ²[j] * x[j] for j in 1:J))
                    @NLconstraint(model, abs2(c_ϵ) * σ² ≤ abs2(ntarget - μ))
                    # @NLconstraint(model, μ + c_ϵ * sqrt(σ²) ≥ ntarget)
                else # solver defaults to SCIP
                    minlp = solver
                    model = Model(minlp)
                    @variable(
                        model,
                        ct.countries[j].l ≤ x[j = 1:J] ≤ ct.countries[j].u,
                        integer = true
                    )
                    @variable(model, μ ≥ ntarget)
                    @variable(model, σ²)
                    @objective(model, Min, sum(c_f[j] * x[j] for j in 1:J))
                    @constraint(model, μ == sum(c_μ[j] * x[j] for j in 1:J))
                    @constraint(model, σ² == sum(c_σ²[j] * x[j] for j in 1:J))
                    @NLconstraint(model, (c_ϵ)^2 * σ² ≤ (ntarget - μ)^2)
                     # @NLconstraint(model, μ + c_ϵ * sqrt(σ²) ≥ ntarget)
                end
            end
            # solvers
            @time optimize!(model)
            @show solution_summary(model)
            ct.centers .= round.(Int, value.(x))
            @show termination_status(model)
            @show primal_status(model)
            @show objective_value(model)
            ct.isoptm[1] = true
            if (termination_status(model) == OPTIMAL || termination_status(model) == LOCALLY_SOLVED)
                ct.solution_status[1] = "An optimal solution has been found."
            else
                ct.solution_status[1] = "The solution is infeasible with this solver."
            end
            ct
        end
    end
end

"""
    lbtest(ct :: ClinicalTrial, ntarget :: Integer; ps :: Real = 0.5)

Tests the lower bound of the center restraints for a clinical trial `ct`.
`ntarget` is the target patient enrollment and `ps` is the desired probability
of success, with a default value of 0.5. The test is calculated using the normal 
approximation of the PoS, and it outputs whether or not the optimal solution is 
the lower bound of the centers.
"""
function lbtest(
    ct :: ClinicalTrial,
    ntarget :: Integer;
    ps :: Real = 0.5
    )
    l = [ct.countries[j].l for j in 1:length(ct.countries)]
    copyto!(ct.centers, l)
    μ, σ² = mean(ct), var(ct)
    PoS = (ccdf(Normal(μ, sqrt(σ²)), ntarget))
    println("Lower bound probability of success: $PoS")
    if PoS > ps
        println("The optimal solution is the lower bound of the centers.")
        return :lb_optimal
    else 
        println("The optimal solution is not the lower bound of the centers.")
        return :lb_not_optimal
    end
end

"""
ubtest(ct :: ClinicalTrial, ntarget :: Integer; ps :: Real = 0.5)

Tests the upper bound of the center restraints for a clinical trial `ct`.
`ntarget` is the target patient enrollment and `ps` is the desired probability
of success, with a default value of 0.5. The test is calculated using the normal 
approximation of the PoS, and it outputs whether or not a solution is feasible.
"""
function ubtest(
    ct :: ClinicalTrial,
    ntarget :: Integer;
    ps :: Real = 0.5
    )
    u = [ct.countries[j].u for j in 1:length(ct.countries)]
    copyto!(ct.centers, u)
    μ, σ² = mean(ct), var(ct)
    PoS = (ccdf(Normal(μ, sqrt(σ²)), ntarget))
    println("Upper bound probability of success: $PoS")
    if PoS < ps
        println("The optimal solution is infeasible.")
        return :ub_infeasible
    else 
        println("The optimal solution is feasible.")
        return :ub_feasible
    end
end
