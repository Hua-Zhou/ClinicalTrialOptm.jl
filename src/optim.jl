"""
    optdes!(ct, ntarget, ps)

Find the optimal clinical trial design, with `ntarget` enrollment and guaranteed 
probability of success `ps`. Overwrite `ct.centers` by the optimal solution.
"""
function optdes!(
    ct :: ClinicalTrial,
    ntarget :: Integer;
    ps :: Real = 0.5
    )
    # test whether the lb is optimal
    (lbtest(ct, ntarget, ps = ps) == :lb_optimal) && (return :opm_successful)
    # test whether the problem is infeasible by checking upper bound
    (ubtest(ct, ntarget, ps = ps) == :ub_infeasible) && (return :opm_infeasible)
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
        # # set up Juniper solver
        # nl_solver = optimizer_with_attributes(
        #     Ipopt.Optimizer, 
        #     MOI.Silent() => true, 
        #     "sb" => "yes", 
        #     "max_iter"   => 9999
        # )
        # mip_solver = optimizer_with_attributes(
        #     HiGHS.Optimizer, 
        #     "output_flag" => false
        # )
        # # mip_solver = optimizer_with_attributes(
        # #     Cbc.Optimizer, 
        # #     "logLevel" => 0
        # # )
        # minlp = optimizer_with_attributes(
        #     Juniper.Optimizer, 
        #     "mip_solver" => mip_solver, 
        #     "nl_solver" => nl_solver,
        #     # "traverse_strategy" => :DBFS,
        #     # "processors" => 4
        # )
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
    end
    # solvers
    @time optimize!(model)
    @show solution_summary(model)
    ct.centers .= round.(Int, value.(x))
    @show termination_status(model)
    @show primal_status(model)
    @show objective_value(model)
    ct.isoptm[1] = true
    ct
end

"""
    TODO
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
    println("Probability of success: $PoS")
    if PoS > ps
        println("The optimal solution is the lower bound of the centers.")
        return :lb_optimal
    else 
        println("The optimal solution is not the lower bound of the centers.")
        return :lb_not_optimal
    end
end

"""
    TODO
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
    println("Probability of success: $PoS")
    if PoS < ps
        println("The optimal solution is infeasible.")
        return :ub_infeasible
    else 
        println("The optimal solution is feasible.")
        return :ub_feasible
    end
end
