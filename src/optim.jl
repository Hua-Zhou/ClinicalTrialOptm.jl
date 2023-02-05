oa_solver = optimizer_with_attributes(HiGHS.Optimizer,
MOI.Silent() => true,
"mip_feasibility_tolerance" => 1e-8,
"mip_rel_gap" => 1e-6,
)
conic_solver = optimizer_with_attributes(SCS.Optimizer, 
MOI.Silent() => true,
)
convex_default_solver = optimizer_with_attributes(Pajarito.Optimizer,
    "time_limit" => 600, 
    "oa_solver" => oa_solver, 
    "conic_solver" => conic_solver,
    "tol_rel_gap" => 0.001
)

nonconvex_default_solver =  optimizer_with_attributes(
    SCIP.Optimizer, "display/verblevel"=>0, "limits/gap"=>0.001 
)

"""
    optdes!(ct :: ClinicalTrial, ntarget :: Integer, ps :: Real, solver)

Find the optimal clinical trial design, with `ntarget` enrollment and guaranteed 
probability of success `ps`.

The solver can be changed to the user's preferance by setting `solver` equal to your solver choice (see documentation page for examples). 
The function overwrites `ct.centers` by the optimal solution.

See also [`lbtest`](@ref), [`ubtest`](@ref), [`ClinicalTrial`](@ref)
"""
function optdes!(
    ct :: ClinicalTrial,
    ntarget :: Integer;
    ps :: Real = 0.5,
    solver = ifelse(ps <= 0.5, convex_default_solver, nonconvex_default_solver)
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
            # convex case
            if c_ϵ ≥ 0
                model = Model(solver)
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
            # non-convex case
            else # c_ϵ < 0
                model = Model(solver)
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

Test the lower bound of the center restraints for a clinical trial `ct` to see if it is already the optimal solution.

The test is calculated using the normal approximation of the PoS, and `ntarget` is the target patient enrollment and `ps` is the desired probability
of success.

See also [`ubtest`](@ref), [`optdes!`](@ref)
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

Test the upper bound of the center restraints for a clinical trial `ct` to see if a solution is feasible.

The test is calculated using the normal approximation of the PoS, and `ntarget` is the target patient enrollment and `ps` is the desired probability
of success.

See also [`lbtest`](@ref), [`optdes!`](@ref)
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
