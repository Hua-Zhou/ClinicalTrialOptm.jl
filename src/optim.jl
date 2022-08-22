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
    # set up solvers
    oa_solver = optimizer_with_attributes(HiGHS.Optimizer,
        MOI.Silent() => true,
        "mip_feasibility_tolerance" => 1e-8,
        "mip_rel_gap" => 1e-6,
    )
    conic_solver = optimizer_with_attributes(Hypatia.Optimizer, 
        MOI.Silent() => true,
    )
    opt = optimizer_with_attributes(Pajarito.Optimizer,
        "time_limit" => 600, 
        "oa_solver" => oa_solver, 
        "conic_solver" => conic_solver,
    )
    # set up model
    model = Model(opt)
    @variable(
        model, 
        ct.countries[i].l ≤ x[i = 1:J] ≤ ct.countries[i].u,
        integer = true
    )
    @variable(model, μ)
    @variable(model, σ²)
    @constraint(model, μ == sum(c_μ[j] * x[j] for j in 1:J))
    @constraint(model, σ² == sum(c_σ²[j] * x[j] for j in 1:J))
    @objective(model, Min, sum(c_f[j] * x[j] for j in 1:J))
    if c_ϵ ≥ 0
        @constraint(model, [σ²; 0.5abs2(c_ϵ); ntarget - μ] in RotatedSecondOrderCone())
    else
        @NLconstraint(model, μ + c_ϵ * sqrt(σ²) - ntarget ≥ 0)
    end
    # solvers
    @time optimize!(model)
    @show solution_summary(model)
    @show ct.centers .= value.(x)
    @show termination_status(model)
    @show primal_status(model)
    @show objective_value(model)
    ct.isoptm[1] = true
    ct
end
