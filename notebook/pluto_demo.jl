### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ c973950e-6c58-449f-a72c-c426b432f61a
begin
	import Pkg
	Pkg.activate(pwd())
	using ClinicalTrialOptm, Distributions, PlutoUI
	import PlutoUI: combine
end

# ╔═╡ 850a6901-e05c-4e90-b8b0-e2f8d4786312
md"""# Optimal Design for Global Clinical Trial
"""

# ╔═╡ 3a68db76-5e79-4991-80c9-371a29e51d8b
TableOfContents(depth = 4)

# ╔═╡ 017c0976-b206-431e-8119-eeaced720a75
md"""## Input
"""

# ╔═╡ 15b1851c-5e1e-4f1d-bfbf-b0b7827ba3a3
md"""### How many countries?
"""

# ╔═╡ e75d917f-2125-4c92-b5c1-9ebc3cd1b4f7
@bind J NumberField(1:20, default=5)

# ╔═╡ 262d09f5-0793-416a-aafd-f2e144104a0f
md"""
### Clinical trial duration (months)
"""

# ╔═╡ 443420d3-b304-468a-b8bc-824c59665c2f
@bind Td NumberField(6.0:60.0, default=15)

# ╔═╡ 3bfc13e5-0b22-409d-9001-db0b66feea33
md"""
### Target enrollment
"""

# ╔═╡ cbd8e610-3142-4f8c-a34a-046ac204a25b
@bind ntarget NumberField(10:10:1_000, default=500)

# ╔═╡ b18a9188-89f2-45ab-823f-d557b544b3da
md"""
### Desired probability of success (actual enrollment ≥ $(ntarget))
"""

# ╔═╡ 1791d352-abee-4e45-b8ac-045cd526f995
@bind ps Slider(0.01:0.01:0.99, default = 0.95, show_value = true)

# ╔═╡ d35d8f48-0140-4c7b-aaec-dc6aa0fcb9e9
md"""
### Country-wise parameters
"""

# ╔═╡ a417af25-d3af-4652-9be8-b9bf8cbc540f
md"Randomize $(@bind randomize_m CheckBox(default = true))"

# ╔═╡ 3345b353-1692-40c4-95ba-0b3d3cf038f0
function m_input(J::Integer, randomize::Bool)
	
	return combine() do Child
		rvar = "randomize_m"
		mvar = "m" .* string.(1:J)
		inputs = [
			md""" Country $(j): $(
				Child(mvar[j], Slider(1.0:0.01:2.0, default = randomize_m ? rand(Uniform(1.0, 2.0)) : 1.5, show_value = true))
			)"""
			
			for j in 1:J
		]
		
		md"""
		#### Mean enrollment rate in each country
		$(inputs)
		"""
	end
end

# ╔═╡ 1656d80f-06d2-437c-b047-9b94401944a0
@bind mtuple m_input(J, randomize_m)

# ╔═╡ 5b226cd4-ce9f-4d84-a9fa-ad826802bb78
m = [mtuple[j] for j in 1:J]

# ╔═╡ 0286e0cc-ffda-48a0-b1a1-1cdcf5e69bec
md"Randomize $(@bind randomize_s² CheckBox(default = true))"

# ╔═╡ 77771f32-a251-43ca-aa3a-fbe4087fb7cf
function s²_input(J::Integer, random::Bool)
	
	return combine() do Child

		s²var = "m" .* string.(1:J)
		inputs = [
			md""" Country $(j): $(
				Child(s²var[j], Slider(1.0:0.01:2.0, default = randomize_s² ? rand(Uniform(1.0, 2.0)) : 1.5, show_value = true))
			)"""
			
			for j in 1:J
		]
		
		md"""
		#### Variance of enrollment rate in each country
		$(inputs)
		"""
	end
end

# ╔═╡ 4d19e48d-3b59-4fa0-b681-de403c8282db
@bind s²tuple s²_input(J, randomize_s²)

# ╔═╡ c8199643-51d4-4c1e-8df8-ac60178c094b
s² = [s²tuple[j] for j in 1:J]

# ╔═╡ 60d0afbb-90b9-4bde-b574-71aa851f97a8
md"Randomize $(@bind randomize_c₀ CheckBox(default = true))"

# ╔═╡ 644b5441-8ca2-4ca2-8a88-8cc18f184a72
function c₀_input(J::Integer, random::Bool)
	
	return combine() do Child

		c₀var = "m" .* string.(1:J)
		inputs = [
			md""" Country $(j): $(
				Child(c₀var[j], Slider(15_000.0:1_000.0:25_000.0, default = randomize_c₀ ? rand(15_000.0:1_000.0:25_000.0) : 20_000.0, show_value = true))
			)"""
			
			for j in 1:J
		]
		
		md"""
		#### Center initialization cost (\$) in each country
		$(inputs)
		"""
	end
end

# ╔═╡ 6b748bf1-5636-44af-9b5a-c60a6a647fe1
@bind c₀tuple c₀_input(J, randomize_c₀)

# ╔═╡ 2b261b6b-98f8-440b-836e-f0c3ee5d63cd
c₀ = [c₀tuple[j] for j in 1:J]

# ╔═╡ bc54b570-d7e1-4838-a09d-527c527ff628
md"Randomize $(@bind randomize_c CheckBox(default = true))"

# ╔═╡ 44c187ec-97f6-48a1-abf6-940e4e03d591
function c_input(J::Integer, random::Bool)
	
	return combine() do Child

		cvar = "m" .* string.(1:J)
		inputs = [
			md""" Country $(j): $(
				Child(cvar[j], Slider(4_000.0:100.0:5_000.0, default = random ? rand(4_000.0:100.0:5_000.0) : 4500, show_value = true))
			)"""
			
			for j in 1:J
		]
		
		md"""
		#### Center maintainance cost (\$/month) in each country
		$(inputs)
		"""
	end
end

# ╔═╡ dbd0274a-ce85-4454-9806-e257bcf4625d
@bind ctuple c_input(J, randomize_c)

# ╔═╡ e5132f3e-d7bd-4733-b8a4-8004bd43c254
c = [ctuple[j] for j in 1:J]

# ╔═╡ fb6092cd-14d4-40b5-bea9-d5279a903e1d
md"""
#### Center intialization time

Assumed to be uniform in the first 6 months, but can be adjusted.
"""

# ╔═╡ 496fba51-bd42-465c-8f16-2e5f2cb6170f
T₀ = [Uniform(0, 6.0) for j in 1:J]

# ╔═╡ 20ca537a-3cb0-4a7a-8613-6ab58f58ddc8
md"Randomize $(@bind randomize_l CheckBox(default = true))"

# ╔═╡ b01eca63-1908-43bd-b371-24d41fc9f6dd
function l_input(J::Integer, random::Bool)
	
	return combine() do Child

		lvar = "m" .* string.(1:J)
		inputs = [
			md""" Country $(j): $(
				Child(lvar[j], Slider(1:5, default = random ? rand(1:3) : 1, show_value = true))
			)"""
			
			for j in 1:J
		]
		
		md"""
		#### Minimum number of centers in each country
		$(inputs)
		"""
	end
end

# ╔═╡ d0b7dd28-7e5f-4a80-ba5d-89b6e0cfeb72
@bind ltuple l_input(J, randomize_l)

# ╔═╡ 02c89e40-d8e6-40a9-b696-077894d8d6ce
l = [ltuple[j] for j in 1:J]

# ╔═╡ 8a03f26d-fcfa-469f-89a6-5158a39a1895
md"Randomize $(@bind randomize_u CheckBox(default = false))"

# ╔═╡ 6fede67a-a7d6-4d39-a563-c2d97eb5a1fc
function u_input(J::Integer, l::Vector, random::Bool)
	
	return combine() do Child

		uvar = "u" .* string.(1:J)
		inputs = [
			md""" Country $(j): $(
				Child(uvar[j], Slider(l[j]:l[j]+99, default = random ? l[j] + rand(50:100) : l[j]+100, show_value = true))
			)"""
			
			for j in 1:J
		]
		
		md"""
		#### Maximum number of centers in each country
		$(inputs)
		"""
	end
end

# ╔═╡ 4d7161a8-85c1-4ab8-8e54-988a94522e3e
@bind utuple u_input(J, l, randomize_u)

# ╔═╡ 9448d5f7-861e-490b-9758-2e1438455892
u = [utuple[j] for j in 1:J]

# ╔═╡ 25d85555-488a-40df-94fb-9633946dabf3
md"Randomize $(@bind randomize_q CheckBox(default = true))"

# ╔═╡ aa2c1140-86d7-4c68-bd21-def6324768a5
function q_input(J::Integer, random::Bool)
	
	return combine() do Child

		qvar = "m" .* string.(1:J)
		inputs = [
			md""" Country $(j): $(
				Child(qvar[j], Slider(1_000:100:2_000, default = random ? rand(1_000:100:2_000) : 1500, show_value = true))
			)"""
			
			for j in 1:J
		]
		
		md"""
		#### Cost (\$) per enrolled patient
		$(inputs)
		"""
	end
end

# ╔═╡ 6e847b99-c823-4839-8595-7c3d0c059816
@bind qtuple q_input(J, randomize_q)

# ╔═╡ 8cc28e2d-913b-44bf-b4e7-7c8cd55f372c
q = [Float64(qtuple[j]) for j in 1:J]

# ╔═╡ 7f175f8c-e0e2-422e-84f2-72c2c1f9cb17
md"Randomize $(@bind randomize_d CheckBox(default = false))"

# ╔═╡ e0ab2d39-bc95-4f31-8cc0-592d65a794d6
function d_input(J::Integer, random::Bool)
	
	return combine() do Child

		dvar = "m" .* string.(1:J)
		inputs = [
			md""" Country $(j): $(
				Child(dvar[j], Slider(0.00:0.01:0.99, default = random ? rand(0.1:0.02:0.2) : 0.0, show_value = true))
			)"""
			
			for j in 1:J
		]
		
		md"""
		#### Drop out rate in each country
		$(inputs)
		"""
	end
end

# ╔═╡ f1d7f2df-abe3-4461-bbb4-39db1794e2d2
@bind dtuple d_input(J, randomize_d)

# ╔═╡ 0a5f89f8-a6d8-4fd9-b643-a9e54b033a04
d = [dtuple[j] for j in 1:J]

# ╔═╡ 36e2afa1-46a3-4843-bccb-885b4938b746
md"""
## Find optimal design by MIP
"""

# ╔═╡ ae15fbd3-8d3d-4095-860a-4c973a18033a
ct = ClinicalTrial(m, s², l, u, c₀, c, q, d, T₀, Td, l)

# ╔═╡ 548f5cb3-f443-4f09-8f29-0f6bdc2ae23f
optdes!(ct, ntarget, ps = ps)

# ╔═╡ Cell order:
# ╟─850a6901-e05c-4e90-b8b0-e2f8d4786312
# ╠═3a68db76-5e79-4991-80c9-371a29e51d8b
# ╠═c973950e-6c58-449f-a72c-c426b432f61a
# ╟─017c0976-b206-431e-8119-eeaced720a75
# ╟─15b1851c-5e1e-4f1d-bfbf-b0b7827ba3a3
# ╟─e75d917f-2125-4c92-b5c1-9ebc3cd1b4f7
# ╟─262d09f5-0793-416a-aafd-f2e144104a0f
# ╟─443420d3-b304-468a-b8bc-824c59665c2f
# ╟─3bfc13e5-0b22-409d-9001-db0b66feea33
# ╟─cbd8e610-3142-4f8c-a34a-046ac204a25b
# ╟─b18a9188-89f2-45ab-823f-d557b544b3da
# ╟─1791d352-abee-4e45-b8ac-045cd526f995
# ╟─d35d8f48-0140-4c7b-aaec-dc6aa0fcb9e9
# ╟─1656d80f-06d2-437c-b047-9b94401944a0
# ╟─a417af25-d3af-4652-9be8-b9bf8cbc540f
# ╟─3345b353-1692-40c4-95ba-0b3d3cf038f0
# ╟─5b226cd4-ce9f-4d84-a9fa-ad826802bb78
# ╟─4d19e48d-3b59-4fa0-b681-de403c8282db
# ╟─0286e0cc-ffda-48a0-b1a1-1cdcf5e69bec
# ╟─77771f32-a251-43ca-aa3a-fbe4087fb7cf
# ╟─c8199643-51d4-4c1e-8df8-ac60178c094b
# ╟─6b748bf1-5636-44af-9b5a-c60a6a647fe1
# ╟─60d0afbb-90b9-4bde-b574-71aa851f97a8
# ╟─644b5441-8ca2-4ca2-8a88-8cc18f184a72
# ╟─2b261b6b-98f8-440b-836e-f0c3ee5d63cd
# ╟─dbd0274a-ce85-4454-9806-e257bcf4625d
# ╟─bc54b570-d7e1-4838-a09d-527c527ff628
# ╟─44c187ec-97f6-48a1-abf6-940e4e03d591
# ╟─e5132f3e-d7bd-4733-b8a4-8004bd43c254
# ╟─fb6092cd-14d4-40b5-bea9-d5279a903e1d
# ╟─496fba51-bd42-465c-8f16-2e5f2cb6170f
# ╟─d0b7dd28-7e5f-4a80-ba5d-89b6e0cfeb72
# ╟─20ca537a-3cb0-4a7a-8613-6ab58f58ddc8
# ╟─b01eca63-1908-43bd-b371-24d41fc9f6dd
# ╟─02c89e40-d8e6-40a9-b696-077894d8d6ce
# ╟─4d7161a8-85c1-4ab8-8e54-988a94522e3e
# ╠═8a03f26d-fcfa-469f-89a6-5158a39a1895
# ╟─6fede67a-a7d6-4d39-a563-c2d97eb5a1fc
# ╟─9448d5f7-861e-490b-9758-2e1438455892
# ╟─6e847b99-c823-4839-8595-7c3d0c059816
# ╟─25d85555-488a-40df-94fb-9633946dabf3
# ╟─aa2c1140-86d7-4c68-bd21-def6324768a5
# ╟─8cc28e2d-913b-44bf-b4e7-7c8cd55f372c
# ╟─f1d7f2df-abe3-4461-bbb4-39db1794e2d2
# ╠═7f175f8c-e0e2-422e-84f2-72c2c1f9cb17
# ╟─e0ab2d39-bc95-4f31-8cc0-592d65a794d6
# ╟─0a5f89f8-a6d8-4fd9-b643-a9e54b033a04
# ╟─36e2afa1-46a3-4843-bccb-885b4938b746
# ╠═ae15fbd3-8d3d-4095-860a-4c973a18033a
# ╠═548f5cb3-f443-4f09-8f29-0f6bdc2ae23f
