### A Pluto.jl notebook ###
# v0.11.5

using Markdown
using InteractiveUtils

# ╔═╡ dfc4e82a-db99-11ea-0166-9b9f74fe0999
md"## SIR model with births and deaths"

# ╔═╡ 3758086a-db9a-11ea-1102-bbac22d924a1
md"visualize solution in phase plane"

# ╔═╡ 9fbb2a38-db9d-11ea-26df-2789729490b3
md"## stochastic simulation of the SIR model"

# ╔═╡ d30011b2-db9e-11ea-0a84-77f66c5de33b
md"## draw a network"

# ╔═╡ 51c07944-db9f-11ea-35e3-d9ab257b6bd6
begin
	# save graph	
# 	using Compose, Cairo
#	Compose.draw(PDF("super_spreader_graph.pdf", 16cm, 16cm), gp)
end

# ╔═╡ 6da1a832-db99-11ea-0098-9b5872e94312
begin
	using DifferentialEquations, PyPlot, Statistics, LightGraphs, Printf, GraphPlot, Colors
	using PyCall
	PyPlot.matplotlib.style.use("grandbudapest.mplstyle")   
	PyPlot.matplotlib.font_manager.fontManager.addfont("OpenSans-Regular.ttf")
	#@pyimport mpl_toolkits.axes_grid1 as mpl_toolkits_axes_grid1
	mpl_toolkits_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
end

# ╔═╡ e6413104-db99-11ea-3d6e-cfa983139ddc
begin
	R₀ = 2.0                # basic reproductive number
	μ_ovr_γ_plus_μ = 0.05   # μ / (μ + γ)
	
	function update_g!(f, u, p, t)
	    # http://www.maths.usyd.edu.au/u/marym/populations/hethcote.pdf
	    s = u[1]
	    i = u[2]
	    
	    # update f
	    f[1] = -R₀ * s * i + μ_ovr_γ_plus_μ - μ_ovr_γ_plus_μ * s
	    f[2] = i * (R₀ * s - 1.0)
	end
	
	# initial condition
	ϵ = 10^-5
	u0 = [1-ϵ; ϵ]
	
	# define the ODE problem
	time_span = (0.0, 250.0)
	prob = ODEProblem(update_g!, u0, time_span)
	
	# numerically solve ODE
	sol = solve(prob)
end

# ╔═╡ 10f1c0da-db9a-11ea-1f3b-41cfa50c451d
begin
	t = range(0.0, time_span[2], length=1500)
	
	s = [sol.(t_j)[1] for t_j in t]
	i = [sol.(t_j)[2] for t_j in t]
	
	figure()
	scatter([1 / R₀], [μ_ovr_γ_plus_μ * (R₀ - 1) / R₀], marker="x", 
		zorder=300, color="r", label="endemic\nequilibrium", s=20, lw=1.5)
	xlabel(L"$[$S$](t)$")
	ylabel(L"$[$I$](t)$")
	
	legend(loc="best")
	cmap = PyPlot.matplotlib.cm.get_cmap("viridis")
	# map a time to a color
	function time_to_color(t::Float64, t_max::Float64)
		if t > t_max
			return cmap(1.0)
		else
			return cmap(t / t_max)
		end
	end
	
	t_max = 100.0
	for k = 1:length(s)-1
	    plot(s[k:k+1], i[k:k+1],
	        color=time_to_color((t[k] + t[k+1]) / 2, t_max)
	        )
	end
	
	plt.gca().set_aspect("equal")	
	ylim([-1e-2, 0.4+1e-2])
	xlim(xmax=1.0+1e-3)
	
	# info
    bbox_props = Dict(:boxstyle=>"round", :fc=>(1.0, 0.98, 0.98), :ec=>"0.5")
        #:ec=>"0.5")#, :alpha=0.9)
    sim_settings = L"$\beta (\mu+\gamma)^{-1}=2$" * "\n" * L"$\mu (\mu+\gamma)^{-1}=0.05$" * "\n" * L"$[$I$]_0=10^{-5}$"
    @assert u0[2] ≈ 10^(-5)
    @assert R₀ ≈ 2.0
	@assert μ_ovr_γ_plus_μ ≈ 0.05
	
    text(0.3, 0.3,
        sim_settings, bbox=bbox_props,
        va="center")

	ax = plt.gca()
	divider = mpl_toolkits_axes_grid1.make_axes_locatable(ax)
	cax = divider.append_axes("right", size="2%", pad=0.05)

	sm = PyPlot.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0.0, vmax=t_max))
	colorbar(sm, cax=cax, label="non-dimensional time" * "\n" * L"$(\gamma+\mu)t$", extend="max")
	tight_layout()
	savefig("sir_with_demographics.pdf", format="pdf", bbox_inches="tight")
	gcf()
end

# ╔═╡ 684f7776-db9e-11ea-2dbc-85773151febb
begin
	γ = 1
	β = R₀ * γ
	
	Δt = 0.0001
	Nt = Int(10.0 / Δt / γ)
	
	t_s = [Δt * (i - 1) for i = 1:Nt]
	
	S = zeros(Nt)
	I = zeros(Nt)
	R = zeros(Nt)
	
	N = 50
	S[1] = N - 1
	I[1] = 1
	R[1] = 0
	
	for j = 2:Nt
	    S[j] = S[j - 1]
	    I[j] = I[j - 1]
	    R[j] = R[j - 1]
	
	    if rand() < β * S[j - 1] * I[j - 1] * Δt / N
	        S[j] -= 1
	        I[j] += 1
	    end
	    if rand() < γ * I[j-1] * Δt
	        I[j] -= 1
	        R[j] += 1
	    end
	end
end

# ╔═╡ 8cfa38e8-db9e-11ea-09a2-a366b1fcc94d
begin
	sir_to_color = Dict("S" => "C3", "I" => "C5", "R" => "C0")
	
	figure()
	plot(t_s, S, label=L"S$(t)$", color=sir_to_color["S"], clip_on=false)
	plot(t_s, I, label=L"I$(t)$", color=sir_to_color["I"], clip_on=false)
	plot(t_s, R, label=L"R$(t)$", color=sir_to_color["R"], clip_on=false)
	
	xlabel(L"non-dimensional time, $\gamma t$")
	ylabel("# individuals")
	legend()
	
	text(6.0, 25.0,
        L"\mathcal{R}_0=2", bbox=bbox_props,
        va="center")
	# dy = 0.04
	# x_pos = 21.0
	# text(x_pos, s[end] + dy, L"$[$S$](t)$", color="C0")
	# text(x_pos, i[end] + dy, L"$[$I$](t)$", color="C2")
	# text(x_pos, r[end] + dy, L"$[$R$](t)$", color="C1")
	# # plot(t, -s .+ log.(s) / R0 .+ 1) # check phase plane formula
	# f(x) = log(x) - R0 * (x - 1)
	# s∞ = fzero(f, 0.1)
	# axhline(y=s∞) # check peak prevalence
	# check final size formula
	# axhline(y=1 .- 1 ./ R0 .* (1 .+ log.(R0))) # check peak prevalence
	# plot(t[t .< 11], exp.(t[t .< 11] * (R0 - 1)) / N, linestyle="--", color="gray")
	xlim([-1e-3, 10])
	ylim([-1e-1, 50])
	title("stochastic SIR model dynamics")
	# legend(title=L"$\mathcal{R}_0=$" * @sprintf("%.2f", R0))
	tight_layout()
	savefig("stochastic_sir.pdf", format="pdf")
	gcf()
end

# ╔═╡ c829589a-db9e-11ea-398e-8139432eedc6
begin
	nv = 15
	ne = 2 * nv
	g = erdos_renyi(nv, ne)
	
	# make superspreader
	for k = 1:10
	    add_edge!(g, 1, rand(1:nv))
	end
	
	nodefillc = [RGB(0.000000, 0.716877, 0.554419) for i = 1:nv]
	nodefillc[1] = RGB(1.000000, 0.403746, 0.397903)
	gp = gplot(g, nodefillc=nodefillc)
end

# ╔═╡ Cell order:
# ╠═6da1a832-db99-11ea-0098-9b5872e94312
# ╟─dfc4e82a-db99-11ea-0166-9b9f74fe0999
# ╠═e6413104-db99-11ea-3d6e-cfa983139ddc
# ╟─3758086a-db9a-11ea-1102-bbac22d924a1
# ╠═10f1c0da-db9a-11ea-1f3b-41cfa50c451d
# ╟─9fbb2a38-db9d-11ea-26df-2789729490b3
# ╠═684f7776-db9e-11ea-2dbc-85773151febb
# ╠═8cfa38e8-db9e-11ea-09a2-a366b1fcc94d
# ╟─d30011b2-db9e-11ea-0a84-77f66c5de33b
# ╠═c829589a-db9e-11ea-398e-8139432eedc6
# ╠═51c07944-db9f-11ea-35e3-d9ab257b6bd6
