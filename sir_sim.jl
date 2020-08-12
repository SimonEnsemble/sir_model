### A Pluto.jl notebook ###
# v0.11.4

using Markdown
using InteractiveUtils

# ╔═╡ 3fc751e4-db87-11ea-1eab-fb78edd0a299
md"# numerical solution to the SIR model

define the basic reproduction number"

# ╔═╡ 2e045442-db88-11ea-0141-d7d72b844d72
md"define the right-hand-side of the ODE, viewed as:

$\frac{d\mathbf{u}}{dt}=f(\mathbf{u}, p, t)$

where $\mathbf{u}=\mathbf{u}(t):=[[S](t), [I](t), [R](t)]$."

# ╔═╡ 67797b4e-db88-11ea-372e-d92979ccd051
md"initial conditions"

# ╔═╡ 84c3a832-db88-11ea-115c-2ff4dd767dfc
md"numerically solve the ODE over a time span"

# ╔═╡ ac2addc8-db88-11ea-0613-571f02ebc56e
md"### visualize solution

first, a traditional visualization of the curves $[S](t)$, $[I](t)$, and $[R](t)$"

# ╔═╡ 2a7584fa-db8c-11ea-2a9c-f3f7b0b7fda6
md"second, a viz that emphasizes $[S](t)+[I](t)+[R](t)=1$, $\forall t$"

# ╔═╡ b78e7f2c-db8c-11ea-2bd6-4b7491b14d65
md"finally, a viz of the dynamics in the phase plane."

# ╔═╡ 1e039f84-db8e-11ea-1e60-b74a403dc896
md"one more thing. let's look at the exponential growth (of $[I](t)$) phase"

# ╔═╡ 2e1599f6-db96-11ea-2a83-e3e8a3213e4e
md"## peak prevalence, final size as a function of $\mathcal{R}_0$"

# ╔═╡ 3dd997b8-db97-11ea-3fdf-33ba59dcdce6
md"## herd immunity"

# ╔═╡ d54ba48e-dc56-11ea-0ac4-8d213422f65f
md"### find time to peak infection"

# ╔═╡ 988e2cdc-db85-11ea-006e-398146591b74
begin
	using DifferentialEquations, PyPlot, Roots, Printf
	PyPlot.matplotlib.style.use("grandbudapest.mplstyle")
	PyPlot.matplotlib.font_manager.fontManager.addfont("OpenSans-Regular.ttf")
	
	rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
	rcParams["mathtext.fontset"] = "dejavusans"
end

# ╔═╡ 263e746a-db88-11ea-0a72-5d23fe7d94be
R₀ = 2.0

# ╔═╡ 4a6df40e-db87-11ea-29a8-6964c5ac203a
function update_f!(f, u, p, t)
    # for clarity, unpack vector u
    s = u[1]
    i = u[2]
    r = u[3]
	# for clarity, upack p[1] as R₀
	R₀ = p[1]
    # update f
    f[1] = -R₀ * s * i
    f[2] = i * (R₀ * s - 1.0)
    f[3] = i
end

# ╔═╡ 70aebf96-db88-11ea-3c06-53267db7cc80
begin
	ϵ = 10^-5             # fraction of population initially infectious
	u0 = [1-ϵ; ϵ; 0.0]    # [S(0), I(0), R(0)]
end

# ╔═╡ 801542c8-db88-11ea-1087-d7da230c7a6f
begin
	# define the ODE problem
	time_span = (0.0, 25.0)
	prob = ODEProblem(update_f!, u0, time_span, [R₀])
	
	# numerically solve ODE
	sol = solve(prob) # now sol(10.0) gives numerical approx. to sol'n at t=10
end

# ╔═╡ b5dd3ccc-db87-11ea-1f18-8571cc964420
begin
	run_checks = false
	sir_to_color = Dict("S" => "C3", "I" => "C5", "R" => "C0")
	
	t = range(0.0, time_span[2], length=500)
	
	s = [sol.(t_j)[1] for t_j in t]
	i = [sol.(t_j)[2] for t_j in t]
	r = [sol.(t_j)[3] for t_j in t]
	
	figure()
	plot(t, s, color=sir_to_color["S"], clip_on=false)
	plot(t, i, color=sir_to_color["I"], clip_on=false, zorder=100)
	plot(t, r, color=sir_to_color["R"], clip_on=false)
	xlabel(L"non-dimensional time, $\gamma t$")
	ylabel("fraction of population")
	
	dy = 0.04
	x_pos = 21.0
	text(x_pos, s[end] + dy, L"$[$S$](t)$", color=sir_to_color["S"])
	text(x_pos, i[end] + dy, L"$[$I$](t)$", color=sir_to_color["I"])
	text(x_pos, r[end] + dy, L"$[$R$](t)$", color=sir_to_color["R"])
	s∞ = fzero(x -> log(x) - R₀ * (x - 1), 0.1)
	if run_checks
		# check final size formula
		axhline(y=s∞, linestyle="--", color="gray")
		# check peak prevalence
		axhline(y=1 .- 1 ./ R₀ .* (1 .+ log.(R₀)), linestyle="--", color="k")
		# check exponential dependence
		plot(t[t .< 11], ϵ * exp.(t[t .< 11] * (R₀ - 1)), linestyle="--", color="r")
	end
	xlim([-1e-3, time_span[2]])
	ylim([-1e-3, 1+1e-3])
	
	# info
	bbox_props = Dict(:boxstyle=>"round", :fc=>(1.0, 0.98, 0.98), :ec=>"0.5")
		#:ec=>"0.5")#, :alpha=0.9)
	sim_settings = L"$\mathcal{R}_0=2$" * "\n" * L"$[$I$]_0=10^{-5}$"
	@assert u0[2] ≈ 10^(-5)
	@assert R₀ ≈ 2.0
	text(2.5, 0.5, 
		sim_settings, bbox=bbox_props,
	    va="center")
	title("SIR model dynamics")
	# legend(title=L"$\mathcal{R}_0=$" * @sprintf("%.2f", R0))
	tight_layout()
	savefig("sir_dynamics.pdf", format="pdf")
	gcf()
end

# ╔═╡ 2ca716b0-db89-11ea-0e4b-91cb095da5bf
begin
	figure()
	fill_between(t, i+r, s+i+r, color=sir_to_color["S"], label=L"$[$S$](t)$")
	fill_between(t, i  , 0.0  , color=sir_to_color["I"], label=L"$[$I$](t)$")
	fill_between(t, i  , i+r  , color=sir_to_color["R"], label=L"$[$R$](t)$")
	xlim([0, 25])
	ylim([0, 1])
	legend(loc="upper left")
	xlabel(L"non-dimensional time, $\gamma t$")
	ylabel("fraction of population")
	text(2.5, 0.5, 
		sim_settings, #bbox=bbox_props,
		va="center")
	title("SIR model dynamics")
	tight_layout()
	savefig("sir_dynamics_sum_1.pdf", format="pdf")
	gcf()
end

# ╔═╡ c03dd014-db8c-11ea-39bf-9d19c901c208
begin
	cmap = PyPlot.matplotlib.cm.get_cmap("viridis") # cmap(x) for x ∈ [0, 1]
	
	# map a time to a color
	function time_to_color(t::Float64, t_max::Float64)
		if t > t_max
			return cmap(1.0)
		else
			return cmap(t / t_max)
		end
	end
	
	t_max = 15.0
	
	figure(figsize=(8, 7))
	plot(s, i, color="C4")
	for k = 1:length(s)-1
	    plot(s[k:k+1], -s[k:k+1] .+ log.(s[k:k+1]) / R₀ .+ 1, 
	        color=time_to_color((t[k] + t[k+1]) / 2, t_max)
	    )
	end
	
	if run_checks
		# check phase plane solution
		plot(s, 1 .- s .+ 1 / R₀ * log.(s / u0[1]))
	end
	
	# initial conditions, peak prevalence, and final conditions
	xys = [(u0[1],  u0[2]), 
		   (1 / R₀, 1 - 1 / R₀ .* (1 + log(R₀ * u0[1]))),
		   (s∞, 0.0)
	]
	labels = [L"$([$S$]_0, [$I$]_0)$", 
			  L"$(\mathcal{R}_0^{-1}, \max_t [$I$](t))$",
		      L"$([$S$]_\infty, 0)$"
	]
	dxs = [-0.14, -0.0, -0.1]
	dys = [0.175, 0.15, 0.15]
	legend_not_annotate = false
	for k = 1:3
		plot(xys[k][1], xys[k][2], marker="X", markersize=7, clip_on=false, 
	    	color="k", label=labels[k], linestyle="None")
		if ! legend_not_annotate
			annotate(labels[k],
				xy=xys[k], xycoords="data",
				xytext=(xys[k][1] + dxs[k], xys[k][2] + dys[k]), textcoords="data",
				arrowprops=Dict(:arrowstyle => "->", #linestyle="dashed",
								:color => "k"
								),
				ha="center",
				bbox=Dict(:boxstyle=>"round", :fc=>(1.0, 0.98, 0.98), :ec=>"0.5",
					:alpha=>0.8),
				zorder=1241
			 )
		end
	end
	
	xticks([i * 0.25 for i = 0:5])
	yticks([i * 0.25 for i = 0:5])
	
	if legend_not_annotate
		text(0.65, 0.5, 
			sim_settings, #bbox=bbox_props,
			va="center")
	else
		text(0.65, 0.85, 
			sim_settings, #bbox=bbox_props,
			va="center")
	end
	
	if legend_not_annotate
		legend(loc="upper right", numpoints=1, prop=Dict(:size => 14))
	end

	
	# color invalid region
	fill_between([0.0, 1.0], [1, 0], [1, 1], color="0.8")
	xlabel(L"$[$S$](t)$")
	ylabel(L"$[$I$](t)$")
	title("phase plane")
	sm = PyPlot.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0.0, vmax=t_max))
	colorbar(sm, label=L"non-dimensional time, $\gamma t$", extend="max")
	
	# set equal axes
	plt.gca().set_aspect("equal")
	xlim([-1e-2, 1])
	ylim([-1e-2, 1])
	
	tight_layout()
	savefig("phase_plane.pdf", format="pdf", bbox_inches="tight")
	gcf()
end

# ╔═╡ fb3e3094-db93-11ea-2ea0-ad854ade3702
begin
	figure()
	plot(t, i, color=sir_to_color["I"], clip_on=false)#, label=L"$[I](t)$")
	plot(t[t .< 11], exp.(t[t .< 11] * (R₀ - 1)) * u0[2], 
		linestyle="--", color="gray")
	text(x_pos, i[end] + dy/3, L"$[$I$](t)$", color=sir_to_color["I"])
	xlim([0, time_span[2]])
	xlabel(L"non-dimensional time, $\gamma t$")
	ylabel(L"$[$I$](t)$")
	ylim([-1e-3, 0.2])
	title("early exponential growth")
	text(1.5, 0.1, 
		sim_settings, bbox=bbox_props,
		va="center")
	tight_layout()
	savefig("exponential.pdf", format="pdf")
	gcf()
end

# ╔═╡ 9b8aba12-db95-11ea-0680-3b2a75a4af3b
begin
	figure()
	plot(t, i, color=sir_to_color["I"])#, label=L"$[I](t)$")
	plot(t[t .< 11], exp.(t[t .< 11] * (R₀ - 1)) * u0[2], 
		linestyle="--", color="gray")
	text(x_pos, i[end] + dy/20, L"$[$I$](t)$", color=sir_to_color["I"])
	xlim([-1e-3, time_span[2]+1e-2])
	yscale("log")
	xlabel(L"non-dimensional time, $\gamma t$")
	ylabel(L"$[$I$](t)$")
	ylim(ymax=0.2)
	title("early exponential growth")
	tight_layout()
	text(5.5, 1e-4, 
		sim_settings, bbox=bbox_props,
		va="center")
	savefig("log_I.pdf", format="pdf")
	gcf()
end

# ╔═╡ 3fab154c-db96-11ea-0ff6-d5252766c2b3
begin
	R₀s = range(1.0, 5.0, length=100)
	figsize = (4.5, 4.0) # used for three plots below.
end

# ╔═╡ 4eb12950-db96-11ea-0797-39eb74fe7087
begin
	max_i = 1 .- 1 ./ R₀s .* (1 .+ log.(R₀s * u0[1]))
	
	figure(figsize=figsize)
	plot(R₀s, max_i, color="C4")
	xlabel(L"$\mathcal{R}_0$")
	ylabel(L"$\max_t$ $[$I$](t)$")
	xlim([1-1e-3, 5+1e-3])
	ylim(ymin=-1e-3)
	title("peak prevalence")
	tight_layout()
	
	text(3.5, 0.1, 
		LaTeXString(split(sim_settings)[2]), bbox=bbox_props,
		va="center")
	# xlim([1.0, maximum(R0)])
	savefig("peak_i.pdf", format="pdf")
	gcf()
end

# ╔═╡ de8c5138-db98-11ea-3942-19e5356d0bcd
LaTeXString(split(sim_settings)[2])

# ╔═╡ 7e6e40b0-db96-11ea-195c-837936eaec5a
begin
	s∞s = similar(R₀s)
	for i = 1:length(R₀s)
	    f(x) = log(x / u0[1]) - R₀s[i] * (x - 1)
	    if R₀s[i] <= 2.0
	        s∞s[i] = fzero(f, 0.1)
	    else
	        s∞s[i] = fzero(f, 0.0, 0.5)
	    end
	end
	
	figure(figsize=figsize)
	plot(R₀s, 1.0 .- s∞s, color="C1", clip_on=false)
	xlabel(L"$\mathcal{R}_0$")
	ylabel(L"$[$R$]_\infty=1-[$S$]_\infty$")
	xlim([1-1e-3, 5+1e-3])
	ylim([-1e-3, 1+1e-3])
	title("final size")
	text(3.5, 0.2, 
		LaTeXString(split(sim_settings)[2]), bbox=bbox_props,
		va="center")
	tight_layout()
	# xlim([1.0, maximum(R0)])
	savefig("rifnty.pdf", format="pdf")
	gcf()
end

# ╔═╡ 43b2803c-db97-11ea-0b9c-7f5d0fdff33d
begin
	v = 1 .- 1 ./ R₀s
	
	figure(figsize=figsize)
	fill_between(R₀s, v, ones(length(v)), color="C6", alpha=0.4, 
		label=L"$v>1-\mathcal{R}_0^{-1}$")
	plot(R₀s, v, color="C6", clip_on=false)
	plot(R₀s, 1 .- s∞s, color="C1", label=L"$[$R$]_\infty=1-[$S$]_\infty$",
		clip_on=false)
	xlabel(L"$\mathcal{R}_0$")
	ylabel("fraction of population")
	xlim([1-1e-3, 5+1e-3])
	ylim([-1e-3, 1+1e-3])
	legend(loc="lower right")
	title("herd immunity")
	tight_layout()
	savefig("herd_immunity.pdf", format="pdf")
	gcf()
end

# ╔═╡ e0ec7372-dc56-11ea-007e-b9e622a0af95
begin
	# define the ODE problem
	new_time_span = (0.0, 65.0)
	
	R₀s_new = range(1.17, 5.0, length=150)
	time_to_peak = similar(R₀s_new)
	
	for (i, R₀) in enumerate(R₀s_new)
		prob = ODEProblem(update_f!, u0, new_time_span, [R₀])
		# numerically solve ODE
		sol = solve(prob)
	
		# find time at peak by grid search
		ts = range(0.0, R₀ < 2 ? new_time_span[2] : 15.0, length=500) # t-grid
		new_i = [sol.(t_j)[2] for t_j in ts]
		time_to_peak[i] = ts[argmax(new_i)]
		
		# to make sure we are finding the time to peak, get max I(t)
		Imax = 1 - 1 / R₀ .* (1 + log(R₀ * u0[1]))
		if i != 1
			continue
		end
		@assert maximum(new_i) <= Imax
		@assert maximum(new_i) >= 0.999 * Imax
	end
	
	ymax = 50.0
	figure()
	ylabel(L"$\alpha$" * " argmax " * L"$[$I$](t)$")
	xlabel(L"$\mathcal{R}_0$")
	plot(R₀s_new, time_to_peak, color="C7")
	title("time to peak prevalence")
	xlim([-1e-3, 5+1e-3])
	ylim([-1e-1, ymax])
	hlines(0.0, 0.0, 1.0, clip_on=false, color="C7")
	vlines(1.0, 0.0, ymax, linestyle="--", color="C6", linewidth=2)
	text(4.0, 30.0, 
		LaTeXString(split(sim_settings)[2]), bbox=bbox_props,
		va="center", ha="center")
	tight_layout()
	savefig("time_to_first_peak.pdf", format="pdf")
	gcf()
end

# ╔═╡ Cell order:
# ╠═988e2cdc-db85-11ea-006e-398146591b74
# ╟─3fc751e4-db87-11ea-1eab-fb78edd0a299
# ╠═263e746a-db88-11ea-0a72-5d23fe7d94be
# ╟─2e045442-db88-11ea-0141-d7d72b844d72
# ╠═4a6df40e-db87-11ea-29a8-6964c5ac203a
# ╟─67797b4e-db88-11ea-372e-d92979ccd051
# ╠═70aebf96-db88-11ea-3c06-53267db7cc80
# ╟─84c3a832-db88-11ea-115c-2ff4dd767dfc
# ╠═801542c8-db88-11ea-1087-d7da230c7a6f
# ╟─ac2addc8-db88-11ea-0613-571f02ebc56e
# ╠═b5dd3ccc-db87-11ea-1f18-8571cc964420
# ╟─2a7584fa-db8c-11ea-2a9c-f3f7b0b7fda6
# ╠═2ca716b0-db89-11ea-0e4b-91cb095da5bf
# ╟─b78e7f2c-db8c-11ea-2bd6-4b7491b14d65
# ╠═c03dd014-db8c-11ea-39bf-9d19c901c208
# ╟─1e039f84-db8e-11ea-1e60-b74a403dc896
# ╠═fb3e3094-db93-11ea-2ea0-ad854ade3702
# ╠═9b8aba12-db95-11ea-0680-3b2a75a4af3b
# ╟─2e1599f6-db96-11ea-2a83-e3e8a3213e4e
# ╠═3fab154c-db96-11ea-0ff6-d5252766c2b3
# ╠═4eb12950-db96-11ea-0797-39eb74fe7087
# ╠═de8c5138-db98-11ea-3942-19e5356d0bcd
# ╠═7e6e40b0-db96-11ea-195c-837936eaec5a
# ╟─3dd997b8-db97-11ea-3fdf-33ba59dcdce6
# ╠═43b2803c-db97-11ea-0b9c-7f5d0fdff33d
# ╟─d54ba48e-dc56-11ea-0ac4-8d213422f65f
# ╠═e0ec7372-dc56-11ea-007e-b9e622a0af95
