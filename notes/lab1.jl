### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ ad0158d2-6934-11f0-0595-d51407442e58
using Pkg;Pkg.activate("/mnt/c/DATA/PopGen25/PopGen25")

# ╔═╡ a2e7d3d4-8145-4c5e-8dfe-6c9e272ef66e
using Plots, PlutoUI, Distributions, Measures, LaTeXStrings, PlutoExtras

# ╔═╡ e5b42b7c-147f-45b8-b4f6-148d2efd9df6
initialize_eqref()

# ╔═╡ 8db1d84c-3116-4da5-9cf6-2bd7093c9225
md"""
# One great reason to like Julia: Pluto

While the goal is to get you all to have scripts in a language you are comfortable working in, I hope you'll at least take a glance at `Julia`. One of the reasons it works well is it has very nice interactive abilities in the form of its `Pluto` package. We'll take a moment to set-up and get this notebook running for yourself, but after this lecture we'll work more free-form and let you use whatever language you desire.

## Let's expand our Wright-Fisher simulations

Today's goal is to take what we have learned in the lecture and discussion, and create Wright-Fisher simulators that we can use as toy models for the rest of the semester. These will necessarily be over-simplified, but they are here just to give you a chance to test your intuition. In Pluto, each code block can be hidden/revealed, but we can also work to include some interactivity to make it work well.

To make this scalable, we'll write two functions. One to randomly mate, one to loop over the Wright-Fisher process for $t$ generations.
"""

# ╔═╡ 7afaaf98-2de8-4797-8640-ed16b1aaef82
md"""
In the cell above, we define how the WrightFisher process should work. We track the allele frequency each generation, starting from $g$ = 1 to $t$. Each generation, the next generation's frequency is just the result of `randomMating` on the previous generation's allele frequency.
"""

# ╔═╡ 73433c15-1a58-4a71-91bc-9377ad739c0c
function randomMating(N,p)
	return(only(rand(Binomial(N,p),1)/N))
end

# ╔═╡ 25747531-b333-4fd8-9b12-0feb528d877b
function WrightFisher(N,p,t,μ,s)
	g=1
    freqs=Float64[]
    append!(freqs,p)
    while g<t
        next = randomMating(N,freqs[g])
        append!(freqs,next)
        g += 1
    end
    return(freqs)
end

# ╔═╡ 753e4dd1-9b35-480b-8585-b0910e946914
md"""
And here we run the actual random sampling function. Recall from last time that it's really just a `Binomial` sample with probability $p$ and $N$ new individuals each generation, but we need the proportion of successes, so we divide the result by $N$.

Let's see if this works! One nice thing about interactive notebooks, is as we change initial parameters the code is quickly re-run. In fact, we can also add fun things like sliders for initial frequency:

_p_ = $(@bind p Slider(0:0.01:1;default=0.5))
"""

# ╔═╡ 7296bd84-3500-435e-9664-4a233bdabd82
md"""
An odd, but in some ways nice, thing about Pluto is that the bindings of the parameters can change in later cells (the whole thing is evaluated for inter-dependencies, and cells are re-run in the necessary order).

_t_ = $(@bind t NumberField(10:10000;default=100))

_N_ = $(@bind N NumberField(10:1000000;default=1000))

"""

# ╔═╡ 365ddd14-8bd4-4e44-a2a0-767b6cb24982
md"""
## Going beyond starndard WrightFisher

Today we talked about how the standard assumptions of the WF model can often be broken. One of the most important ways WF is wrong is in not accounting for `mutation`. This is a fairly simple thing to include in our broad model, since we can make mutation quite general.

### Mutation

Mutation can be thought of as a random event during the production of gametes (so before formation of the next generation). A gamete which _should_ carry the `a` allele will instead carry an `A`, or vice-versa. In the most general model, the probability of both events is equal and is given as $\mu$. To add it to our code, we need to define one new functions, and modify our WrightFisher function as well.

We'll save computational time by simplifying mutation - recall that we produce _infinite_ gametes, so we don't need to do random sampling for them. The allele frequencies can change deterministically.

Try to fill out the code below yourself, there should be helpful hints that pop up.
"""

# ╔═╡ ad025dd9-41a2-4ae3-8e9d-8f7dd04ca626
function mutation(p,μ)
	next = p
	return(next)
end

# ╔═╡ c9e74853-606c-4331-a768-6cc3ca323b68
md"""

We can then define our mutation rate for simulations to be:

_μ_ = $(@bind μ NumberField(0:0.0001:0.5;default=0.000))

Just so you don't have to scroll back up, here's the plot of the simulations again:
"""

# ╔═╡ c1fa92f7-b86b-44bb-a2c5-d318365dd2e8
md"""
### Selection

Selection can be implemented in a similar fashion: Let's suppose the number of gametes is still infinite, but the proportion of gametes that each individual contributes now depends on its genotype. The fitness of `aa` will be considered the reference fitness, so will always be set to 1. We can then specify the relative fitness of `Aa` and `AA` individuals, and we know that the frequency of the `A` allele among the gametes will be:  
"""

# ╔═╡ 5908f3ea-77ad-44f4-a985-22fb752b2b93
texeq"p` = \frac{(p_{AA}W_{AA}+1/2p_{Aa}W{Aa})}{(p_{AA}W_{AA}+p_{Aa}W_{Aa}+p_{aa}W_{aa})} = \frac{p_{AA}w_{AA}+1/2p_{Aa}w_{Aa}}{1+p_{AA}(w_{AA}-1)+p_{Aa}(w_{Aa}-1)} \label{sel}"

# ╔═╡ aa82cb4b-e23a-4748-9d1f-5153535e0978
md"""
Let's make a function that can modify the allele frequencies due to selection. Just fill out the below with the quantity in Eq $(eqref("sel")).
"""

# ╔═╡ 98d5e5a3-f7fc-4d7a-9cea-6dbf9d23d074
function sel(wAA,wAa,pAA,pAa) 
	next = (pAA*wAA .+ 0.5 .* pAa .* wAa)./(1 .+ pAA .* (wAA .- 1) .- pAa .* (wAa .-1))
	return(next)
end

# ╔═╡ b95d6076-aaf0-4df7-9b32-d98aaebdb14b
md"""
We next need to modify the WrightFisher equation to include selection.

There's different ways of approaching this, but we'll go for a relatively simple one that will only require a single parameter by focusing on directional selection.

Specifically, we'll think of the case where `AA` has twice the effect of `Aa`, and the fitness of `Aa` is _1+s_.

_s_ = $(@bind s NumberField(-0.5:0.01:0.5;default=0.0))

"""

# ╔═╡ 0a6aa8d1-8a6f-41c3-a308-d94992029056
begin
	sims = [WrightFisher(N,p,t,μ,s) for _ in 1:1000]
	p1 = plot(sims,leg=false, linecolor=:black,linealpha=0.2,ylims=(0,1),
		xlabel="Generations",ylabel=L"p_t",margins=5mm)
	plot!(p1,mean.(eachcol(sims)),linecolor=:red,linewidth=2)
end

# ╔═╡ f2938983-93e3-4b88-a6c9-f95b2b26f212
p1

# ╔═╡ 23ea68e3-0a30-4ac8-8380-962463ee8036
	function dirsel(s,p)
		return(sel(1.0+2.0*s,1.0+s,p*p,2*p*(1-p)))
	end

# ╔═╡ a4dba1da-c1a5-4592-b4bd-daa42bdb16d3
begin
	p2 = bar(["aa","Aa","AA"],[1.0,1.0+s,1.0+2s],leg=false,y="Relative Fitness",margins=5mm,ylim=(0.8,1.2))
	plot(p2,p1)
end

# ╔═╡ 6c458fe4-2e11-46c0-b143-4255609f6108
md"""
### Putting it together

One nice thing about the way that `Julia` evaluates functions is that we can do some nice generalizations over parameter space more broadly. For instance, you might want to know what the _equilibrium_ allele frequency is in the long run. 

There are a few approaches here: you could use some modeling package (`SciML`, for instance) to ask what the long term allele frequency is. Because we are including a _stochastic_ element, in the form of drift, those solutions will be unattainable for a bunch of the parameters. Here we'll plot the average result after 50 simulations in their final generation.

_N_ = $(@bind N2 NumberField(10:1000000;default=1000))

_p_ = $(@bind pt Slider(0:0.01:1;default=0.5) )

_t_ = $(@bind tt Slider(10:1:1000;default =100) )
"""

# ╔═╡ fe18db09-2d67-4aa1-ae3e-5e2aee293978
begin
	s_range = 10 .^ range(-6,-1, length=25)
	mu_range = 10 .^ range(-6,-1, length=25)
	function lastWF(μ,s;reps=10)
		mean([only(last(WrightFisher(N2,pt,tt,μ,s))) for _ in 1:reps])
	end
	Z = @. lastWF(mu_range',s_range;reps=50)
	contour(mu_range,s_range,Z,fill=true,title="Frequency after $tt generations",xlabel=L"\mu",ylabel=L"s",scale=:log10,clims=(0,1))
end

# ╔═╡ Cell order:
# ╠═ad0158d2-6934-11f0-0595-d51407442e58
# ╠═e5b42b7c-147f-45b8-b4f6-148d2efd9df6
# ╟─8db1d84c-3116-4da5-9cf6-2bd7093c9225
# ╠═25747531-b333-4fd8-9b12-0feb528d877b
# ╟─7afaaf98-2de8-4797-8640-ed16b1aaef82
# ╠═73433c15-1a58-4a71-91bc-9377ad739c0c
# ╟─753e4dd1-9b35-480b-8585-b0910e946914
# ╠═0a6aa8d1-8a6f-41c3-a308-d94992029056
# ╟─7296bd84-3500-435e-9664-4a233bdabd82
# ╟─365ddd14-8bd4-4e44-a2a0-767b6cb24982
# ╠═ad025dd9-41a2-4ae3-8e9d-8f7dd04ca626
# ╟─c9e74853-606c-4331-a768-6cc3ca323b68
# ╟─f2938983-93e3-4b88-a6c9-f95b2b26f212
# ╟─c1fa92f7-b86b-44bb-a2c5-d318365dd2e8
# ╟─5908f3ea-77ad-44f4-a985-22fb752b2b93
# ╟─aa82cb4b-e23a-4748-9d1f-5153535e0978
# ╠═98d5e5a3-f7fc-4d7a-9cea-6dbf9d23d074
# ╟─b95d6076-aaf0-4df7-9b32-d98aaebdb14b
# ╠═23ea68e3-0a30-4ac8-8380-962463ee8036
# ╟─a4dba1da-c1a5-4592-b4bd-daa42bdb16d3
# ╟─6c458fe4-2e11-46c0-b143-4255609f6108
# ╠═fe18db09-2d67-4aa1-ae3e-5e2aee293978
# ╠═a2e7d3d4-8145-4c5e-8dfe-6c9e272ef66e
