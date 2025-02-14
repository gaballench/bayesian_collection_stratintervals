using StratIntervals
using Distributions
using Random
using Turing
using StatsPlots

Random.seed!(1986)

# stratinterval param true values, we will vary lambda and then show its effect on the shape of the posterior of thetas
λs = -1.5:0.5:1.5
θ1_true = 550.0
θ2_true = 500.0

# simulation and plotting params
ndata = 75
iters = 100000
xlimits = [460.0, 580.0]
ylimits = [0.0, 0.15]

# priors, quite uninformative
θ1_prior = Normal(θ1_true, 10)
θ2_prior = Normal(θ2_true, 10)

# we will simulate 50 data points from the true ThreeParBeta with varying lambda, and then sample from the posterior
for i in λs
    data_i = rand(ThreeParBeta(θ1_true, θ2_true, i), ndata)
    λ_prior = Normal(i, 2.0)
    stratint_i = StratInterval(data_i, θ1_prior, θ2_prior, λ_prior)
    prior_stratint_i = sample_stratinterval(stratint_i, iters, NUTS(), true, false)
    posterior_stratint_i = sample_stratinterval(stratint_i, iters, NUTS(), false, false)
    histogram(data_i, normalize=:pdf, legend=false, ylim=ylimits, xlim=xlimits,
              bins=range(θ2_true, stop = θ1_true, length = 8))
    plot!(ThreeParBeta(θ1_true, θ2_true, i), color=:darkorchid)
    density!(prior_stratint_i[:, :θ1, :], color=:darkblue)
    density!(posterior_stratint_i[:, :θ1, :], color=:black)
    density!(prior_stratint_i[:, :θ2, :], color=:darkblue)
    density!(posterior_stratint_i[:, :θ2, :], color=:black)
    savefig("plots/lambda_"*string(i)*"_on_thetas.svg")
end

