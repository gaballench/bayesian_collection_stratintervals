using StratIntervals
using Turing
using Distributions
using CSV
using DataFrames
using Random
using StatsPlots

"""
command line arguments are, in the same order:

clean_uniques a 0-1 value to be parsed as boolean for activating unique values or including all duplicates
theta1_prior the actual string _between quotation marks_ of Julia code with the distribution specification, e.g. "Normal(5.0,1.0)"

usage

julia origin_barracudas.jl clean_uniques theta1_prior

example

julia origin_barracudas.jl 1 "Normal(100.0, 2.0)"
"""

# command-line argument parsing
clean_uniques = parse(Bool, ARGS[1])
theta1_prior = eval(Meta.parse(ARGS[2]))
#seed_n = parse(Int64, ARGS[3])

# tell us what you are doing this time
println("Starting sampling with the following arguments:\n"*ARGS[1]*", "*ARGS[2])

# set plot titles
main_title = "clean="*string(clean_uniques)*", prior="*string(ARGS[2])

# set xlim and ylim for all plots for easily composing later
xlim_interval = [0, 200]
ylim_interval = [0, 0.075]

# seed set automatically but saved to filename
seed_n = rand(1:1000)

# concatenate cleaning bool, prior, and seed as part of the filename
filename = ARGS[1]*"_"*replace(ARGS[2], "," => "")*"_"*string(seed_n)

Random.seed!(seed_n)

# do not print out the progress bar as it will be parallel exec
setprogress!(false)


### read the data
dataset = CSV.read("Sphyraenidae.csv", DataFrame)

data = dataset.midpoint_ma
if clean_uniques
    unique!(data)
end

### priors

#theta1_prior already specified through command line args
theta2_prior = 0.0
lambda_prior = Normal(0.0, 2.0)

### Using all the occurrences as they are

#histogram(data, bins=range(minimum(data), stop = maximum(data), length = 18), normalize=true, plot_title=main_title, legend=false, xlim=xlim_interval, ylim=ylim_interval)
#plot!(theta1_prior)
#savefig("figures/occs_"*filename*".svg")

# stratinterval object specification and MCMC sampling
barracuda_interval = StratInterval(data, theta1_prior, theta2_prior, lambda_prior)

post_sample = sample_stratinterval(barracuda_interval, 10000, NUTS(), false, false)

plot(post_sample, plot_title=main_title)
savefig("figures/post_sampling_"*filename*".svg")

# plot everything, data, prior, and posterior
histogram(data, bins=range(minimum(data), stop = maximum(data), length = 18), normalize=true, plot_title=main_title, xlim=xlim_interval, ylim=ylim_interval, label="Data")
plot!(theta1_prior, plot_title=main_title, label="Prior", linewidth=3)
density!(post_sample[:,"θ1",:], label="Posterior", linewidth=3)
savefig("figures/prior_vs_posterior_"*filename*".svg")

# save the summary statistics
CSV.write("summary/posterior_"*filename*"stats.txt", DataFrame(describe(post_sample)[1]))
CSV.write("summary/posterior_"*filename*"quantiles.txt", DataFrame(describe(post_sample)[2]))
