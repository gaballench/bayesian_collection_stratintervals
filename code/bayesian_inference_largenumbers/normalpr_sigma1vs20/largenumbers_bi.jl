println("Script starting...")

using StratIntervals
using Distributions
using Random
using StatsBase # this one provides stuff like var and iqr
using Turing 
using DataFrames
using CSV

### large number behaviours for the Bayesian version

"""
command line arguments are, in the same order:

sigma sigma in the normal prior, I'm using 1 or 20 for informative vs. uninformative priors
ndata how many data in the strat interval
nsims how many simulations, I am using 1000
mcmciters MCMC generations, I'm using 10000 samples
#seed the seed for reproducing analysis
simulation one of "theta1", "theta2", or "both"

usage

julia largenumbers_bi.jl sigma ndata nsims mcmciters simulation
"""

sigma = parse(Float64, ARGS[1])
ndata = parse(Int64, ARGS[2])
nsims = parse(Int64, ARGS[3])
mcmciters = parse(Int64, ARGS[4])
#seed = parse(Int64, ARGS[5])
sampling = ARGS[5]

# the seed is set itself at random at the beginning and stored in the results filename

seed = rand(1:1000)

Random.seed!(seed)

setprogress!(false)

# true params, and all times are before present (minus)
theta2 = 100.0
theta1 = 150.0
lambda = 0.0

#############################################################
######### We want to repeat the same approach below       ###
######### but many times for assessing coverage, IQR, etc ###
#############################################################

function simulate_posteriors(true_theta1, true_theta2, true_lambda, prior_theta1, prior_theta2, prior_lambda, ndata, nsims, mcmciters, sigma)
    # repeat iterations, for a given fixed sample size
    sampling_theta = Array{Chains}(undef, nsims)
    for i in 1:nsims
        # if sampling from prior
        if ndata == 0
            simdata = rand(ThreeParBeta(true_theta1, true_theta2, true_lambda), 2) # fake data
            sampling_theta[i] = sample_stratinterval(StratInterval(simdata, prior_theta1, prior_theta2, prior_lambda), mcmciters, NUTS(), true, false)
        end
        # if sampling from posterior
        if ndata > 0
            simdata = rand(ThreeParBeta(true_theta1, true_theta2, true_lambda), ndata)
            sampling_theta[i] = sample_stratinterval(StratInterval(simdata, prior_theta1, prior_theta2, prior_lambda), mcmciters, NUTS(), false, false)
        end
    end
    # decide whether we are sampling one or two params and return accordingly
    if sum(map(x -> x isa Float64, (prior_theta1, prior_theta2))) == 1
        # assign to true_param the value of whatever we are sampling, t1 or t2
        # according to the type
        if prior_theta1 isa Float64
            #println("Sampled param is theta2")
            true_param = theta2
        end
        if prior_theta2 isa Float64
            #println("Sampled param is theta1")
            true_param = theta1
        end        
        # summarise the nsimulations posteriors of theta2        
        # collect the densities for theta2
        thetae_array = map(x -> sampling_theta[x][:, 1, :], 1:length(sampling_theta))        
        # calculate posterior medians
        theta_median = map(x -> median(x), thetae_array)        
        # calculate the iqr of each estimate of theta2, element-wise
        theta_iqrs = map(x -> iqr(x), thetae_array)        
        # calculate the coverage
        quartiles1 = map(x -> quantile(x, [0.25])[1], thetae_array)
        quartiles3 = map(x -> quantile(x, [0.75])[1], thetae_array)
        #println(quartiles1)
        #println(quartiles3)
        #println(true_param)
        # calculate the posterior coverage
        theta_coverage = map((x,y) -> x <= true_param <= y, quartiles1, quartiles3)
        return DataFrame(ndata="$ndata", sigma="$sigma", theta_median=theta_median, theta_iqrs=theta_iqrs, quartiles1=quartiles1, quartiles3=quartiles3, theta_coverage=theta_coverage)
    end
    if sum(map(x -> x isa Float64, [prior_theta1, prior_theta2])) == 0        
        # assign to true_param the value of whatever we are sampling, t1 or t2
        # according to the type
        true_theta1 = theta1
        true_theta2 = theta2
        # summarise the nsimulations posteriors of theta2        
        # collect the densities for theta1
        thetae1_array = map(x -> sampling_theta[x][:, 1, :], 1:length(sampling_theta))        
        # collect the densities for theta2
        thetae2_array = map(x -> sampling_theta[x][:, 2, :], 1:length(sampling_theta))        
        # calculate posterior medians
        theta1_median = map(x -> median(x), thetae1_array)        
        theta2_median = map(x -> median(x), thetae2_array)
        # calculate the iqr of each estimate of theta2, element-wise
        theta1_iqrs = map(x -> iqr(x), thetae1_array)        
        theta2_iqrs = map(x -> iqr(x), thetae2_array)
        # calculate the coverage
        quartiles1_1 = map(x -> quantile(x, [0.25])[1], thetae1_array)
        quartiles3_1 = map(x -> quantile(x, [0.75])[1], thetae1_array)
        quartiles1_2 = map(x -> quantile(x, [0.25])[1], thetae2_array)
        quartiles3_2 = map(x -> quantile(x, [0.75])[1], thetae2_array)
        #println(quartiles13)
        #println(true_param)
        # calculate the posterior coverage
        theta1_coverage = map((x,y) -> x <= true_theta1 <= y, quartiles1_1, quartiles3_1)
        theta2_coverage = map((x,y) -> x <= true_theta2 <= y, quartiles1_2, quartiles3_2)
        return DataFrame(ndata="$ndata", sigma="$sigma", theta1_median=theta1_median, theta2_median=theta2_median, theta1_iqrs=theta1_iqrs, theta2_iqrs=theta2_iqrs, quartiles1_1=quartiles1_1, quartiles3_1=quartiles3_1, quartiles1_2=quartiles1_2, quartiles3_2=quartiles3_2, theta1_coverage=theta1_coverage, theta2_coverage=theta2_coverage)
    end
    error("Failed catching number of sampled parameters != 1")
end

###########################################################
######### sampling proper                         #########
###########################################################

println("Sampling arguments:\nsigma=$(sigma), ndata=$(ndata), nsims=$(nsims), mcmciters=$(mcmciters), seed=$(seed), and sampling=$(sampling) were passed as command line arguments")

if sampling == "theta2"
    println("Sampling theta2")
    # simulation 1: sample theta2

    # priors
    prior_theta1 = theta1 # fixed to the true value
    prior_theta2 = Normal(100.0, sigma)
    prior_lambda = lambda # fixed to the true value
    sim_results = simulate_posteriors(theta1, theta2, lambda, prior_theta1, prior_theta2, prior_lambda, ndata, nsims, mcmciters, sigma)
end

if sampling == "theta1"
    println("Sampling theta1")
    # simulation 2: sample theta1

    # priors
    prior_theta1 = Normal(150.0, sigma)
    prior_theta2 = theta2 # fixed to the true value
    prior_lambda = lambda # fixed to the true value
    sim_results = simulate_posteriors(theta1, theta2, lambda, prior_theta1, prior_theta2, prior_lambda, ndata, nsims, mcmciters, sigma)
end

if sampling == "both"
    println("Sampling both theta1 and theta2")
    # simulation 2: sample theta1

    # priors
    prior_theta1 = Normal(150.0, sigma)
    prior_theta2 = Normal(100.0, sigma)
    prior_lambda = 0.0 # fixed to the true value
    sim_results = simulate_posteriors(theta1, theta2, lambda, prior_theta1, prior_theta2, prior_lambda, ndata, nsims, mcmciters, sigma)
end

CSV.write("$(sampling)_$(ndata)_$(sigma)_$(seed).csv", sim_results, header=true)

println("Bye bye!")

exit()