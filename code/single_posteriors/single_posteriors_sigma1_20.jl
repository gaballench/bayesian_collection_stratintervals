using StratIntervals
using Distributions
using Random
using StatsBase # this one provides stuff like var and iqr
using Turing 
using DataFrames
using CSV
using StatsPlots

### large number behaviours for the Bayesian version

mypalette = palette([:darkblue, :red, :black]) # only way to pick group colours in groupedbar

Random.seed!(1985)

setprogress!(false)

# true params, and all times are before present (minus)
theta2 = 100.0
theta1 = 150.0
lambda = 0.0
ndata = [10, 50, 100, 200, 500]
iters=10000

###########################################################
######### Informative priors with sd=1.0 #########
###########################################################

### Simulation 1: Fixed theta2 (at 100.0) and uniform lambda (0.0)

# priors
theta1_prior = Normal(150.0, 1.0) # somewhat informative prior
theta2_prior = theta2 # fixed to the true value
lambda_prior = lambda # fixed to the true value

# iterate over sample sizes
sampling_theta1 = Array{Chains}(undef, length(ndata))
for i in 1:length(ndata)
    # if sampling from prior
    if ndata[i] == 0
        simdata = rand(ThreeParBeta(theta1, theta2, lambda), 2) # fake data
        sampling_theta1[i] = sample_stratinterval(StratInterval(simdata, theta1_prior, theta2_prior, lambda_prior), iters, NUTS(), true, false)
    end
    # if sampling from posterior
    if ndata[i] > 0
        simdata = rand(ThreeParBeta(theta1, theta2, lambda), ndata[i])
        sampling_theta1[i] = sample_stratinterval(StratInterval(simdata, theta1_prior, theta2_prior, lambda_prior), iters, NUTS(), false, false)
    end
end

# summarise and plot the posterior of theta1

# collect the densities for theta2 and do the violin plot
thetae1_array = map(x -> sampling_theta1[x][:, 1, :], 1:length(sampling_theta1))
violin(thetae1_array, primary=false, color="darkblue")
plot!(xticks = (1:length(ndata), ndata))
hline!([theta1], color="black", legend=false, linestyle=:dash)
# save the fig
savefig("thetae1_sd1.svg")

# calculate the iqr of each estimate of theta2, element-wise
iqrs_theta1 = map(x -> iqr(x), thetae1_array)
bar(iqrs_theta1, legend=false, color="darkblue")
plot!(xticks = (1:length(ndata), ndata))
# save the fig
savefig("thetae1_sd1_iqrs.svg")

### Simulation 2: Fixed theta1 (at 150.0) and uniform lambda (0.0)

# priors
theta1_prior = theta1 # fixed to the true value
theta2_prior = Normal(100.0, 1.0) # somewhat informative prior
lambda_prior = lambda # fixed to the true value

# iterate over sample sizes
sampling_theta2 = Array{Chains}(undef, length(ndata))
for i in 1:length(ndata)
    # if sampling from prior
    if ndata[i] == 0
        simdata = rand(ThreeParBeta(theta1, theta2, lambda), 2) # fake data
        sampling_theta2[i] = sample_stratinterval(StratInterval(simdata, theta1_prior, theta2_prior, lambda_prior), iters, NUTS(), true, false)
    end
    # if sampling from posterior
    if ndata[i] > 0
        simdata = rand(ThreeParBeta(theta1, theta2, lambda), ndata[i])
        #println("FAD-LAD are $(minimum(simdata)) and $(maximum(simdata))")
        sampling_theta2[i] = sample_stratinterval(StratInterval(simdata, theta1_prior, theta2_prior, lambda_prior), iters, NUTS(), false, false)
    end
end

# summarise and plot the posterior of theta1

# collect the densities for theta2 and do the violin plot
thetae2_array = map(x -> sampling_theta2[x][:, 1, :], 1:length(sampling_theta2))
violin(thetae2_array, primary=false, color="darkblue")
plot!(xticks = (1:length(ndata), ndata))
hline!([theta2], color="black", legend=false, linestyle=:dash)
# save the fig
savefig("thetae2_sd1.svg")

# calculate the iqr of each estimate of theta1, element-wise
iqrs_theta2 = map(x -> iqr(x), thetae2_array)
bar(iqrs_theta2, legend=false, color="darkblue")
plot!(xticks = (1:length(ndata), ndata))
# save the fig
savefig("thetae2_sd1_iqrs.svg")

### Simulation 3: Fixed uniform lambda (0.0), both thetae estimated

# priors
theta2_prior = Normal(100.0, 1.0) # somewhat informative prior
theta1_prior = Normal(150.0, 1.0) # somewhat informative prior
lambda_prior = 0.0 # fixed to the true value

# iterate over sample sizes
sampling_theta12 = Array{Chains}(undef, length(ndata))
for i in 1:length(ndata)
    # if sampling from prior
    if ndata[i] == 0
        simdata = rand(ThreeParBeta(theta1, theta2, lambda), 2) # fake data
        sampling_theta12[i] = sample_stratinterval(StratInterval(simdata, theta1_prior, theta2_prior, lambda_prior), iters, NUTS(), true, false)
    end
    # if sampling from posterior
    if ndata[i] > 0
        simdata = rand(ThreeParBeta(theta1, theta2, lambda), ndata[i])
        sampling_theta12[i] = sample_stratinterval(StratInterval(simdata, theta1_prior, theta2_prior, lambda_prior), iters, NUTS(), false, false)
    end
end

# summarise and plot the posterior of coestimated theta1,2

# collect the densities for theta2 and do the violin plot
thetae12_1_array = map(x -> sampling_theta12[x][:, 1, :], 1:length(sampling_theta12))
violin(thetae12_1_array, primary=false, color="darkblue")
plot!(xticks = (1:length(ndata), ndata))
hline!([theta1], color="black", legend=false, linestyle=:dash)
# save the fig
savefig("thetae12_1_sd1.svg")

thetae12_2_array = map(x -> sampling_theta12[x][:, 2, :], 1:length(sampling_theta12))
violin(thetae12_2_array, primary=false, color="darkblue")
plot!(xticks = (1:length(ndata), ndata))
hline!([theta2], color="black", legend=false, linestyle=:dash)
# save the fig
savefig("thetae12_2_sd1.svg")

# calculate the iqr of each estimate of theta1, element-wise
iqrs_theta12_1 = map(x -> iqr(x), thetae12_1_array)
iqrs_theta12_2 = map(x -> iqr(x), thetae12_2_array)

plot(x=1:length(ndata), y=iqrs_theta2)

iqrs_df = DataFrame(theta1=iqrs_theta1, theta2=iqrs_theta2, theta12_1=iqrs_theta12_1, theta12_2=iqrs_theta12_2)

@df iqrs_df plot(1:length(ndata), :theta1, label="\$\\theta_1\$ estimated")
@df iqrs_df plot!(1:length(ndata), :theta2, label="\$\\theta_2\$ estimated")
@df iqrs_df plot!(1:length(ndata), :theta12_1, label="\$\\theta_1\$ when both coestimated")
@df iqrs_df plot!(1:length(ndata), :theta12_2, label="\$\\theta_2\$ when both coestimated")
plot!()
plot!(xticks = (1:length(ndata), ndata))
# save the fig
savefig("thetae_all_sd1_iqrs.svg")

iqrs_df_grouped = DataFrame(group=[repeat(["fixed"], length(iqrs_theta1)); repeat(["coestimated"], length(iqrs_theta12_1))],
                            theta1=[iqrs_theta1; iqrs_theta12_1],
                            theta2=[iqrs_theta2; iqrs_theta12_2])

@df iqrs_df_grouped groupedbar(:theta1, palette=mypalette, group=:group)
plot!(xticks = (1:length(ndata), ndata))
savefig("theta1_sd1_grouped_iqrs.svg")

@df iqrs_df_grouped groupedbar(:theta2, palette=mypalette, group=:group)
plot!(xticks = (1:length(ndata), ndata))
savefig("theta2_sd1_grouped_iqrs.svg")

###########################################################
######### Uninformative priors with sd=20.0         #########
###########################################################

### Simulation 1: Fixed theta2 (at 100.0) and uniform lambda (0.0)

# priors
theta1_prior = Normal(150.0, 20.0) # somewhat informative prior
theta2_prior = theta2 # fixed to the true value
lambda_prior = lambda # fixed to the true value

# iterate over sample sizes
sampling_theta1 = Array{Chains}(undef, length(ndata))
for i in 1:length(ndata)
    # if sampling from prior
    if ndata[i] == 0
        simdata = rand(ThreeParBeta(theta1, theta2, lambda), 2) # fake data
        sampling_theta1[i] = sample_stratinterval(StratInterval(simdata, theta1_prior, theta2_prior, lambda_prior), iters, NUTS(), true, false)
    end
    # if sampling from posterior
    if ndata[i] > 0
        simdata = rand(ThreeParBeta(theta1, theta2, lambda), ndata[i])
        sampling_theta1[i] = sample_stratinterval(StratInterval(simdata, theta1_prior, theta2_prior, lambda_prior), iters, NUTS(), false, false)
    end
end

# summarise and plot the posterior of theta1

# collect the densities for theta2 and do the violin plot
thetae1_array = map(x -> sampling_theta1[x][:, 1, :], 1:length(sampling_theta1))
violin(thetae1_array, primary=false, color="darkblue")
plot!(xticks = (1:length(ndata), ndata))
hline!([theta1], color="black", legend=false, linestyle=:dash)
# save the fig
savefig("thetae1_sd20.svg")

# calculate the iqr of each estimate of theta2, element-wise
iqrs_theta1 = map(x -> iqr(x), thetae1_array)
bar(iqrs_theta1, legend=false, color="darkblue")
plot!(xticks = (1:length(ndata), ndata))
# save the fig
savefig("thetae1_sd20_iqrs.svg")

### Simulation 2: Fixed theta1 (at 150.0) and uniform lambda (0.0)

# priors
theta1_prior = theta1 # fixed to the true value
theta2_prior = Normal(100.0, 20.0) # somewhat informative prior
lambda_prior = lambda # fixed to the true value

# iterate over sample sizes
sampling_theta2 = Array{Chains}(undef, length(ndata))
for i in 1:length(ndata)
    # if sampling from prior
    if ndata[i] == 0
        simdata = rand(ThreeParBeta(theta1, theta2, lambda), 2) # fake data
        sampling_theta2[i] = sample_stratinterval(StratInterval(simdata, theta1_prior, theta2_prior, lambda_prior), iters, NUTS(), true, false)
    end
    # if sampling from posterior
    if ndata[i] > 0
        simdata = rand(ThreeParBeta(theta1, theta2, lambda), ndata[i])
        #println("FAD-LAD are $(minimum(simdata)) and $(maximum(simdata))")
        sampling_theta2[i] = sample_stratinterval(StratInterval(simdata, theta1_prior, theta2_prior, lambda_prior), iters, NUTS(), false, false)
    end
end

# summarise and plot the posterior of theta1

# collect the densities for theta2 and do the violin plot
thetae2_array = map(x -> sampling_theta2[x][:, 1, :], 1:length(sampling_theta2))
violin(thetae2_array, primary=false, color="darkblue")
plot!(xticks = (1:length(ndata), ndata))
hline!([theta2], color="black", legend=false, linestyle=:dash)
# save the fig
savefig("thetae2_sd20.svg")

# calculate the iqr of each estimate of theta1, element-wise
iqrs_theta2 = map(x -> iqr(x), thetae2_array)
bar(iqrs_theta2, legend=false, color="darkblue")
plot!(xticks = (1:length(ndata), ndata))
# save the fig
savefig("thetae2_sd20_iqrs.svg")

### Simulation 3: Fixed uniform lambda (0.0), both thetae estimated

# priors
theta2_prior = Normal(100.0, 20.0) # somewhat informative prior
theta1_prior = Normal(150.0, 20.0) # somewhat informative prior
lambda_prior = 0.0 # fixed to the true value

# iterate over sample sizes
sampling_theta12 = Array{Chains}(undef, length(ndata))
for i in 1:length(ndata)
    # if sampling from prior
    if ndata[i] == 0
        simdata = rand(ThreeParBeta(theta1, theta2, lambda), 2) # fake data
        sampling_theta12[i] = sample_stratinterval(StratInterval(simdata, theta1_prior, theta2_prior, lambda_prior), iters, NUTS(), true, false)
    end
    # if sampling from posterior
    if ndata[i] > 0
        simdata = rand(ThreeParBeta(theta1, theta2, lambda), ndata[i])
        sampling_theta12[i] = sample_stratinterval(StratInterval(simdata, theta1_prior, theta2_prior, lambda_prior), iters, NUTS(), false, false)
    end
end

# summarise and plot the posterior of coestimated theta1,2

# collect the densities for theta2 and do the violin plot
thetae12_1_array = map(x -> sampling_theta12[x][:, 1, :], 1:length(sampling_theta12))
violin(thetae12_1_array, primary=false, color="darkblue")
plot!(xticks = (1:length(ndata), ndata))
hline!([theta1], color="black", legend=false, linestyle=:dash)
# save the fig
savefig("thetae12_1_sd20.svg")

thetae12_2_array = map(x -> sampling_theta12[x][:, 2, :], 1:length(sampling_theta12))
violin(thetae12_2_array, primary=false, color="darkblue")
plot!(xticks = (1:length(ndata), ndata))
hline!([theta2], color="black", legend=false, linestyle=:dash)
# save the fig
savefig("thetae12_2_sd20.svg")

# calculate the iqr of each estimate of theta1, element-wise
iqrs_theta12_1 = map(x -> iqr(x), thetae12_1_array)
iqrs_theta12_2 = map(x -> iqr(x), thetae12_2_array)

plot(x=1:length(ndata), y=iqrs_theta2)

iqrs_df = DataFrame(theta1=iqrs_theta1, theta2=iqrs_theta2, theta12_1=iqrs_theta12_1, theta12_2=iqrs_theta12_2)

@df iqrs_df plot(1:length(ndata), :theta1, label="\$\\theta_1\$ estimated")
@df iqrs_df plot!(1:length(ndata), :theta2, label="\$\\theta_2\$ estimated")
@df iqrs_df plot!(1:length(ndata), :theta12_1, label="\$\\theta_1\$ when both coestimated")
@df iqrs_df plot!(1:length(ndata), :theta12_2, label="\$\\theta_2\$ when both coestimated")
plot!()
plot!(xticks = (1:length(ndata), ndata))
# save the fig
savefig("thetae_all_sd20_iqrs.svg")

iqrs_df_grouped = DataFrame(group=[repeat(["fixed"], length(iqrs_theta1)); repeat(["coestimated"], length(iqrs_theta12_1))],
                            theta1=[iqrs_theta1; iqrs_theta12_1],
                            theta2=[iqrs_theta2; iqrs_theta12_2])

@df iqrs_df_grouped groupedbar(:theta1, palette=mypalette, group=:group)
plot!(xticks = (1:length(ndata), ndata))
savefig("theta1_sd20_grouped_iqrs.svg")

@df iqrs_df_grouped groupedbar(:theta2, palette=mypalette, group=:group)
plot!(xticks = (1:length(ndata), ndata))
savefig("theta2_sd20_grouped_iqrs.svg")
