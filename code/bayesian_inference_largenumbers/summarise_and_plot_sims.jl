using StatsPlots
using DataFrames
using CSV
using CategoricalArrays
using Statistics

# true values
theta1 = 150.0
thet2 = 100.0
lambda = 0.0

# read dir paths and get the idx for each type of simulation
sims = readdir("batchsims", join=true)
theta1_idx = findall(x -> occursin(r"theta1.*.csv", x), sims)
theta2_idx = findall(x -> occursin(r"theta2.*.csv", x), sims)
both_idx = findall(x -> occursin(r"both.*.csv", x), sims)

# initialise the dataframes
theta1_df = CSV.read(sims[theta1_idx][1], DataFrame);
theta2_df = CSV.read(sims[theta2_idx][1], DataFrame);
both_df = CSV.read(sims[both_idx][1], DataFrame);

# read iteratively the theta1 csvs
for i in sims[theta1_idx][2:end]
    theta1_df = vcat(theta1_df, CSV.read(i, DataFrame));
end

# read iteratively the theta2 csvs
for i in sims[theta2_idx][2:end]
    theta2_df = vcat(theta2_df, CSV.read(i, DataFrame));
end

# read iteratively the both csvs
for i in sims[both_idx][2:end]
    both_df = vcat(both_df, CSV.read(i, DataFrame));
end

# sort the dataframes by sigma and then ndata
sort!(theta1_df, [:sigma, :ndata])
sort!(theta2_df, [:sigma, :ndata])
sort!(both_df, [:sigma, :ndata])

# convert to categorical the ndata column so that it can be sorted correctly when using boxplot or density
theta1_df.ndata = categorical(theta1_df.ndata)
theta2_df.ndata = categorical(theta2_df.ndata)
both_df.ndata = categorical(both_df.ndata)

# making StatsPlots restrict the colours so that I can pick the easier
mypalette = palette([:darkblue, :red, :black]) # only way to pick group colours in groupedbar

#######################################################
### plotting the densities for theta1 with sigma=1  ###
#######################################################

# plotting the densities for theta1
@df theta1_df[theta1_df.sigma .== 1.0, :] boxplot(:ndata, :theta_median, primary=false, palette=mypalette, alpha=1.0, legend=false)
hline!([theta1], palette=mypalette, color=3, line=(1, :dash))

# plotting the densities for theta2
@df theta2_df[theta2_df.sigma .== 1.0, :] boxplot(:ndata, :theta_median, primary=false, palette=mypalette, alpha=1.0, legend=false)
hline!([theta2], palette=mypalette, color=3, line=(1, :dash))

# plotting the densities for theta1 when both are coestimated
@df both_df[both_df.sigma .== 1.0, :] boxplot(:ndata, :theta1_median, primary=false, palette=mypalette, alpha=1.0, legend=false)
hline!([theta1], palette=mypalette, color=3, line=(1, :dash))

# plotting the densities for theta2 when both are coestimated
@df both_df[both_df.sigma .== 1.0, :] boxplot(:ndata, :theta2_median, primary=false, palette=mypalette, alpha=1.0, legend=false)
hline!([theta2], palette=mypalette, color=3, line=(1, :dash))

#######################################################
### plotting the densities for theta1 with sigma=20 ###
#######################################################

# plotting the densities for theta1
@df theta1_df[theta1_df.sigma .== 20.0, :] boxplot(:ndata, :theta_median, primary=false, palette=mypalette, alpha=20.0, legend=false)
hline!([theta1], palette=mypalette, color=3, line=(1, :dash))

# plotting the densities for theta2
@df theta2_df[theta2_df.sigma .== 20.0, :] boxplot(:ndata, :theta_median, primary=false, palette=mypalette, alpha=20.0, legend=false)
hline!([theta2], palette=mypalette, color=3, line=(1, :dash))

# plotting the densities for theta1 when both are coestimated
@df both_df[both_df.sigma .== 20.0, :] boxplot(:ndata, :theta1_median, primary=false, palette=mypalette, alpha=20.0, legend=false)
hline!([theta1], palette=mypalette, color=3, line=(1, :dash))

# plotting the densities for theta2 when both are coestimated
@df both_df[both_df.sigma .== 20.0, :] boxplot(:ndata, :theta2_median, primary=false, palette=mypalette, alpha=20.0, legend=false)
hline!([theta2], palette=mypalette, color=3, line=(1, :dash))

#######################################################
### plotting the IQR for theta1 with sigma=1        ###
#######################################################

# plotting the IRQ for theta1
@df theta1_df[theta1_df.sigma .== 1.0, :] boxplot(:ndata, :theta_iqrs, primary=false, palette=mypalette, alpha=20.0, legend=false)

# plotting the IRQ for theta2
@df theta2_df[theta2_df.sigma .== 1.0, :] boxplot(:ndata, :theta_iqrs, primary=false, palette=mypalette, alpha=20.0, legend=false)

# plotting the IRQ for theta1 when both are coestimated
@df both_df[both_df.sigma .== 1.0, :] boxplot(:ndata, :theta1_iqrs, primary=false, palette=mypalette, alpha=20.0, legend=false)

# plotting the IRQ for theta1 when both are coestimated
@df both_df[both_df.sigma .== 1.0, :] boxplot(:ndata, :theta2_iqrs, primary=false, palette=mypalette, alpha=20.0, legend=false)

#######################################################
### plotting the IQR for theta1 with sigma=20       ###
#######################################################

# plotting the IRQ for theta1
@df theta1_df[theta1_df.sigma .== 1.0, :] boxplot(:ndata, :theta_iqrs, primary=false, palette=mypalette, alpha=20.0, legend=false)

# plotting the IRQ for theta2
@df theta2_df[theta2_df.sigma .== 1.0, :] boxplot(:ndata, :theta_iqrs, primary=false, palette=mypalette, alpha=20.0, legend=false)

# plotting the IRQ for theta1 when both are coestimated
@df both_df[both_df.sigma .== 1.0, :] boxplot(:ndata, :theta1_iqrs, primary=false, palette=mypalette, alpha=20.0, legend=false)

# plotting the IRQ for theta1 when both are coestimated
@df both_df[both_df.sigma .== 1.0, :] boxplot(:ndata, :theta2_iqrs, primary=false, palette=mypalette, alpha=20.0, legend=false)

#######################################################
### The coverage for the median is by definition 0.5###
### exactly as seen below, this is nothing new nor  ###
### useful                                          ###
#######################################################

#######################################################
### plotting the coverage for theta1 with sigma=1   ###
#######################################################
"""
# plotting the coverage for theta1
@df combine(groupby(theta1_df[theta1_df.sigma .== 1.0, :], :ndata),
            :theta_coverage => x -> sum(x)/1000) bar(:ndata, :theta_coverage_function, primary=false, palette=mypalette, alpha=20.0, legend=false)

# plotting the coverage for theta2
@df combine(groupby(theta2_df[theta2_df.sigma .== 1.0, :], :ndata),
            :theta_coverage => x -> sum(x)/1000) bar(:ndata, :theta_coverage_function, primary=false, palette=mypalette, alpha=20.0, legend=false)

# plotting the coverage for theta1 when both are coestimated
@df combine(groupby(both_df[both_df.sigma .== 1.0, :], :ndata),
            :theta1_coverage => x -> sum(x)/1000) bar(:ndata, :theta1_coverage_function, primary=false, palette=mypalette, alpha=20.0, legend=false)

# plotting the coverage for theta2 when both are coestimated
@df combine(groupby(both_df[both_df.sigma .== 1.0, :], :ndata),
            :theta2_coverage => x -> sum(x)/1000) bar(:ndata, :theta2_coverage_function, primary=false, palette=mypalette, alpha=20.0, legend=false)

#######################################################
### plotting the coverage for theta1 with sigma=20  ###
#######################################################

# plotting the coverage for theta1
@df combine(groupby(theta1_df[theta1_df.sigma .== 20.0, :], :ndata),
            :theta_coverage => x -> sum(x)/1000) bar(:ndata, :theta_coverage_function, primary=false, palette=mypalette, alpha=20.0, legend=false)

# plotting the coverage for theta2
@df combine(groupby(theta2_df[theta2_df.sigma .== 20.0, :], :ndata),
            :theta_coverage => x -> sum(x)/1000) bar(:ndata, :theta_coverage_function, primary=false, palette=mypalette, alpha=20.0, legend=false)

# plotting the coverage for theta1 when both are coestimated
@df combine(groupby(both_df[both_df.sigma .== 20.0, :], :ndata),
            :theta1_coverage => x -> sum(x)/1000) bar(:ndata, :theta1_coverage_function, primary=false, palette=mypalette, alpha=20.0, legend=false)

# plotting the coverage for theta2 when both are coestimated
@df combine(groupby(both_df[both_df.sigma .== 20.0, :], :ndata),
            :theta2_coverage => x -> sum(x)/1000) bar(:ndata, :theta2_coverage_function, primary=false, palette=mypalette, alpha=20.0, legend=false)

######### examining coverage

#@df theta1_df[theta1_df.sigma .== 1.0, :] bar(:ndata, string.(:theta_coverage), primary=false, palette=mypalette, alpha=20.0, legend=false)
"""

####################################
### Mean squared error is useful ###
####################################

########## mean squared error
@df combine(groupby(theta1_df[theta1_df.sigma .== 1.0, :], :ndata),
            :theta_median => x -> sum((x .- theta1).^2)/1000) bar(:ndata, :theta_median_function, primary=false, palette=mypalette, alpha=20.0, legend=false)

@df combine(groupby(theta2_df[theta2_df.sigma .== 1.0, :], :ndata),
            :theta_median => x -> sum((x .- theta2).^2)/1000) bar(:ndata, :theta_median_function, primary=false, palette=mypalette, alpha=20.0, legend=false)

@df combine(groupby(theta1_df[theta1_df.sigma .== 20.0, :], :ndata),
            :theta_median => x -> sum((x .- theta1).^2)/1000) bar(:ndata, :theta_median_function, primary=false, palette=mypalette, alpha=20.0, legend=false)

@df combine(groupby(theta2_df[theta2_df.sigma .== 20.0, :], :ndata),
            :theta_median => x -> sum((x .- theta2).^2)/1000) bar(:ndata, :theta_median_function, primary=false, palette=mypalette, alpha=20.0, legend=false)
