using StratIntervals
using CSV
using DataFrames
using StatsPlots
using Distributions
using Turing
using LaTeXStrings
using Random
# for calculating probs on conflation PDFs
using Interpolations
using KernelDensity
using QuadGK
using Optim

Random.seed!(2005)

dataset = CSV.read("../data/cerrejon_palynomorphs.tsv", DataFrame)

#The following cases will be used for fixing theta2, decreasing in sample size
fixt2 = dataset[dataset.fix .== "t2", :species]

#The following cases will be used for fixing none, decreasing in sample size
fixnone = dataset[dataset.fix .== "none", :species]

"""
The best example of each fixing class
t2: Psilatricolporites_sp
none: Arecipites_regio
"""

# Psilatricolporites_sp

psilatricolporites_data = collect(dataset[dataset.species .== "Psilatricolporites_sp", Cols(x -> startswith(x, "d"))][1,:])

filter!(!ismissing, psilatricolporites_data)

psilatricolporites_θ1 = Normal(maximum(psilatricolporites_data)+200, 140)
psilatricolporites_θ2 = 0.0
psilatricolporites_λ = Normal(0, 2)

psilatricolporites_stratint = StratInterval(psilatricolporites_data,
                                            psilatricolporites_θ1,
                                            psilatricolporites_θ2,
                                            psilatricolporites_λ)

psilatricolporites_prior = sample_stratinterval(psilatricolporites_stratint,
                                                    100000,
                                                    NUTS(),
                                                    true, false)

psilatricolporites_posterior = sample_stratinterval(psilatricolporites_stratint,
                                                    100000,
                                                    NUTS(),
                                                    false, false)

CSV.write("../output/psilatricolporites_n71_posterior_stats.txt", DataFrame(describe(psilatricolporites_posterior)[1]))
CSV.write("../output/psilatricolporites_n71_posterior_quantiles.txt", DataFrame(describe(psilatricolporites_posterior)[2]))

plot(psilatricolporites_posterior)
savefig("../output/psilatricolporites_n71_posterior.svg")

histogram(psilatricolporites_data,
          normalize=:pdf,
          label=L"\textit{Psilatricolporites}"*" data",
          color=:lightgray,
          bins=range(minimum(psilatricolporites_data),
                     stop = maximum(psilatricolporites_data),
                     length = 15),
          xlabel="Depth (m)",
          ylabel="Density",
          ylim=[0, 0.04],
          xlim=[0, 1000])
density!(psilatricolporites_prior[:,:θ1,:], label="Prior θ1", color=:blue)
density!(psilatricolporites_posterior[:,:θ1,:], label="Posterior θ1", color=:black)
vline!([collect(hpd(psilatricolporites_posterior)[:θ1,:lower])[1],
        median(psilatricolporites_posterior[:,:θ1,:]),
        collect(hpd(psilatricolporites_posterior)[:θ1,:upper])[1]],
       line=:dash,
       color=:darkgray,
       label="Posterior quantiles\n0.25, 0.50, 0.75")
savefig("../output/psilatricolporites_n71_distibutions.svg")

density(psilatricolporites_prior[:,:λ,:], label="Prior λ", color=:blue)
density!(psilatricolporites_posterior[:,:λ,:], label="Posterior λ", color=:black)
savefig("../output/psilatricolporites_n71_lambda.svg")

# Arecipites_regio

arecipites_data = collect(dataset[dataset.species .== "Arecipites_regio", Cols(x -> startswith(x, "d"))][1,:])

filter!(!ismissing, arecipites_data)

arecipites_θ1 = Normal(maximum(arecipites_data)+200, 140)
arecipites_θ2 = Exponential(minimum(arecipites_data)/2)
arecipites_λ = Normal(0, 2)

arecipites_stratint = StratInterval(arecipites_data,
                                            arecipites_θ1,
                                            arecipites_θ2,
                                            arecipites_λ)

arecipites_prior = sample_stratinterval(arecipites_stratint,
                                                    100000,
                                                    NUTS(),
                                                    true, false)

arecipites_posterior = sample_stratinterval(arecipites_stratint,
                                                    100000,
                                                    NUTS(),
                                                    false, false)


CSV.write("../output/arecipites_n46_posterior_stats.txt", DataFrame(describe(arecipites_posterior)[1]))
CSV.write("../output/arecipites_n46_posterior_quantiles.txt", DataFrame(describe(arecipites_posterior)[2]))

plot(arecipites_posterior)
savefig("../output/arecipites_n46_sposterior.svg")

histogram(arecipites_data,
          normalize=:pdf,
          label=L"\textit{Arecipites}"*" "*L"\textit{regio}"*" data",
          color=:lightgray,
          bins=range(minimum(arecipites_data),
                     stop = maximum(arecipites_data),
                     length = 15),
          xlabel="Depth (m)",
          ylabel="Density",
          ylim=[0, 0.04],
          xlim=[0, 1000])
# plot the densities of theta1
density!(arecipites_prior[:,:θ1,:], label="Prior θ1 (right),\nθ2 (left)", color=:blue)
density!(arecipites_posterior[:,:θ1,:], label="Posterior θ1 (right),\nθ2 (left)", color=:black)
vline!([collect(hpd(arecipites_posterior)[:θ1,:lower])[1],
        median(arecipites_posterior[:,:θ1,:]),
        collect(hpd(arecipites_posterior)[:θ1,:upper])[1]],
       line=:dash,
       color=:darkgray,
       label="Posterior quantiles\n0.25, 0.50, 0.75")
# plot the densities of theta2
density!(arecipites_prior[:,:θ2,:], label="", color=:blue)
density!(arecipites_posterior[:,:θ2,:], label="", color=:black)
vline!([collect(hpd(arecipites_posterior)[:θ2,:lower])[1],
        median(arecipites_posterior[:,:θ2,:]),
        collect(hpd(arecipites_posterior)[:θ2,:upper])[1]],
       line=:dash,
       color=:darkgray,
       label="")
savefig("../output/arecipites_n46_distibutions.svg")

density(arecipites_prior[:,:λ,:], label="Prior λ", color=:blue)
density!(arecipites_posterior[:,:λ,:], label="Posterior λ", color=:black)
savefig("../output/arecipites_n46_lambda.svg")

"""
A medium-sized example of each fixing class
t2: Ischyosporites_problematicus
none: Scabratriporites_triangularis
"""

# Ischyosporites_problematicus

ischyosporites_data = collect(dataset[dataset.species .== "Ischyosporites_problematicus", Cols(x -> startswith(x, "d"))][1,:])

filter!(!ismissing, ischyosporites_data)

ischyosporites_θ1 = Normal(maximum(ischyosporites_data)+200, 140)
ischyosporites_θ2 = 0.0
ischyosporites_λ = Normal(0, 2)

ischyosporites_stratint = StratInterval(ischyosporites_data,
                                            ischyosporites_θ1,
                                            ischyosporites_θ2,
                                            ischyosporites_λ)

ischyosporites_prior = sample_stratinterval(ischyosporites_stratint,
                                                    100000,
                                                    NUTS(),
                                                    true, false)

ischyosporites_posterior = sample_stratinterval(ischyosporites_stratint,
                                                    100000,
                                                    NUTS(),
                                                    false, false)

CSV.write("../output/ischyosporites_n33_posterior_stats.txt", DataFrame(describe(ischyosporites_posterior)[1]))
CSV.write("../output/ischyosporites_n33_posterior_quantiles.txt", DataFrame(describe(ischyosporites_posterior)[2]))

plot(ischyosporites_posterior)
savefig("../output/ischyosporites_n33_posterior.svg")

histogram(ischyosporites_data,
          normalize=:pdf,
          label=L"\textit{Ischyosporites}"*" "*L"\textit{problematicus}"*" data",
          color=:lightgray,
          bins=range(minimum(ischyosporites_data),
                     stop = maximum(ischyosporites_data),
                     length = 15),
          xlabel="Depth (m)",
          ylabel="Density",
          ylim=[0, 0.04],
          xlim=[0, 1000])
density!(ischyosporites_prior[:,:θ1,:], label="Prior θ1", color=:blue)
density!(ischyosporites_posterior[:,:θ1,:], label="Posterior θ1", color=:black)
vline!([collect(hpd(ischyosporites_posterior)[:θ1,:lower])[1],
        median(ischyosporites_posterior[:,:θ1,:]),
        collect(hpd(ischyosporites_posterior)[:θ1,:upper])[1]],
       line=:dash,
       color=:darkgray,
       label="Posterior quantiles\n0.25, 0.50, 0.75")
savefig("../output/ischyosporites_n33_distibutions.svg")

density(ischyosporites_prior[:,:λ,:], label="Prior λ", color=:blue)
density!(ischyosporites_posterior[:,:λ,:], label="Posterior λ", color=:black)
savefig("../output/ischyosporites_n33_lambda.svg")

# Scabratriporites_triangularis

scabratriporites_data = collect(dataset[dataset.species .== "Scabratriporites_triangularis", Cols(x -> startswith(x, "d"))][1,:])

filter!(!ismissing, scabratriporites_data)

scabratriporites_θ1 = Normal(maximum(scabratriporites_data)+200, 140)
scabratriporites_θ2 = Exponential(minimum(scabratriporites_data)/2)
scabratriporites_λ = Normal(0, 2)

scabratriporites_stratint = StratInterval(scabratriporites_data,
                                            scabratriporites_θ1,
                                            scabratriporites_θ2,
                                            scabratriporites_λ)

scabratriporites_prior = sample_stratinterval(scabratriporites_stratint,
                                                    100000,
                                                    NUTS(),
                                                    true, false)

scabratriporites_posterior = sample_stratinterval(scabratriporites_stratint,
                                                    100000,
                                                    NUTS(),
                                                    false, false)

CSV.write("../output/scabratriporites_n26_posterior_stats.txt", DataFrame(describe(scabratriporites_posterior)[1]))
CSV.write("../output/scabratriporites_n26_posterior_quantiles.txt", DataFrame(describe(scabratriporites_posterior)[2]))

plot(scabratriporites_posterior)
savefig("../output/scabratriporites_n26_posterior.svg")

histogram(scabratriporites_data,
          normalize=:pdf,
          label=L"\textit{Scabratriporites}"*" "*L"\textit{triangularis}"*" data",
          color=:lightgray,
          bins=range(minimum(scabratriporites_data),
                     stop = maximum(scabratriporites_data),
                     length = 15),
          xlabel="Depth (m)",
          ylabel="Density",
          ylim=[0, 0.04],
          xlim=[0, 1000])
# plot the densities of theta1
density!(scabratriporites_prior[:,:θ1,:], label="Prior θ1 (right),\nθ2 (left)", color=:blue)
density!(scabratriporites_posterior[:,:θ1,:], label="Posterior θ1 (right),\nθ2 (left)", color=:black)
vline!([collect(hpd(scabratriporites_posterior)[:θ1,:lower])[1],
        median(scabratriporites_posterior[:,:θ1,:]),
        collect(hpd(scabratriporites_posterior)[:θ1,:upper])[1]],
       line=:dash,
       color=:darkgray,
       label="Posterior quantiles\n0.25, 0.50, 0.75")
# plot the densities of theta2
density!(scabratriporites_prior[:,:θ2,:], label="", color=:blue)
density!(scabratriporites_posterior[:,:θ2,:], label="", color=:black)
vline!([collect(hpd(scabratriporites_posterior)[:θ2,:lower])[1],
        median(scabratriporites_posterior[:,:θ2,:]),
        collect(hpd(scabratriporites_posterior)[:θ2,:upper])[1]],
       line=:dash,
       color=:darkgray,
       label="")
savefig("../output/scabratriporites_n26_distibutions.svg")

density(scabratriporites_prior[:,:λ,:], label="Prior λ", color=:blue)
density!(scabratriporites_posterior[:,:λ,:], label="Posterior λ", color=:black)
savefig("../output/scabratriporites_n26_lambda.svg")

"""
The lost sample of coal seam 115
We have a coal sample labeled as 115, but without stratigraphic information, what is the
best estimate of stratigraphic position we can have?

palynomorphs 1,2,5,6,12,13 are found in the sample

  1 │ Psilatricolporites_sp                 71  t2
  2 │ Retitricolporites_sp                  62  t2
  5 │ Psilamonocolpites_grandis             45  none
  6 │ Echinatisporis_minutus                42  none
 12 │ Retitricolpites_sp                    24  none
 13 │ Matonisporites_sp                     22  t2

conflate the posterior predictive stratigraphic intervals of these six taxa

but first remove the position corresponding to this layer, 258.88, column d72

everything needs to be recalculated, even if the first species had already been sampled, because we need to remove one occurrence
"""

# Psilatricolporites_sp t2
deleteat!(psilatricolporites_data, psilatricolporites_data .== 258.88)

psilatricolporites_θ1 = Normal(maximum(psilatricolporites_data)+200, 140)
psilatricolporites_θ2 = 0.0
psilatricolporites_λ = Normal(0, 2)

psilatricolporites_stratint = StratInterval(psilatricolporites_data,
                                            psilatricolporites_θ1,
                                            psilatricolporites_θ2,
                                            psilatricolporites_λ)

# Retitricolporites_sp t2
retitricolporites_data = collect(dataset[dataset.species .== "Retitricolporites_sp", Cols(x -> startswith(x, "d"))][1,:])

filter!(!ismissing, retitricolporites_data)

deleteat!(retitricolporites_data, retitricolporites_data .== 258.88)

retitricolporites_θ1 = Normal(maximum(retitricolporites_data)+200, 140)
retitricolporites_θ2 = 0.0
retitricolporites_λ = Normal(0, 2)

retitricolporites_stratint = StratInterval(retitricolporites_data,
                                            retitricolporites_θ1,
                                            retitricolporites_θ2,
                                            retitricolporites_λ)

# Psilamonocolpites_grandis none
psilamonocolpites_data = collect(dataset[dataset.species .== "Psilamonocolpites_grandis", Cols(x -> startswith(x, "d"))][1,:])

filter!(!ismissing, psilamonocolpites_data)

deleteat!(psilamonocolpites_data, psilamonocolpites_data .== 258.88)

psilamonocolpites_θ1 = Normal(maximum(psilamonocolpites_data)+200, 140)
psilamonocolpites_θ2 = Exponential(minimum(psilamonocolpites_data)/2)
psilamonocolpites_λ = Normal(0, 2)

psilamonocolpites_stratint = StratInterval(psilamonocolpites_data,
                                            psilamonocolpites_θ1,
                                            psilamonocolpites_θ2,
                                            psilamonocolpites_λ)

# Echinatisporis_minutus none
echinatisporis_data = collect(dataset[dataset.species .== "Echinatisporis_minutus", Cols(x -> startswith(x, "d"))][1,:])

filter!(!ismissing, echinatisporis_data)

deleteat!(echinatisporis_data, echinatisporis_data .== 258.88)

echinatisporis_θ1 = Normal(maximum(echinatisporis_data)+200, 140)
echinatisporis_θ2 = Exponential(minimum(echinatisporis_data)/2)
echinatisporis_λ = Normal(0, 2)

echinatisporis_stratint = StratInterval(echinatisporis_data,
                                            echinatisporis_θ1,
                                            echinatisporis_θ2,
                                            echinatisporis_λ)

# Retitricolpites_sp none
retitricolpites_data = collect(dataset[dataset.species .== "Retitricolpites_sp", Cols(x -> startswith(x, "d"))][1,:])

filter!(!ismissing, retitricolpites_data)

deleteat!(retitricolpites_data, retitricolpites_data .== 258.88)

retitricolpites_θ1 = Normal(maximum(retitricolpites_data)+200, 140)
retitricolpites_θ2 = Exponential(minimum(retitricolpites_data)/2)
retitricolpites_λ = Normal(0, 2)

retitricolpites_stratint = StratInterval(retitricolpites_data,
                                            retitricolpites_θ1,
                                            retitricolpites_θ2,
                                            retitricolpites_λ)

# Matonisporites_sp t2
matonisporites_data = collect(dataset[dataset.species .== "Matonisporites_sp", Cols(x -> startswith(x, "d"))][1,:])

filter!(!ismissing, matonisporites_data)

deleteat!(matonisporites_data, matonisporites_data .== 258.88)

matonisporites_θ1 = Normal(maximum(matonisporites_data)+200, 140)
matonisporites_θ2 = 0.0
matonisporites_λ = Normal(0, 2)

matonisporites_stratint = StratInterval(matonisporites_data,
                                            matonisporites_θ1,
                                            matonisporites_θ2,
                                            matonisporites_λ)

# combine all stratints in a vector
vecinterval = [psilatricolporites_stratint,
               retitricolporites_stratint,
               psilamonocolpites_stratint,
               echinatisporis_stratint,
               retitricolpites_stratint,
               matonisporites_stratint]

# sample the intervals and calculate the posterior predictives
# intervals 3, 4, and 5 had poor convergence, run for longer! 100000
coalseam115_postpredict_vec = sample_stratinterval(vecinterval, 100000, NUTS(), false, true)

# from tau_collection
postpredict_matrix = reduce(hcat, getindex.(coalseam115_postpredict_vec, 2))
# from conflate(Array{Float64,2})
interpolations = map(i -> InterpKDE(kde(postpredict_matrix[:, i])), 1:size(postpredict_matrix,2))
# numerat not needed beforehand
#evaluation of the product pdf is pdf(numerator, [vector of the same value, one per distribution])
denominat, err = quadgk(x -> prod(map(dd -> pdf(dd, x), interpolations)), -Inf, Inf, rtol=1e-8)

function wrapper_conflate(x)    
    prod(map(dd -> pdf(dd, x), interpolations))/denominat
end

# these two were found manually, but are what numerical optimisation should be doing
quadgk(x -> wrapper_conflate(x), -Inf, 54.5)[1]
quadgk(x -> wrapper_conflate(x), -Inf, 193.784)[1]
quadgk(x -> wrapper_conflate(x), 341.8, Inf)[1]
qlower = 54.4
qmedian = 193.784
qhigher = 341.8


# prepare to plot
xx = 0:1.0:800

yy = map(x -> tau_collection(coalseam115_postpredict_vec, x), xx)

density(coalseam115_postpredict_vec[1][2], label=L"\textit{Psilatricolporites}")
density!(coalseam115_postpredict_vec[2][2], label=L"\textit{Retitricolporites}")
density!(coalseam115_postpredict_vec[3][2], label=L"\textit{Psilamonocolpites}"*" "*L"\textit{grandis}")
density!(coalseam115_postpredict_vec[4][2], label=L"\textit{Echinatisporis}"*" "*L"\textit{minutus}")
density!(coalseam115_postpredict_vec[5][2], label=L"\textit{Retitricolpites}")
density!(coalseam115_postpredict_vec[6][2], label=L"\textit{Matonisporites}")
plot!(xx, yy, label="Coal seam 115 PDF", color=:black, width=2)
vline!([258.88], line=:dash, color=:black, label="Coal seal 115 true position")
vline!([qlower, qmedian, qhigher], line=:dash, color=:darkgrey, label="HPD interval and median")

savefig("../output/tau_coalseam115.svg");
