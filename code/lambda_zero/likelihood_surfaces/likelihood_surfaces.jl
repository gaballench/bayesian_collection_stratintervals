using StratIntervals
using Distributions
using Random
using StatsPlots
using Optim
using StatsBase
using Turing 
using ColorSchemes

Random.seed!(1985)

# true params
theta2 = 100.0
theta1 = 150.0
lambda = 0.0
ndata = [2, 10, 50, 100, 1000, 10000]

data2 = rand(ThreeParBeta(theta1, theta2, lambda), ndata[1])
data10 = rand(ThreeParBeta(theta1, theta2, lambda), ndata[2])
data50 = rand(ThreeParBeta(theta1, theta2, lambda), ndata[3])
data100 = rand(ThreeParBeta(theta1, theta2, lambda), ndata[4])
data1000 = rand(ThreeParBeta(theta1, theta2, lambda), ndata[5])
data10000 = rand(ThreeParBeta(theta1, theta2, lambda), ndata[6])
#histogram(data10000)

### contour plots for pairs of params

function theta12_surface(theta1, theta2, lambda, data)
    sum(logpdf.(ThreeParBeta(theta1, theta2, lambda), data))
end

### lnL surfaces as a function of sample size

### All of these seem to be for the supplementary file, except for maybe the n=1000 example

"""
### data2
# theta1 and lambda
xs = range(maximum(data2), maximum(data2)+1.2, length=30)
ys = range(-2,2, length=30)
zs = [theta12_surface(x, theta2, y, data2) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\lambda\$", plot_title="n = 2", cbar=false, clabels=true)
scatter!([theta1], [lambda], color="black", legend=false)
savefig("n2_theta1_lambda_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\lambda\$", zlabel="lnL", plot_title="n = 2", cbar=false)
scatter!([theta1], [lambda], [theta12_surface(theta1, theta2, lambda, data2)], color="black", legend=false)
savefig("n2_theta1_lambda_surface.svg")

# theta2 and lambda
xs = range(minimum(data2)-40, minimum(data2), length=30)
ys = range(-2,2, length=30)
zs = [theta12_surface(theta1, x, y, data2) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_2\$", ylabel="\$\\lambda\$", plot_title="n = 2", cbar=false, clabels=true)
scatter!([theta2], [lambda], color="black", legend=false)
savefig("n2_theta2_lambda_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_2\$", ylabel="\$\\lambda\$", zlabel="lnL", plot_title="n = 2", cbar=false)
scatter!([theta2], [lambda], [theta12_surface(theta1, theta2, lambda, data2)], color="black", legend=false)
savefig("n2_theta2_lambda_surface.svg")

# theta1 and theta2
xs = range(maximum(data2), maximum(data2)+1.2, length=30)
ys = range(minimum(data2)-40, minimum(data2), length=30)
zs = [theta12_surface(x, y, lambda, data2) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\theta_2\$", plot_title="n = 2", cbar=false, clabels=true)
scatter!([theta1], [theta2], color="black", legend=false)
savefig("n2_theta1_theta2_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\theta_2\$", zlabel="lnL", plot_title="n = 2", cbar=false)
scatter!([theta1], [theta2], [theta12_surface(theta1, theta2, lambda, data2)], color="white", legend=false)
savefig("n2_theta1_theta2_surface.svg")
"""

### data10
# theta1 and lambda
xs = range(maximum(data10), maximum(data10)+2, length=30)
ys = range(-2,2, length=30)
zs = [theta12_surface(x, theta2, y, data10) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\lambda\$", plot_title="n = 10", cbar=false, clabels=true)
scatter!([theta1], [lambda], color="black", legend=false)
savefig("n10_theta1_lambda_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\lambda\$", zlabel="lnL", plot_title="n = 10", cbar=false)
scatter!([theta1], [lambda], [theta12_surface(theta1, theta2, lambda, data10)], color="black", legend=false)
savefig("n10_theta1_lambda_surface.svg")

# theta2 and lambda
xs = range(minimum(data10)-12, minimum(data10), length=30)
ys = range(-2,2, length=30)
zs = [theta12_surface(theta1, x, y, data10) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_2\$", ylabel="\$\\lambda\$", plot_title="n = 10", cbar=false, clabels=true)
scatter!([theta2], [lambda], color="black", legend=false)
savefig("n10_theta2_lambda_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_2\$", ylabel="\$\\lambda\$", zlabel="lnL", plot_title="n = 10", cbar=false)
scatter!([theta2], [lambda], [theta12_surface(theta1, theta2, lambda, data10)], color="black", legend=false)
savefig("n10_theta2_lambda_surface.svg")

# theta1 and theta2
xs = range(maximum(data10), maximum(data10)+2, length=30)
ys = range(minimum(data10)-12, minimum(data10), length=30)
zs = [theta12_surface(x, y, lambda, data10) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\theta_2\$", plot_title="n = 10", cbar=false, clabels=true)
scatter!([theta1], [theta2], color="black", legend=false)
savefig("n10_theta1_theta2_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\theta_2\$", zlabel="lnL", plot_title="n = 10", cbar=false)
scatter!([theta1], [theta2], [theta12_surface(theta1, theta2, lambda, data10)], color="white", legend=false)
savefig("n10_theta1_theta2_surface.svg")

### data50
# theta1 and lambda
xs = range(maximum(data50), maximum(data50)+1.2, length=30)
ys = range(-2,2, length=30)
zs = [theta12_surface(x, theta2, y, data50) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\lambda\$", plot_title="n = 50", cbar=false, clabels=true)
scatter!([theta1], [lambda], color="black", legend=false)
savefig("n50_theta1_lambda_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\lambda\$", zlabel="lnL", plot_title="n = 50", cbar=false)
scatter!([theta1], [lambda], [theta12_surface(theta1, theta2, lambda, data50)], color="black", legend=false)
savefig("n50_theta1_lambda_surface.svg")

# theta2 and lambda
xs = range(minimum(data50)-1, minimum(data50), length=30)
ys = range(-2,2, length=30)
zs = [theta12_surface(theta1, x, y, data50) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_2\$", ylabel="\$\\lambda\$", plot_title="n = 50", cbar=false, clabels=true)
scatter!([theta2], [lambda], color="black", legend=false)
savefig("n50_theta2_lambda_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_2\$", ylabel="\$\\lambda\$", zlabel="lnL", plot_title="n = 50", cbar=false)
scatter!([theta2], [lambda], [theta12_surface(theta1, theta2, lambda, data50)], color="black", legend=false)
savefig("n50_theta2_lambda_surface.svg")

# theta1 and theta2
xs = range(maximum(data50), maximum(data50)+1.2, length=30)
ys = range(minimum(data50)-1, minimum(data50), length=30)
zs = [theta12_surface(x, y, lambda, data50) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\theta_2\$", plot_title="n = 50", cbar=false, clabels=true)
scatter!([theta1], [theta2], color="black", legend=false)
savefig("n50_theta1_theta2_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\theta_2\$", zlabel="lnL", plot_title="n = 50", cbar=false)
scatter!([theta1], [theta2], [theta12_surface(theta1, theta2, lambda, data50)], color="white", legend=false)
savefig("n50_theta1_theta2_surface.svg")

### data100
# theta1 and lambda
xs = range(maximum(data100), maximum(data100)+0.2, length=30)
ys = range(-2,2, length=30)
zs = [theta12_surface(x, theta2, y, data100) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\lambda\$", plot_title="n = 100", cbar=false, clabels=true)
scatter!([theta1], [lambda], color="black", legend=false)
savefig("n100_theta1_lambda_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\lambda\$", zlabel="lnL", plot_title="n = 100", cbar=false)
scatter!([theta1], [lambda], [theta12_surface(theta1, theta2, lambda, data100)], color="black", legend=false)
savefig("n100_theta1_lambda_surface.svg")

# theta2 and lambda
xs = range(minimum(data100)-0.2, minimum(data100), length=30)
ys = range(-2,2, length=30)
zs = [theta12_surface(theta1, x, y, data100) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_2\$", ylabel="\$\\lambda\$", plot_title="n = 100", cbar=false, clabels=true)
scatter!([theta2], [lambda], color="black", legend=false)
savefig("n100_theta2_lambda_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_2\$", ylabel="\$\\lambda\$", zlabel="lnL", plot_title="n = 100", cbar=false)
scatter!([theta2], [lambda], [theta12_surface(theta1, theta2, lambda, data100)], color="black", legend=false)
savefig("n100_theta2_lambda_surface.svg")

# theta1 and theta2
xs = range(maximum(data100), maximum(data100)+0.2, length=30)
ys = range(minimum(data100)-0.2, minimum(data100), length=30)
zs = [theta12_surface(x, y, lambda, data100) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\theta_2\$", plot_title="n = 100", cbar=false, clabels=true)
scatter!([theta1], [theta2], color="black", legend=false)
savefig("n100_theta1_theta2_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\theta_2\$", zlabel="lnL", plot_title="n = 100", cbar=false)
scatter!([theta1], [theta2], [theta12_surface(theta1, theta2, lambda, data100)], color="white", legend=false)
savefig("n100_theta1_theta2_surface.svg")

### data1000
# theta1 and lambda
xs = range(maximum(data1000), maximum(data1000)+0.2, length=30)
ys = range(-2,2, length=30)
zs = [theta12_surface(x, theta2, y, data1000) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\lambda\$", plot_title="n = 1000", cbar=false, clabels=true)
scatter!([theta1], [lambda], color="black", legend=false)
savefig("n1000_theta1_lambda_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\lambda\$", zlabel="lnL", plot_title="n = 1000", cbar=false)
scatter!([theta1], [lambda], [theta12_surface(theta1, theta2, lambda, data1000)], color="black", legend=false)
savefig("n1000_theta1_lambda_surface.svg")

# theta2 and lambda
xs = range(minimum(data1000)-0.2, minimum(data1000), length=30)
ys = range(-2,2, length=30)
zs = [theta12_surface(theta1, x, y, data1000) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_2\$", ylabel="\$\\lambda\$", plot_title="n = 1000", cbar=false, clabels=true)
scatter!([theta2], [lambda], color="black", legend=false)
savefig("n1000_theta2_lambda_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_2\$", ylabel="\$\\lambda\$", zlabel="lnL", plot_title="n = 1000", cbar=false)
scatter!([theta2], [lambda], [theta12_surface(theta1, theta2, lambda, data1000)], color="black", legend=false)
savefig("n1000_theta2_lambda_surface.svg")

# theta1 and theta2
xs = range(maximum(data1000), maximum(data1000)+0.2, length=30)
ys = range(minimum(data1000)-0.2, minimum(data1000), length=30)
zs = [theta12_surface(x, y, lambda, data1000) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\theta_2\$", plot_title="n = 1000", cbar=false, clabels=true)
scatter!([theta1], [theta2], color="black", legend=false)
savefig("n1000_theta1_theta2_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\theta_2\$", zlabel="lnL", plot_title="n = 1000", cbar=false)
scatter!([theta1], [theta2], [theta12_surface(theta1, theta2, lambda, data1000)], color="white", legend=false)
savefig("n1000_theta1_theta2_surface.svg")

### data10000
# theta1 and lambda
xs = range(maximum(data10000), maximum(data10000)+0.2, length=30)
ys = range(-2,2, length=30)
zs = [theta12_surface(x, theta2, y, data10000) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\lambda\$", plot_title="n = 10000", cbar=false, clabels=true)
scatter!([theta1], [lambda], color="black", legend=false)
savefig("n10000_theta1_lambda_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\lambda\$", zlabel="lnL", plot_title="n = 10000", cbar=false)
scatter!([theta1], [lambda], [theta12_surface(theta1, theta2, lambda, data10000)], color="black", legend=false)
savefig("n10000_theta1_lambda_surface.svg")

# theta2 and lambda
xs = range(minimum(data10000)-0.2, minimum(data10000), length=30)
ys = range(-2,2, length=30)
zs = [theta12_surface(theta1, x, y, data10000) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_2\$", ylabel="\$\\lambda\$", plot_title="n = 10000", cbar=false, clabels=true)
scatter!([theta2], [lambda], color="black", legend=false)
savefig("n10000_theta2_lambda_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_2\$", ylabel="\$\\lambda\$", zlabel="lnL", plot_title="n = 10000", cbar=false)
scatter!([theta2], [lambda], [theta12_surface(theta1, theta2, lambda, data10000)], color="black", legend=false)
savefig("n10000_theta2_lambda_surface.svg")

# theta1 and theta2
xs = range(maximum(data10000), maximum(data10000)+0.2, length=30)
ys = range(minimum(data10000)-0.2, minimum(data10000), length=30)
zs = [theta12_surface(x, y, lambda, data10000) for y in ys, x in xs]

contourf(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\theta_2\$", plot_title="n = 10000", cbar=false, clabels=true)
scatter!([theta1], [theta2], color="black", legend=false)
savefig("n10000_theta1_theta2_contour.svg")
surface(xs, ys, zs, xlabel="\$\\theta_1\$", ylabel="\$\\theta_2\$", zlabel="lnL", plot_title="n = 10000", cbar=false)
scatter!([theta1], [theta2], [theta12_surface(theta1, theta2, lambda, data10000)], color="white", legend=false)
savefig("n10000_theta1_theta2_surface.svg")

########## lnL surfaces for individual parameters as a function of sample size

# lambda as a function of sample size
plot(x -> sum(logpdf.(ThreeParBeta(theta1, theta2, x), data10)), -4, 4, xlabel="\$\\lambda\$", ylabel="lnL", palette=:rainbow, legend=:right, label="n=10", plot_title="Marginal \$\\lambda\$")
#plot!(x -> sum(logpdf.(ThreeParBeta(theta1, theta2, x), data10)), -4, 4, label="n=10")
plot!(x -> sum(logpdf.(ThreeParBeta(theta1, theta2, x), data50)), -4, 4, label="n=50")
plot!(x -> sum(logpdf.(ThreeParBeta(theta1, theta2, x), data100)), -4, 4, label="n=100")
plot!(x -> sum(logpdf.(ThreeParBeta(theta1, theta2, x), data1000)), -4, 4, label="n=1000")
#plot!(x -> sum(logpdf.(ThreeParBeta(theta1, theta2, x), data10000)), -4, 4)
vline!([lambda], color="black", linestyle=:dash, label="true \$\\lambda\$")
savefig("marginal_lambda_vs_n.svg")

# are thetae unidentifiable? wtf? theta1 isn't, theta2 is only closer to 0.0!

# theta1 as a function of sample size
plot(x -> sum(logpdf.(ThreeParBeta(x, theta2, lambda), data10)), maximum(data10), 200, xlabel="\$\\theta_1\$", ylabel="lnL", palette=:rainbow, legend=:right, label="n=10", plot_title="Marginal \$\\theta_1\$")
#plot!(x -> sum(logpdf.(ThreeParBeta(x, theta2, lambda), data10)), maximum(data10), 200, label="n=10")
plot!(x -> sum(logpdf.(ThreeParBeta(x, theta2, lambda), data50)), maximum(data50), 200, label="n=50")
plot!(x -> sum(logpdf.(ThreeParBeta(x, theta2, lambda), data100)), maximum(data100), 200, label="n=100")
plot!(x -> sum(logpdf.(ThreeParBeta(x, theta2, lambda), data1000)), maximum(data1000), 200, label="n=1000")
#plot(x -> sum(logpdf.(ThreeParBeta(x, theta2, lambda), data10000)), maximum(data10000), 200)
vline!([theta1], color="black", linestyle=:dash, label="true \$\\theta_1\$")
savefig("marginal_theta1_vs_n.svg")

# theta2 as a function of sample size
plot(x -> sum(logpdf.(ThreeParBeta(theta1, x, lambda), data10)), 50.0, minimum(data10), xlabel="\$\\theta_2\$", ylabel="lnL", palette=:rainbow, legend=:right, label="n=10", plot_title="Marginal \$\\theta_2\$")
#plot!(x -> sum(logpdf.(ThreeParBeta(theta1, x, lambda), data10)), 50.0, minimum(data10), label="n=10")
plot!(x -> sum(logpdf.(ThreeParBeta(theta1, x, lambda), data50)), 50.0, minimum(data50), label="n=50")
plot!(x -> sum(logpdf.(ThreeParBeta(theta1, x, lambda), data100)), 50.0, minimum(data100), label="n=100")
plot!(x -> sum(logpdf.(ThreeParBeta(theta1, x, lambda), data1000)), 50.0, minimum(data1000), label="n=1000")
#plot!(x -> sum(logpdf.(ThreeParBeta(theta1, x, lambda), data10000)), 50.0, minimum(data10000))
vline!([theta2], color="black", linestyle=:dash, label="true \$\\theta_2\$")
savefig("marginal_theta2_vs_n.svg")

###### subfigures for the composite plot of the lnL surface of the smallest data size 

plot(x -> sum(logpdf.(ThreeParBeta(theta1, theta2, x), data10)), -4, 4, xlabel="\$\\lambda\$", ylabel="lnL", plot_title="n=10")
vline!([lambda], color="black", linestyle=:dash, legend=false)
savefig("marginal_lambda_vs_small_n.svg")

plot(x -> sum(logpdf.(ThreeParBeta(x, theta2, lambda), data2)), maximum(data10), 151, xlabel="\$\\theta_1\$", ylabel="lnL", plot_title="n=10")
vline!([theta1], color="black", linestyle=:dash, legend=false)
savefig("marginal_theta1_vs_small_n.svg")

plot(x -> sum(logpdf.(ThreeParBeta(theta1, x, lambda), data2)), 99, minimum(data10), xlabel="\$\\theta_2\$", ylabel="lnL", plot_title="n=10")
vline!([theta2], color="black", linestyle=:dash, legend=false)
savefig("marginal_theta2_vs_small_n.svg")

###### subfigures for the composite plot of the lnL surface of the largest data size 

plot(x -> sum(logpdf.(ThreeParBeta(theta1, theta2, x), data1000)), -4, 4, xlabel="\$\\lambda\$", ylabel="lnL", plot_title="n=1000")
vline!([lambda], color="black", linestyle=:dash, legend=false)
savefig("marginal_lambda_vs_large_n.svg")

plot(x -> sum(logpdf.(ThreeParBeta(x, theta2, lambda), data1000)), maximum(data1000), 151, xlabel="\$\\theta_1\$", ylabel="lnL", plot_title="n=1000")
vline!([theta1], color="black", linestyle=:dash, legend=false)
savefig("marginal_theta1_vs_large_n.svg")

plot(x -> sum(logpdf.(ThreeParBeta(theta1, x, lambda), data1000)), 99, minimum(data1000), xlabel="\$\\theta_2\$", ylabel="lnL", plot_title="n=1000")
vline!([theta2], color="black", linestyle=:dash, legend=false)
savefig("marginal_theta2_vs_large_n.svg")
