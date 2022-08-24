import Pkg
Pkg.add("Plots")
Pkg.add("PyPlot")

using Plots, Random
path   = "C:/Users/Reza/Desktop/juliatesting/media/" # Absolute, change this

Random.seed!(2018)

x = randn(1000)
y = randn(1000)
z = randn(1000)

histogram(x, bins=20, alpha=0.4, label="A")
histogram!(y, bins=20, alpha=0.6, label="B")
hist_1 = histogram!(z, bins=0:8, alpha=1, label="C")

# Makie
using CairoMakie
CairoMakie.activate!(type = "png")

data = randn(1000)
his = hist(data, normalization = :pdf, bar_labels = :values,
     label_formatter=x-> round(x, digits=2), label_size = 15,
     strokewidth = 0.5, strokecolor = (:white, 0.5), color = :values)

save(path * "histogram.png"   , his   )
save(path * "histogram_2.png", hist_1)
