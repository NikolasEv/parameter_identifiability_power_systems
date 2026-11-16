import Pkg
Pkg.activate("Plots"; shared = true)
using LinearAlgebra
using LaTeXStrings
import DelimitedFiles: readdlm

using CairoMakie
CairoMakie.activate!()

# using GLMakie
# GLMakie.activate!()

data = readdlm(joinpath(@__DIR__, "..", "data", "geo09.txt"))

# Figure dimensions in inches
ϕ = 0.5*(sqrt(5) - 1) # An appealing ratio
# Typographic units: 1 pt = 1/72 of an inch
scale = 4 # Visual scale for screen
dpi = 300 
w = scale * round(Int, (3 + 3/8)*72) # Singule column of revtex
h = round(Int, w * ϕ) # An appealing ratio
dpi = 300

f = with_theme(merge(theme_latexfonts(), Theme(fontsize = scale*8))) do
    f = Figure(size=(w, h))
    ax = Axis(f[1,1], xlabel = L"\tau")
    dx1 = lines!(ax, data[:,2], data[:,3], label = L"\delta x_1", linewidth = 4)
    dx2 = lines!(ax, data[:,2], data[:,4], label = L"\delta x_2", linewidth = 4)
    dx3 = lines!(ax, data[:,2], data[:,5], label = L"\delta x_3", linewidth = 4)
    dx4 = lines!(ax, data[:,2], data[:,6], label = L"\delta x_4", linewidth = 4)
    x2d = lines!(ax, data[:,2], data[:,7], label = L"x''_d", linewidth = 4)
    dTd = lines!(ax, data[:,2], data[:,8], label = L"\delta T_d", linewidth = 4)
    dTq = lines!(ax, data[:,2], data[:,9], label = L"\delta T_q", linewidth = 4)
    Td0 = lines!(ax, data[:,2], data[:,10], label = L"T''_{d0}", linewidth = 4, linestyle = :dash)
    Td0 = lines!(ax, data[:,2], data[:,11], label = L"T''_{q0}", linewidth = 4, linestyle = :dash)
    # xlims!(ax, low = 0, high = 6)
    # f[1,2] = Legend(f, ax, framevisible = false)
    axislegend(ax, position = :lb, nbanks = 2, patchsize = (40, 20))
    f
end

save("SG09_geodesic.png", f, px_per_unit = dpi*(3 + 3/8)/w)
