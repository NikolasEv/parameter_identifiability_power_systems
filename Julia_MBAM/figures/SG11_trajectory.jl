import Pkg
Pkg.activate("Plots"; shared = true)
using LinearAlgebra
using LaTeXStrings
import DelimitedFiles: readdlm

using CairoMakie
CairoMakie.activate!()

# using GLMakie
# GLMakie.activate!()

data = readdlm(joinpath(@__DIR__, "..", "data", "SG11_trajectory.txt"))

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
    ax = Axis(f[1,1], xlabel = L"t\ (sec)")
    delta = lines!(ax, data[:,1], data[:,2], label = L"\delta", linewidth = 4)
    omega = lines!(ax, data[:,1], data[:,3], label = L"\omega", linewidth = 4)
    e1q = lines!(ax, data[:,1], data[:,4], label = L"e'_q", linewidth = 4)
    e1d = lines!(ax, data[:,1], data[:,5], label = L"e'_d", linewidth = 4)
    e2q = lines!(ax, data[:,1], data[:,6], label = L"e''_q", linewidth = 4)
    e2d = lines!(ax, data[:,1], data[:,7], label = L"e''_d", linewidth = 4)
    xlims!(ax, low = 0, high = 6)
    axislegend(ax)
    f
end

save("SG11_trajectory.png", f, px_per_unit = dpi*(3 + 3/8)/w)
