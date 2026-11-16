import Pkg
Pkg.activate("Plots"; shared = true)
using LinearAlgebra
using LaTeXStrings
import DelimitedFiles: readdlm

using CairoMakie
CairoMakie.activate!()

# using GLMakie
# GLMakie.activate!()

j = readdlm(joinpath(@__DIR__, "..", "data", "SG06_jacobian.txt"))
u,s,v = svd(j)

# Figure dimensions in inches
ϕ = 0.5*(sqrt(5) - 1) # An appealing ratio
# Typographic units: 1 pt = 1/72 of an inch
scale = 4 # Visual scale for screen
dpi = 300 
w = scale * round(Int, (3 + 3/8)*72) # Singule column of revtex
h = round(Int, w * ϕ) # An appealing ratio
dpi = 300

function plot_information_spectrum(s, v)
    f = Figure(size=(w, h))
    gl = GridLayout()
    ax1 = Axis(f[1,1],
               yscale = Makie.log10,
               yticks = (exp10.([-8, -5, -2, 1]),
                         [L"10^{-8}",
                          L"10^{-5}",
                          L"10^{-2}",
                          L"10^1",
                          ])
               )
    for si ∈  s      
        lines!(ax1, [0.15, 0.85], [si^2, si^2], color = "black", linewidth = 3)
    end    
    xlims!(ax1, 0..1)
    ylims!(ax1, 2.4207086839812423e-10, 700.9710575274788) # Same limits as for 11 parameters
    hidexdecorations!(ax1)
    ax2 = Axis(f[1,2][1,1],
               aspect = 1,
               xticks = 1:length(s),
               yticks = (reverse(1:length(s)), [L"\delta x_2", L"\delta x_3", L"\delta x_4", L"x''_d", L"\delta T_d", L"\delta T_q"]),
               xticklabelsvisible = false,
               xlabel = L"stiff $\leftarrow$ Eigenmodes $\rightarrow$ sloppy"
               )
    heatmap!(1:length(s), reverse(1:length(s)), v'.^2)
    Colorbar(f[1,2][1,2], limits = (0, 1), ticks = [0, 0.25, 0.5, 0.75, 1.0])
    colsize!(f.layout, 1, Fixed(round(Int, w*0.1)))    
    colsize!(f.layout, 2, Fixed(round(Int, w*0.7)))
    colgap!(f.layout, 0)
    f
end    
f = with_theme(()->plot_information_spectrum(s,v), merge(theme_latexfonts(), Theme(fontsize = scale*8)))
save("SG06_FIM.png", f, px_per_unit = dpi*(3 + 3/8)/w)
