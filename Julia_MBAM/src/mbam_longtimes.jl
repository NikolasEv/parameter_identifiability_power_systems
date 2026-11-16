# MBAM of the infinite-bus synchronous generator, long times (t in [3, 5]).
# Extracted from the org-babel notebook used for the paper. Run from the src/ directory.
# Each section removes one parameter; the comment at the end of each section records
# the limit chosen from the geodesic.

######################################################################
# Setup Environment
######################################################################
# --- Package Activate
import Pkg; Pkg.activate("..")
# --- Imports
using LinearAlgebra
using Plots
using Models
using ParametricModels
using DataFrames
import DelimitedFiles: writedlm
using OrdinaryDiffEq
using Printf: @sprintf
import Geometry
import Geometry.Geodesics
# --- my_include
myinclude(path...) = include(joinpath(pwd(), path...))
# --- xform
xform(x) = x > 1 ? 1 / ( 1 + 1/x) : x/(1 + x)
mkpath("MBAM")  # diagnostic plots are saved here

######################################################################
# SG11
######################################################################
# --- SG11_modelfunctions
myinclude("pvalues.jl")
name = "SG11"
myinclude(@sprintf("%s.jl", name))
myinclude("ic.jl")
myinclude("rhs.jl")
myinclude("obs.jl")
t = range(3, 5, step = 0.02)
# --- SG11_define_f
function f(ps::SG, t)
    y = solve_dae(ps, ic, res, obs, t, differential_vars; abstol = 1e-12, reltol = 1e-12)
    return y
end

# --- SG11_save_trajectory
times = collect(0:0.01:6)
y = f(SG(), times)
writedlm(joinpath("..", "data", "SG11_trajectory.txt"), hcat(times, y))
# --- SG11_longtimes
pmodel = PModel(SG,
		  parameter_transforms,
		  f,
		  OLSData("SG", ModelArgs(t), f(SG{Float64}(), t)))
model = Models.Model(pmodel, "SG")
# --- SG_11xvalues
x = xvalues(pmodel)
pnames = ParametricModels.pnames(pmodel)
# --- SG11_Jacobians
j = model.jacobian(x)
writedlm(joinpath("..", "data", "SG11_jacobian.txt"), j)
g = j'j;
U,S,Vs = svd(j)
# --- SG11_Plot_MBAM_velocity
plot( Geometry.Geodesics.vi(x, j, model.Avv; i = 0), label="Singular Values")
xticks!(1:length(pnames), String.(pnames))
savefig(joinpath("MBAM", "SG11_velocity.png"))
# --- SG11_calc_geodesic
geo = Geometry.geodesic(model, x, Geometry.Geodesics.vi(x, j, model.Avv, i = 0), minratio = 0.01)
writedlm(joinpath("..", "data", "geo11.txt"), hcat(geo.ts, geo.τs, geo.xs, geo.vs))
# --- SG11_plot_ts_taus
plot(geo.ts, geo.τs)
savefig(joinpath("MBAM", "SG11_ts_taus.png"))

# --- SG11_plot_geodesic
plot(geo.τs, geo.xs)
savefig(joinpath("MBAM", "SG11_path.png"))
# Based on this, we are going to take y_8 to 0 (D to zero)

######################################################################
# SG10b
######################################################################
# --- SG10b_modelfunctions
myinclude("pvalues.jl")
name = "SG10b"
myinclude(@sprintf("%s.jl", name))
myinclude("ic.jl")
myinclude("rhs.jl")
myinclude("obs.jl")
t = range(3, 5, step = 0.02)
# --- SG10b_define_f
function f(ps::SG, t)
    y = solve_dae(ps, ic, res, obs, t, differential_vars; abstol = 1e-12, reltol = 1e-12)
    return y
end
# --- SG10b_longtimes
pmodel = PModel(SG,
		  parameter_transforms,
		  f,
		  OLSData("SG", ModelArgs(t), f(SG{Float64}(), t)))
model = Models.Model(pmodel, "SG")
# --- SG10b_xvalues
x = xvalues(pmodel)
pnames = ParametricModels.pnames(pmodel)
# --- SG10b_Jacobians
j = model.jacobian(x)
g = j'j;
U,S,Vs = svd(j)
# --- SG10b_Plot_MBAM_velocity
plot( Geometry.Geodesics.vi(x, j, model.Avv; i = 0, forward=true), label="Singular Values")
xticks!(1:length(pnames), String.(pnames))
savefig(joinpath("MBAM", "SG10b_velocity.png"))
# --- SG10b_calc_geodesic
geo = Geometry.geodesic(model, x, Geometry.Geodesics.vi(x, j, model.Avv, i = 0, forward=true), minratio = 1e-2)
writedlm(joinpath("..", "data", "geo10b.txt"), hcat(geo.ts, geo.τs, geo.xs, geo.vs))
# --- SG10b_plot_ts_taus
plot(geo.ts, geo.τs)
savefig(joinpath("MBAM", "SG10b_ts_taus.png"))

# --- SG10b_plot_geodesic
plot(geo.ts, geo.xs)
savefig(joinpath("MBAM", "SG10b_path.png"))
# Based on this, we are going to take y_7 to 0 (H to zero)

######################################################################
# SG09
######################################################################
# --- SG09_modelfunctions
myinclude("pvalues.jl")
name = "SG09"
myinclude(@sprintf("%s.jl", name))
myinclude("ic.jl")
myinclude("rhs.jl")
myinclude("obs.jl")
t = range(3, 5, step = 0.02)
# --- SG09_define_f
function f(ps::SG, t)
    y = solve_dae(ps, ic, res, obs, t, differential_vars; abstol = 1e-12, reltol = 1e-12)
    return y
end
# --- SG09_longtimes
pmodel = PModel(SG,
		  parameter_transforms,
		  f,
		  OLSData("SG", ModelArgs(t), f(SG{Float64}(), t)))
model = Models.Model(pmodel, "SG")
# --- SG09_xvalues
x = xvalues(pmodel)
pnames = ParametricModels.pnames(pmodel)
# --- SG09_Jacobians
j = model.jacobian(x)
g = j'j;
U,S,Vs = svd(j)
# --- SG09_Plot_MBAM_velocity
plot( Geometry.Geodesics.vi(x, j, model.Avv; i = 0), label="Singular Values")
xticks!(1:length(pnames), String.(pnames))
savefig(joinpath("MBAM", "SG09_velocity.png"))
# --- SG09_calc_geodesic
geo = Geometry.geodesic(model, x, Geometry.Geodesics.vi(x, j, model.Avv, i = 0), minratio = 1e-2)
writedlm(joinpath("..", "data", "geo09.txt"), hcat(geo.ts, geo.τs, geo.xs, geo.vs))
# --- SG09_plot_ts_taus
plot(geo.ts, geo.τs)
savefig(joinpath("MBAM", "SG09_ts_taus.png"))

# --- SG09_plot_geodesic
plot(geo.ts, geo.xs)
savefig(joinpath("MBAM", "SG09_path.png"))
# --- SG09_Plot_MBAM_velocityf
plot( geo.vs[end,:] )
xticks!(1:length(pnames), String.(pnames))
savefig(joinpath("MBAM", "SG09_velocityf.png"))
# Based on this, we are going to take y_7 to 0 (T2d0 to zero)

######################################################################
# SG08
######################################################################
# --- SG08_modelfunctions
myinclude("pvalues.jl")
name = "SG08"
myinclude(@sprintf("%s.jl", name))
myinclude("ic08.jl")
myinclude("rhs.jl")
myinclude("obs.jl")
t = range(3, 5, step = 0.02)
# --- SG08_define_f
function f(ps::SG, t)
    y = solve_dae(ps, ic, res, obs, t, differential_vars; abstol = 1e-12, reltol = 1e-12)
    return y
end
# --- SG08_longtimes
pmodel = PModel(SG,
		  parameter_transforms,
		  f,
		  OLSData("SG", ModelArgs(t), f(SG{Float64}(), t)))
model = Models.Model(pmodel, "SG")
# --- SG08_xvalues
x = xvalues(pmodel)
pnames = ParametricModels.pnames(pmodel)
# --- SG08_Jacobians
j = model.jacobian(x)
g = j'j;
U,S,Vs = svd(j)
# --- SG08_Plot_MBAM_velocity
plot( Geometry.Geodesics.vi(x, j, model.Avv; i = 0), label="Singular Values")
xticks!(1:length(pnames), String.(pnames))
savefig(joinpath("MBAM", "SG08_velocity.png"))
# --- SG08_calc_geodesic
geo = Geometry.geodesic(model, x, Geometry.Geodesics.vi(x, j, model.Avv, i = 0))
# --- SG08_plot_ts_taus
plot(geo.ts, geo.τs)
savefig(joinpath("MBAM", "SG08_ts_taus.png"))

# --- SG08_plot_geodesic
plot(geo.ts, geo.xs)
savefig(joinpath("MBAM", "SG08_path.png"))
# --- SG08_Plot_MBAM_velocityf
plot( geo.vs[end,:] )
xticks!(1:length(pnames), String.(pnames))
savefig(joinpath("MBAM", "SG08_velocityf.png"))
# Based on this, we are going to take y_7 to 0 (T2q0 to zero)

######################################################################
# SG07
######################################################################
# --- SG07_modelfunctions
myinclude("pvalues.jl")
name = "SG07"
myinclude(@sprintf("%s.jl", name))
myinclude("ic07.jl")
myinclude("rhs.jl")
myinclude("obs.jl")
t = range(3, 5, step = 0.02)
# --- SG07_define_f
function f(ps::SG, t)
    y = solve_dae(ps, ic, res, obs, t, differential_vars; abstol = 1e-12, reltol = 1e-12)
    return y
end
# --- SG07_longtimes
pmodel = PModel(SG,
		  parameter_transforms,
		  f,
		  OLSData("SG", ModelArgs(t), f(SG{Float64}(), t)))
model = Models.Model(pmodel, "SG")
# --- SG07_xvalues
x = xvalues(pmodel)
pnames = ParametricModels.pnames(pmodel)
# --- SG07_Jacobians
j = model.jacobian(x)
g = j'j;
U,S,Vs = svd(j)
# --- SG07_Plot_MBAM_velocity
plot( Geometry.Geodesics.vi(x, j, model.Avv; i = 0), label="Singular Values")
xticks!(1:length(pnames), String.(pnames))
savefig(joinpath("MBAM", "SG07_velocity.png"))
# --- SG07_calc_geodesic
geo = Geometry.geodesic(model, x, Geometry.Geodesics.vi(x, j, model.Avv, i = 0))
# --- SG07_plot_ts_taus
plot(geo.ts, geo.τs)
savefig(joinpath("MBAM", "SG07_ts_taus.png"))
# --- SG07_plot_geodesic
plot(geo.ts, geo.xs)
savefig(joinpath("MBAM", "SG07_path.png"))
# --- SG07_Plot_MBAM_velocityf
plot( geo.vs[end,:] )
xticks!(1:length(pnames), String.(pnames))
savefig(joinpath("MBAM", "SG07_velocityf.png"))
# Based on this, we are going to take y_1 to 0 (deltax1 to zero)

######################################################################
# SG06
######################################################################
# --- SG06_modelfunctions
myinclude("pvalues.jl")
name = "SG06"
myinclude(@sprintf("%s.jl", name))
myinclude("ic07.jl")
myinclude("rhs.jl")
myinclude("obs.jl")
t = range(3, 5, step = 0.02)
# --- SG06_define_f
function f(ps::SG, t)
    y = solve_dae(ps, ic, res, obs, t, differential_vars; abstol = 1e-12, reltol = 1e-12)
    return y
end
# --- SG06_longtimes
pmodel = PModel(SG,
		  parameter_transforms,
		  f,
		  OLSData("SG", ModelArgs(t), f(SG{Float64}(), t)))
model = Models.Model(pmodel, "SG")
# --- SG09_save_trajectory
times = collect(0:0.01:6)
y = f(SG(), times)
writedlm(joinpath("..", "data", "SG06_trajectory.txt"), hcat(times, y))
# --- SG06_xvalues
x = xvalues(pmodel)
pnames = ParametricModels.pnames(pmodel)
# --- SG06_Jacobians
j = model.jacobian(x)
writedlm(joinpath("..", "data", "SG06_jacobian.txt"), j)
g = j'j;
U,S,Vs = svd(j)
