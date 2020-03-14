using Printf
using Statistics
using LinearAlgebra
using MDBM
using Interpolations
using Plots
pyplot()
using PyPlot
pygui(true);
using DelimitedFiles
using Random
using DifferentialEquations
using Revise
includet("ISIM_DE_sol.jl")
includet("ISIM_LMS_RK.jl")
includet("ISIM.jl")
using Main.SIDE
using Main.SILMRK
using Main.SI

rng = MersenneTwister(1234)

####### Multi-Dimensional Bisection Method #########

function foo(x,y)
    return(norm(BAS.normmax(ISIM([omega kappa x epsilon y tau0],(n,gmax,mult,BAS.BRK4))[:,end]))-1)
end


function foo(x,y)
    return(norm(normmax(ISIMcont([omega kappa x epsilon y tau0])[:,end]))-1)
end

function foo(x,y)
    return(norm(normmax(ISIM_LSM([omega kappa x epsilon y tau0])[:,end]))-1)
end

ax1=Axis(-1.0:0.5:5,"delta") # initial grid in x direction
ax2=Axis(-1.5:0.5:1.5,"b") # initial grid in y direction

ax1=Axis(0.1:0.25:5.0,"tau") # initial grid in x direction
ax2=Axis(0.05:0.2:1.5,"b") # initial grid in y direction

mymdbm=MDBM_Problem(foo,[ax1,ax2])
iteration=3 #number of refinements (resolution doubling)
@time MDBM.solve!(mymdbm,iteration)

x_eval,y_eval=getevaluatedpoints(mymdbm)
x_sol,y_sol=getinterpolatedsolution(mymdbm)

fig = figure(1);clf()
PyPlot.scatter(x_eval,y_eval,s=5)
PyPlot.scatter(x_sol,y_sol,s=5)
