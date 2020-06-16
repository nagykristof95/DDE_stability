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
includet("ISIM_DE_append.jl")
includet("SDM.jl")
includet("SDM_append.jl")
using Main.SIDE
using Main.SILMRK
using Main.SI
using Main.SIDEA
using Main.SD
using Main.SDapp

rng = MersenneTwister(1234)
myred=RGB(214/255,34/255,36/255)
myblue=RGB(14/255,97/255,134/255)
####### Multi-Dimensional Bisection Method #########

function foo(x,y)
    return(norm(BAS.normmax(SDM_full([omega kappa x epsilon y tau0],(n,gmax,mult,"Full"))[:,end]))-1)
end

function foo(x,y)
    return(norm(BAS.normmax(SDM_full_app([omega kappa x epsilon y tau0],(n,gmax,mult))[:,end]))-1)
end

function foo(x,y)
    return(norm(normmax(ISIMcont([omega kappa x epsilon y tau0])[:,end]))-1)
end

function foo(x,y)
    return(norm(normmax(ISIM_LSM([omega kappa x epsilon y tau0])[:,end]))-1)
end

function foo(x,y)
    return(norm(BAS.normmax(ISIM_DE([omega kappa x epsilon y tau0],(400,12,12,ALG0),(1e-8,1e-8))[:,end]))-1)
end

ax1=Axis(0.0:0.5:20.0,"delta") # initial grid in x direction
ax2=Axis(-1.0:0.25:1.0,"b") # initial grid in y direction

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

T_FULL=120.921822
x_FULL,y_FULL

T_1=91.822125
x_1,y_1

T_2=94.311763
x_2,y_2

T_3=94.691673
x_3,y_3

fig = figure(1);clf()
Plots.scatter(x_FULL,y_FULL,markersize=:3,markercolor=myred,markerstrokewidth = 0)
Plots.scatter!(x_1,y_1,markersize=:3,markercolor=myblue,markerstrokewidth = 0)

PyPlot.savefig("1.pdf")

ax1=Axis(1000.0:100.0:4000.0,"delta") # initial grid in x direction
ax2=Axis(0.0:0.0005:0.005,"b") # initial grid in y direction

function foo(x,y)
    return(norm(BAS.normmax(ISIM_DE([mx my kx ky sx sy kt kr d ae z strat x y],(100,12,12,ALG0),(1e-8,1e-8))[:,end]))-1)
end

function foo(x,y)
    return(norm(BAS.normmax(ISIM_DE([mx my kx ky sx sy kt kr d ae z strat x y],(100,12,12,ALG0),(1e-4,1e-4))[:,end]))-1)
end

mymdbm2=MDBM_Problem(foo,[ax1,ax2])
iteration=3 #number of refinements (resolution doubling)
@time MDBM.solve!(mymdbm2,iteration)

x_eval2,y_eval2=getevaluatedpoints(mymdbm2)
x_sol2,y_sol2=getinterpolatedsolution(mymdbm2)

fig = figure(1);clf()
PyPlot.scatter(x_eval2,y_eval2,s=5)
PyPlot.scatter(x_sol2,y_sol2,s=5)
