using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using Random
using DifferentialEquations
using Revise
includet("ISIM_DE_sol.jl")
using Main.SIDE

omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
v=[omega kappa delta epsilon b tau0]

@time SIDE.ISIM_DE(v)
