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
rng = MersenneTwister(1234)

valrefA_lin=-0.645047283999984+0.884500092450041*im
valrefA_quad=-0.645047338561817+0.884500265852046*im

#graph parameters
gmaxtest=collect(20:1:20)
ntest=vcat([10,30,50,70,100,150,250],collect(500:1000:8000))
multtest=collect(20:1:20)
rep1=3

#Runge-Kutta methods
meth="RK"
method="RK"
#point A LIN
inttyp0=BSpline(Linear())
BR=BRK4
A_RK4_lin=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)
BR=BRK3
A_RK3_lin=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)
BR=BRK2
A_RK2_lin=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)
BR=BEE
A_EE_lin=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)

#point A QUAD
inttyp0=BSpline(Quadratic(Line(OnGrid())))
BR=BRK4
A_RK4_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)
BR=BRK3
A_RK3_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)
BR=BRK2
A_RK2_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)
BR=BEE
A_EE_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)

#Multistep methods
meth="LMS"
#point A LIN
inttyp0=BSpline(Linear())
BLM=BLM4
A_LM4_lin=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)
BLM=BLM3
A_LM3_lin=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)
BLM=BLM2
A_LM2_lin=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)
BLM=BLM1
A_LM1_lin=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)

#point A QUAD
inttyp0=BSpline(Quadratic(Line(OnGrid())))
BLM=BLM4
A_LM4_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)
BLM=BLM3
A_LM3_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)
BLM=BLM2
A_LM2_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)
BLM=BLM1
A_LM1_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)

#Julia
meth="RK"
method="Julia"
#point A LIN
inttyp0=BSpline(Linear())
alg1=MethodOfSteps(BS3())
A_BS3_l=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)
alg1=MethodOfSteps(Tsit5())
A_T5_l=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)
alg1=MethodOfSteps(Vern6())
A_V6_l=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)

#point A QUAD
inttyp0=BSpline(Quadratic(Line(OnGrid())))
alg1=MethodOfSteps(BS3())
A_BS3_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)
alg1=MethodOfSteps(Tsit5())
A_T5_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)
alg1=MethodOfSteps(Vern6())
A_V6_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)






# # # RK plots # # #
#linear
Plots.plot(A_RK4_lin[:,3],hcat(A_RK4_lin[:,4],A_RK3_lin[:,4],A_RK2_lin[:,4],A_EE_lin[:,4]),xscale=:log10, yscale=:log10,title="Computation time A (linear ip.)",label=["RK4" "RK3" "RK2" "RK1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_RK4_lin[:,3],hcat(A_RK4_lin[:,5],A_RK3_lin[:,5],A_RK2_lin[:,5],A_EE_lin[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A (linear ip.)",label=["RK4" "RK3" "RK2" "RK1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4_lin[:,4],A_RK3_lin[:,4],A_RK2_lin[:,4],A_EE_lin[:,4]),hcat(A_RK4_lin[:,5],A_RK3_lin[:,5],A_RK2_lin[:,5],A_EE_lin[:,5]),xscale=:log10, yscale=:log10,title="Computation time A (linear ip.)",label=["RK4" "RK3" "RK2" "RK1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

#quadratic
Plots.plot(A_RK4_q[:,3],hcat(A_RK4_q[:,4],A_RK3_q[:,4],A_RK2_q[:,4],A_EE_q[:,4]),xscale=:log10, yscale=:log10,title="Computation time A (quadratic ip.)",label=["RK4" "RK3" "RK2" "RK1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_RK4_q[:,3],hcat(A_RK4_q[:,5],A_RK3_q[:,5],A_RK2_q[:,5],A_EE_q[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A (quadratic ip.)",label=["RK4" "RK3" "RK2" "RK1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4_q[:,4],A_RK3_q[:,4],A_RK2_q[:,4],A_EE_q[:,4]),hcat(A_RK4_q[:,5],A_RK3_q[:,5],A_RK2_q[:,5],A_EE_q[:,5]),xscale=:log10, yscale=:log10,title="Computation time A (quadratic ip.)",label=["RK4" "RK3" "RK2" "RK1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

# # # LMS plots # # #
#linear
Plots.plot(A_LM4_lin[:,3],hcat(A_LM4_lin[:,4],A_LM3_lin[:,4],A_LM2_lin[:,4],A_LM1_lin[:,4]),xscale=:log10, yscale=:log10,title="Computation time A (linear ip.)",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_LM4_lin[:,3],hcat(A_LM4_lin[:,5],A_LM3_lin[:,5],A_LM2_lin[:,5],A_LM1_lin[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A (linear ip.)",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_LM4_lin[:,4],A_LM3_lin[:,4],A_LM2_lin[:,4],A_EE_lin[:,4]),hcat(A_LM4_lin[:,5],A_LM3_lin[:,5],A_LM2_lin[:,5],A_EE_lin[:,5]),xscale=:log10, yscale=:log10,title="Computation time A (linear ip.)",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

#quadratic
Plots.plot(A_LM4_q[:,3],hcat(A_LM4_q[:,4],A_LM3_q[:,4],A_LM2_q[:,4],A_LM1_q[:,4]),xscale=:log10, yscale=:log10,title="Computation time A (quadratic ip.)",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_LM4_q[:,3],hcat(A_LM4_q[:,5],A_LM3_q[:,5],A_LM2_q[:,5],A_LM1_q[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A (quadratic ip.)",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_LM4_q[:,4],A_LM3_q[:,4],A_LM2_q[:,4],A_EE_q[:,4]),hcat(A_LM4_q[:,5],A_LM3_q[:,5],A_LM2_q[:,5],A_LM1_q[:,5]),xscale=:log10, yscale=:log10,title="Computation time A (quadratic ip.)",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

# # # Julia plots # # #
#linear
Plots.plot(A_V6_l[:,3],hcat(A_V6_l[:,4],A_T5_l[:,4],A_BS3_l[:,4]),xscale=:log10, yscale=:log10,title="Computation time A (linear ip.)",label=["V6" "T5" "BS3"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_V6_l[:,3],hcat(A_V6_l[:,5],A_T5_l[:,5],A_BS3_l[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A (linear ip.)",label=["V6" "T5" "BS3"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_V6_l[:,4],A_T5_l[:,4],A_BS3_l[:,4]),hcat(A_V6_l[:,5],A_T5_l[:,5],A_BS3_l[:,5]),xscale=:log10, yscale=:log10,title="Computation time A (linear ip.)",label=["V6" "T5" "BS3"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

#quadratic
Plots.plot(A_V6_q[:,3],hcat(A_V6_q[:,4],A_T5_q[:,4],A_BS3_q[:,4]),xscale=:log10, yscale=:log10,title="Computation time A (quadratic ip.)",label=["V6" "T5" "BS3"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_V6_q[:,3],hcat(A_V6_q[:,5],A_T5_q[:,5],A_BS3_q[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A (quadratic ip.)",label=["V6" "T5" "BS3"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_V6_q[:,4],A_T5_q[:,4],A_BS3_q[:,4]),hcat(A_V6_q[:,5],A_T5_q[:,5],A_BS3_q[:,5]),xscale=:log10, yscale=:log10,title="Computation time A (quadratic ip.)",label=["V6" "T5" "BS3"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Atest=transpose([0.0 0.5 1.0 1.5 2.0 2.5; 1.3 2.3 -5.0 -3.1 -0.5 0.3])
inttyp0=BSpline(Linear())
it_l=it(Atest)
inttyp0=BSpline(Quadratic(Line(OnGrid())))
it_q=it(Atest)

it_lf(x)=sub(it_l,x)[1]
it_qf(x)=sub(it_q,x)[1]
Plots.plot(it_lf,0.0,2.5,label="linear")
it_lf(x)=sub(it_l,x)[1]
Plots.plot!(it_qf,0.0,2.5,label="quadratic")
