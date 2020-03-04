using Printf
using Statistics
using LinearAlgebra
using GenericSchur
using Interpolations
using DelimitedFiles
using Random
using DifferentialEquations
using Plots
using TimerOutputs
pyplot()
using PyPlot
pygui(true);
using Revise
includet("ISIM_DE_sol.jl")
includet("ISIM_LMS.jl")
includet("ISIM.jl")
using Main.SIDE
using Main.SILMS
using Main.SI

omega=1.0; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
v=[omega kappa delta epsilon b tau0]

#USUAL PARAMETERS
abstol0=1e-8
ALG0=MethodOfSteps(RK4())
gmax=6
n=500
mult=8
#END OF USUAL PARAMETERS

ISIM_DE(v,(n,gmax,mult,ALG0))

ISIM_LMS(v,(n,gmax,mult,BAS.BLM4))

ISIM(v,(n,gmax,mult,BAS.BRK4))

ISIM(v,(n,gmax,mult,ALG0))

function test(gmaxtest1,ntest1,multtest1,rep,valref)
    gmaxn=size(gmaxtest1)[1]
    ntestn=size(ntest1)[1]
    multtestn=size(multtest1)[1]
    table=zeros(Float64,multtestn*ntestn*gmaxn,6)

    tabtempcomp=zeros(Float64,rep)
    tabtemperr=zeros(ComplexF64,rep)
    tabtempeval=zeros(Float64,rep)
    reset_timer!(BAS.to)
    reset_timer!(SI.toISIM)
    reset_timer!(SILMS.toSILMS)
    reset_timer!(SIDE.toSIDE)
    for j1=1:multtestn
        for j2=1:gmaxn
            for j3=1:ntestn
                for j4=1:rep
                    if meth=="RK"
                        valtemp=@timed ISIM(v,(ntest1[j3],gmaxtest1[j2],multtest1[j1],ALG0))
                        if typeof(ALG0) <: MethodOfSteps
                            tempeval=TimerOutputs.ncalls(SI.toISIM["fmult_Julia"])
                            reset_timer!(SI.toISIM)
                        elseif typeof(ALG0) <: Tuple
                            tempeval=TimerOutputs.ncalls(BAS.to["fmult"])
                            reset_timer!(BAS.to)
                        end

                        elseif meth=="LMS"
                            valtemp=@timed ISIM_LMS(v,(ntest1[j3],gmaxtest1[j2],multtest1[j1],ALG0))
                            tempeval=TimerOutputs.ncalls(SILMS.toSILMS["fmult_LMS"])
                            reset_timer!(SILMS.toSILMS)

                        elseif meth=="DE"
                            valtemp=@timed ISIM_DE(v,(ntest1[j3],gmaxtest1[j2],multtest1[j1],ALG0))
                            tempeval=TimerOutputs.ncalls(SIDE.toSIDE["fmult_DE"])
                            reset_timer!(SIDE.toSIDE)
                end
                    tabtemperr[j4]=BAS.normmax((valtemp[1])[:,end])
                    tabtempcomp[j4]=valtemp[2]
                    tabtempeval[j4]=tempeval
                end
                vecttempabs=real.(tabtemperr)+im*abs.(imag.(tabtemperr))
                tableerr=real(abs(mean(vecttempabs-valref*ones(ComplexF64,rep,1))))
                tablecomp=mean(tabtempcomp)
                tableeval=mean(tabtempeval)
                table[(j1-1)*ntestn*gmaxn+(j2-1)*ntestn+j3,:]=[convert(Float64,multtest1[j1]) convert(Float64,gmaxtest1[j2]) convert(Float64,ntest1[j3]) tablecomp tableerr tableeval]
            end
        end
    end
    return(table)
end

#graph parameters
gmaxtest=collect(20:1:20)
ntest=vcat([10,30,50,70,100,150,250],collect(500:1000:4000))
ntest=collect(100:200:1000)
ntest=collect(20:15:95)
multtest=collect(20:1:20)
rep1=1

#Calculating reference values
ALG0=MethodOfSteps(RK4())
gmax=20
mult=20
n=100000

#point A
omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
vA=[omega kappa delta epsilon b tau0]
valrefAn=BAS.normmax((ISIM(v,(n,gmax,mult,BAS.BRK4)))[:,end])
valrefA=real(valrefAn)+im*abs(imag.(valrefAn))
valrefA_lin=-0.645047283999984+0.884500092450041*im
valrefA_quad=-0.645047338561817+0.884500265852046*im
valrefA=-0.645047283999984+0.884500092450041*im

#point B
omega=1; kappa=0.2; delta=3.246; epsilon=1.0; b=-0.81; tau0=2*pi;
vB=[omega kappa delta epsilon b tau0]
valrefBn=BAS.normmax((ISIM(v,(n,gmax,mult,BAS.BRK4)))[:,end])
valrefB=real(valrefBn)+im*abs(imag.(valrefBn))

################################################################################
#Runge-Kutta methods
meth="RK"
#pointA LIN
v=vA
ALG0=BAS.BRK4
A_RK4=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=BAS.BRK3
A_RK3=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=BAS.BRK2
A_RK2=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=BAS.BEE
A_RK1=test(gmaxtest,ntest,multtest,rep1,valrefA)

Plots.plot(A_RK4[:,3],hcat(A_RK4[:,5],A_RK3[:,5],A_RK2[:,5],A_RK1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK4" "RK3" "RK2" "RK1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")
Plots.plot(hcat(A_RK4[:,6],A_RK3[:,6],A_RK2[:,6],A_RK1[:,6]),hcat(A_RK4[:,5],A_RK3[:,5],A_RK2[:,5],A_RK1[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["RK4" "RK3" "RK2" "RK1"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")


meth="LMS"
#point A
v=vA
ALG0=BAS.BLM4
A_LM4=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=BAS.BLM3
A_LM3=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=BAS.BLM2
A_LM2=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=BAS.BLM1
A_LM1=test(gmaxtest,ntest,multtest,rep1,valrefA)

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")
Plots.plot(hcat(A_LM4[:,6],A_LM3[:,6],A_LM2_J[:,6],A_LM1[:,6]),hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")


#DifferentialEquations.jl limited step + internal interpolation
meth="DE"

#point A
v=vA
ALG0=MethodOfSteps(Euler())
A_EE_J_=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=MethodOfSteps(BS3())
A_BS3_=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=MethodOfSteps(RK4())
A_RK4_J_=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=MethodOfSteps(Tsit5())
A_T5_=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=MethodOfSteps(Vern6())
A_V6_=test(gmaxtest,ntest,multtest,rep1,valrefA)

Plots.plot(A_EE_Jq[:,3],hcat(A_EE_Jq[:,5],A_BS3q[:,5],A_RK4_Jq[:,5],A_T5q[:,5],A_V6q[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["EE DE q" "BS3 DE q" "RK4 DE q" "T5 DE q" "V6 DE q"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")
Plots.plot(hcat(A_EE_Jq[:,6],A_BS3q[:,6],A_RK4_Jq[:,6],A_T5q[:,6],A_V6q[:,6]),hcat(A_EE_Jq[:,5],A_BS3q[:,5],A_RK4_Jq[:,5],A_T5q[:,5],A_V6q[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["EE DE" "BS3 DE" "RK4 DE" "T5 DE" "V6 DE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")
