using Printf
using Statistics
using LinearAlgebra
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

omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
v=[omega kappa delta epsilon b tau0]



#USUAL PARAMETERS
abstol0=1e-8
ALG0=MethodOfSteps(RK4())
gmax=20
n=500
mult=4
#END OF USUAL PARAMETERS

ISIM_DE(v,(n,gmax,mult,ALG0))

ISIM_LMS(v,(n,gmax,mult,BAS.BLM3))

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
                table[(j1-1)*ntestn*gmaxn+(j2-1)*ntestn+j3,:]=[convert(Float64,multtest1[j1]) convert(Float64,gmaxtest1[j2]) convert(Float64,ntest1[j3]) tablecomp tableerr tabtempeval]
            end
        end
    end
    return(table)
end

#graph parameters
gmaxtest=collect(20:1:20)
ntest=vcat([10,30,50,70,100,150,250],collect(500:1000:8000))
ntest=collect(100:200:1000)
multtest=collect(20:1:20)
rep1=1

#Calculating reference values
ALG0=MethodOfSteps(RK4())
gmax=20
mult=20
n=20000

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
A_EE=test(gmaxtest,ntest,multtest,rep1,valrefA)

Plots.plot(A_RK4[:,3],hcat(A_RK4[:,5],A_RK3[:,5],A_RK2[:,5],A_EE[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4[:,6],A_RK3[:,6],A_RK2[:,6],A_EE[:,6]),hcat(A_RK4[:,5],A_RK3[:,5],A_RK2[:,5],A_EE[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#point B
v=vB
ALG0=BAS.BRK4
B_RK4=test(gmaxtest,ntest,multtest,rep1,valrefB)
ALG0=BAS.BRK3
B_RK3=test(gmaxtest,ntest,multtest,rep1,valrefB)
ALG0=BAS.BRK2
B_RK2=test(gmaxtest,ntest,multtest,rep1,valrefB)
ALG0=BAS.BEE
B_EE=test(gmaxtest,ntest,multtest,rep1,valrefB)

Plots.plot(B_RK4[:,3],hcat(B_RK4[:,5],B_RK3[:,5],B_RK2[:,5],B_EE[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error B",label=["RK4" "RK3" "RK2" "RK1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(B_RK4[:,6],B_RK3[:,6],B_RK2[:,6],B_EE[:,6]),hcat(B_RK4[:,5],B_RK3[:,5],B_RK2[:,5],B_EE[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods B",label=["RK4" "RK3" "RK2" "RK1"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#Linear-Multistep methods

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

Plots.plot(hcat(A_LM4[:,6],A_LM3[:,6],A_LM2[:,6],A_LM1[:,6]),hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")


#point B
v=vB
ALG0=BAS.BLM4
A_LM4=test(gmaxtest,ntest,multtest,rep1,valrefB)
ALG0=BAS.BLM3
A_LM3=test(gmaxtest,ntest,multtest,rep1,valrefB)
ALG0=BAS.BLM2
A_LM2=test(gmaxtest,ntest,multtest,rep1,valrefB)
ALG0=BAS.BLM1
A_LM1=test(gmaxtest,ntest,multtest,rep1,valrefB)

Plots.plot(B_LM4[:,3],hcat(B_LM4[:,5],B_LM3[:,5],B_LM2[:,5],B_LM1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error B",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(B_LM4[:,6],B_LM3[:,6],B_LM2[:,6],B_LM1[:,6]),hcat(B_LM4[:,5],B_LM3[:,5],B_LM2[:,5],B_LM1[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods B",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")


#DifferentialEquations.jl limited step
meth="DE"

#point A
v=vA
ALG0=MethodOfSteps(Euler())
A_EE_J=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=MethodOfSteps(BS3())
A_BS3=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=MethodOfSteps(RK4())
A_RK4_J=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=MethodOfSteps(Tsit5())
A_T5=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=MethodOfSteps(Vern6())
A_V6=test(gmaxtest,ntest,multtest,rep1,valrefA)


Plots.plot(A_EE_J[:,3],hcat(A_EE_J[:,5],A_BS3[:,5],A_RK4_J[:,5],A_T5[:,5],A_V6[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["EE DE" "BS3 DE" "RK4 DE" "T5 DE" "V6 DE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_EE_J[:,6],A_BS3[:,6],A_RK4_J[:,6],A_T5[:,6],A_V6[:,6]),hcat(A_EE_J[:,5],A_BS3[:,5],A_RK4_J[:,5],A_T5[:,5],A_V6[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["EE DE" "BS3 DE" "RK4 DE" "T5 DE" "V6 DE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")


#RK3 & LM4 & BS3 & T5 compare

Plots.plot(A_RK3[:,3],hcat(A_RK3[:,4],A_LM4[:,4],A_BS3[:,4],A_T5[:,4]),xscale=:log10, yscale=:log10,title="Computation time A",label=["RK3" "LM4" "BS3" "T5"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_RK3[:,3],hcat(A_RK3[:,5],A_LM4[:,5],A_BS3[:,5],A_T5[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK3" "LM4" "BS3" "T5"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK3[:,4],A_LM4[:,4],A_BS3[:,4],A_T5[:,4]),hcat(A_RK3[:,5],A_LM4[:,5],A_BS3[:,5],A_T5[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["RK3" "LM4" "BS3" "T5"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")
