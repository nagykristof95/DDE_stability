using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using Plots
using TimerOutputs
using DifferentialEquations
pyplot()
using PyPlot
pygui(true);
using Revise
includet("cont_integration_pw.jl")
includet("ISIM_integration_cont.jl")
using Main.ISCPW
using Main.ISCI

omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
v=[omega kappa delta epsilon b tau0]

gmax=20
mult=20
tol=1e-15
ALG0=MethodOfSteps(Tsit5())

sajt=ISIM_CI(v,(tol,gmax,mult,ALG0))
print(BAS.normmax(sajt[:,end]))
sajt=ISIM_CPW(v,(tol,gmax,mult,ALG0))

function test_cont(gmaxtest1,toltest1,multtest1,rep,valref)
    gmaxn=size(gmaxtest1)[1]
    toltestn=size(toltest1)[1]
    multtestn=size(multtest1)[1]
    table=zeros(Float64,multtestn*toltestn*gmaxn,6)

    tabtempcomp=zeros(Float64,rep)
    tabtemperr=zeros(ComplexF64,rep)
    tabtempeval=zeros(Float64,rep)
    reset_timer!(ISCI.toISCI)
    reset_timer!(ISCPW.toISCPW)
    for j1=1:multtestn
        for j2=1:gmaxn
            for j3=1:toltestn
                for j4=1:rep
                    if meth=="CI"
                        valtemp=@timed ISIM_CI(v,(toltest1[j3],gmaxtest1[j2],multtest1[j1],ALG0))
                            tempeval=TimerOutputs.ncalls(ISCI.toISCI["fmult_CI"])
                            reset_timer!(ISCI.toISCI)
                        elseif meth=="CPW"
                            valtemp=@timed ISIM_CPW(v,(toltest1[j3],gmaxtest1[j2],multtest1[j1],ALG0))
                            tempeval=TimerOutputs.ncalls(ISCPW.toISCPW["fmult_PW"])
                            reset_timer!(ISCPW.toISCPW)
                end
                    tabtemperr[j4]=BAS.normmax((valtemp[1])[:,end])
                    tabtempcomp[j4]=valtemp[2]
                    tabtempeval[j4]=tempeval
                end
                vecttempabs=real.(tabtemperr)+im*abs.(imag.(tabtemperr))
                tableerr=abs(mean(vecttempabs-valref*ones(ComplexF64,rep,1)))
                tablecomp=mean(tabtempcomp)
                tableeval=mean(tabtempeval)
                table[(j1-1)*toltestn*gmaxn+(j2-1)*toltestn+j3,:]=[convert(Float64,multtest1[j1]) convert(Float64,gmaxtest1[j2]) convert(Float64,toltest1[j3]) tablecomp tableerr tableeval]
            end
        end
    end
    return(table)
end

#graph parameters
gmaxtest=collect(20:1:20)
toltest=vec([1e-3 1e-5 1e-7 1e-9])
multtest=collect(20:1:20)
rep1=1

valCWn=BAS.normmax(sajt2[:,end])
valCW=real(valCWn)+im*abs(imag.(valCWn))
abs(valrefA-valCW)
print(valCW)
valrefCI= -0.6450472841322661 + 0.8845000929131368*im
valrefCPW=-0.6450472841322712 + 0.8845000929131275*im



meth="CPW"
ALG0=MethodOfSteps(BS3())
BS3_CPW=test_cont(gmaxtest,toltest,multtest,rep1,valrefCPW)
ALG0=MethodOfSteps(RK4())
RK4_CPW=test_cont(gmaxtest,toltest,multtest,rep1,valrefCPW)
ALG0=MethodOfSteps(Tsit5())
T5_CPW=test_cont(gmaxtest,toltest,multtest,rep1,valrefCPW)

meth="CI"
ALG0=MethodOfSteps(BS3())
BS3_CI=test_cont(gmaxtest,toltest,multtest,rep1,valrefCI)
ALG0=MethodOfSteps(RK4())
RK4_CI=test_cont(gmaxtest,toltest,multtest,rep1,valrefCI)
ALG0=MethodOfSteps(Tsit5())
T5_CI=test_cont(gmaxtest,toltest,multtest,rep1,valrefCI)

ALG0=MethodOfSteps(BS3())
BS3_CPW=test_cont(gmaxtest,toltest,multtest,rep1,valrefCPW)
