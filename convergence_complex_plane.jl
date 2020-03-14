using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using Plots
using TimerOutputs
pyplot()
using PyPlot
pygui(true);
using Revise
includet("ISIM_DE_sol.jl")
includet("ISIM_LMS_RK.jl")
includet("ISIM.jl")
using Main.SIDE
using Main.SILMRK
using Main.SI

omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
v=[omega kappa delta epsilon b tau0]


#USUAL PARAMETERS
abstol0=1e-8
ALG0=MethodOfSteps(RK4())
gmax=20
n=1000
mult=20

function convergence(ntest0,rep)
    ns=size(ntest0)[1]
    table=zeros(ComplexF64,ns,2)
    tabletemp=zeros(ComplexF64,rep)
    for j=1:ns
        for jj=1:rep
            if meth=="RK"
                valtemp=BAS.normmax(ISIM(v,(ntest0[j],gmax,mult,ALG0))[:,end])
            elseif meth=="LMS"
                valtemp=BAS.normmax(ISIM_LMS(v,(ntest0[j],gmax,mult,ALG0))[:,end])
            elseif meth=="DE"
                valtemp=BAS.normmax(ISIM_DE(v,(ntest0[j],gmax,mult,ALG0))[:,end])
            end
            valtemp=real(valtemp)+im*abs(imag(valtemp))
            tabletemp[jj]=valtemp
        end
        table[j,:]=[ntest0[j] mean(tabletemp)]
    end
    return(table)
end

convRKl=convergence(ntest,2)
convEEl=convergence(ntest,2)
convLMl=convergence(ntest,2)
refAlast=refA[:,end]

BAS.evplot(refAlast,1)
plot!(real.(hcat(convRKl[:,2],convRKq[:,2],convRKc[:,2],convEEl[:,2],convLMl[:,2],)),imag.(hcat(convRKl[:,2],convRKq[:,2],convRKc[:,2],convEEl[:,2],convLMl[:,2])),markershapes=[:circle],label=["RK4 l" "RK4 q" "RK4 c" "EE" "LM4 l"])
