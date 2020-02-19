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

function test_tol(gmaxtest1,toltest1,multtest1,rep)
    valref=valrefA
    gmaxn=size(gmaxtest1)[1]
    toltestn=size(toltest1)[1]
    multtestn=size(multtest1)[1]
    table=zeros(Float64,multtestn*toltestn*gmaxn,5)

    tabtempcomp=zeros(Float64,rep)
    tabtemperr=zeros(ComplexF64,rep)

    for j1=1:multtestn
        global mult=multtest1[j1]
        for j2=1:gmaxn
            global gmax=gmaxtest1[j2]
            for j3=1:toltestn
                global abstol0=toltest1[j3]
                for j4=1:rep
                    if meth=="RK"
                    valtemp=@timed ISIM(v)
                elseif meth=="LMS"
                    valtemp=@timed ISIM_LSM(v)
                end
                    tabtemperr[j4]=normmax((valtemp[1])[:,end])
                    tabtempcomp[j4]=valtemp[2]
                end
                vecttempabs=real.(tabtemperr)+im*abs.(imag.(tabtemperr))
                tableerr=real(mean(abs.(vecttempabs-valref*ones(ComplexF64,rep,1))))
                tablecomp=mean(tabtempcomp)
                table[(j1-1)*toltestn*gmaxn+(j2-1)*toltestn+j3,:]=[convert(Float64,multtest1[j1]) convert(Float64,gmaxtest1[j2]) convert(Float64,toltest1[j3]) tablecomp tableerr]
            end
        end
    end
    return(table)
end

#DifferentialEquations.jl limited step
meth="RK"
#point A
omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
v=[omega kappa delta epsilon b tau0]
gmax=20
mult=20
n=20000
method="RK"
BR=BRK4
valrefAn=normmax((ISIM(v))[:,end])
valrefA=real(valrefAn)+im*abs(imag.(valrefAn))


gmaxtest=collect(20:1:20)
toltest=[1e-4, 1e-6, 1e-8, 1e-10, 1e-12]
multtest=collect(20:1:20)
n=500

method="Julia"
alg1=MethodOfSteps(BS3())
A_BS3_tol=test_tol(gmaxtest,toltest,multtest,5)
alg1=MethodOfSteps(Tsit5())
A_T5_tol=test_tol(gmaxtest,toltest,multtest,5)

Plots.plot(A_BS3_tol[:,3],hcat(A_BS3_tol[:,4],A_T5_tol[:,4],),xscale=:log10, yscale=:log10,title="Computation time A",label=["BS3 tol" "T5 tol"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_BS3_tol[:,3],hcat(A_BS3_tol[:,5],A_T5_tol[:,5],),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 tol" "T5 tol"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_BS3[:,4],A_T5[:,4],A_RB23[:,4],A_V6[:,4]),hcat(A_BS3[:,5],A_T5[:,5],A_RB23[:,5],A_V6[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")
