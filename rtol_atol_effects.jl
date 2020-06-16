using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using DifferentialEquations
using Plots
using TimerOutputs
pyplot()
using PyPlot
pygui(true);
using Revise
includet("ISIM_DE_sol.jl")
using Main.SIDE

setprecision(BigFloat,128)
setrounding(BigFloat, RoundUp)

prec=Float64
precF=ComplexF64

# prec=BigFloat
# precF=Complex{BigFloat}

omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
v=prec.([omega kappa delta epsilon b tau0])
vA=v

function test_tol(ntest1,rtoltest1,atoltest1,rep,valref)
    ntestn=size(ntest1)[1]
    rtoltestn=size(rtoltest1)[1]
    atoltestn=size(atoltest1)[1]
    table0=zeros(Float64,rtoltestn*atoltestn*ntestn,6)

    tabtempcomp=zeros(Float64,rep)
    tabtemperr=zeros(ComplexF64,rep)
    tabtempeval=zeros(Float64,rep)
    reset_timer!(SIDE.toSIDE)
    for j1=1:rtoltestn
        for j2=1:atoltestn
            for j3=1:ntestn
                for j4=1:rep
                    if meth=="DE"
                        #valtemp=@timed ISIM_DE(v,(ntest1[j3],gmax,mult,ALG0),(rtoltest1[j1],atoltest1[j2]))
                        valtemp=@timed ISIM_DE(v,(ntest1[j3],gmax,mult,ALG0),(atoltest1[j2],atoltest1[j2]))
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
                tableeval=mult*mean(tabtempeval)
                table0[(j1-1)*ntestn*atoltestn+(j2-1)*ntestn+j3,:]=[convert(Float64,rtoltest1[j1]) convert(Float64,atoltest1[j2]) convert(Float64,ntest1[j3]) tablecomp tableerr tableeval]
            end
        end
    end
    return(table0)
end

ALG0=MethodOfSteps(RK4())
gmax=12
n=2
mult=12
nlim=20000

rep1=2

rtoltest=[1e-1; 1e-2; 1e-3; 1e-4; 1e-5; 1e-6; 1e-7; 1e-8; 1e-9; 1e-10; 1e-11; 1e-12]
atoltest=[1e-1; 1e-2; 1e-3; 1e-4; 1e-5; 1e-6; 1e-7; 1e-8; 1e-9; 1e-10; 1e-11; 1e-12]
ntest=collect(10:15:300)
ntest=collect(10:250:3000)
ntest=[2; 3; 4; 5; 6; 7; 8; 9; 10; 25; 50; 100; 200; 300; 400; 500; 1000; 1500; 2000; 2500; 3000]
ntest=[100]

#point A
omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
v=prec.([omega kappa delta epsilon b tau0])
vA=v
valrefA_DEn=BAS.normmax((ISIM_DE(vA,(nlim,gmax,mult,ALG0),(1e1,1e1)))[:,end])
valrefA_DE=real(valrefA_DEn)+im*abs(imag.(valrefA_DEn))
valrefA_DE=-0.6450472841316337 + 0.8845000929140697*im
valrefA_DE_stf=-0.6899534713343562 + 0.8551923368854594*im
valrefA_DE_sgn=-0.6896571878614839 + 0.8555411241268935*im
valrefA_DE_sgn_ad=-0.6896644288107385 + 0.8555322147671132*im
#point B
omega=1; kappa=0.2; delta=3.246; epsilon=1.0; b=-0.81; tau0=2*pi;
v=prec.([omega kappa delta epsilon b tau0])
vB=v
valrefB_DEn=BAS.normmax((ISIM_DE(vB,(nlim,gmax,mult,ALG0),(1e1,1e1)))[:,end])
valrefB_DE=real(valrefB_DEn)+im*abs(imag.(valrefB_DEn))
valrefB_DE=0.778812710770659 + 0.2160496901859506*im
valrefB_DE_stf=-0.8150668231207973 + 0.14036267365239075*im

#heatmap
meth="DE"
ALG0=MethodOfSteps(RK4())
A_RK4_DE_ad=test_tol(ntest,rtoltest,atoltest,rep1,valrefA_DE)
B_RK4_DE_ad=test_tol(ntest,rtoltest,atoltest,rep1,valrefB_DE)

Plots.contour(log10.(A_RK4_DE_ad[:,1]),log10.(A_RK4_DE_ad[:,2]),log10.(A_RK4_DE_ad[:,4]),st=:surface,color=:inferno,camera=(0,90))

Plots.contour(log10.(A_RK4_DE_ad[:,1]),log10.(A_RK4_DE_ad[:,2]),log10.(A_RK4_DE_ad[:,5]),st=:surface,color=:inferno,camera=(0,90))

#reltol fixed
meth="DE"
ntest=[100]
rtoltest=[1e-6]
atoltest=[1e-6]
atoltest=[1e-5; 1e-6; 1e-7; 1e-8; 1e-9; 1e-10]

#point A
probaA=ISIM_DE(vA,(5,gmax,mult,ALG0),(1e-8,1e-8))
println(proba)
ALG0=MethodOfSteps(BS3())
A_BS3_DE_ad=test_tol(ntest,rtoltest,atoltest,rep1,valref)
ALG0=MethodOfSteps(RK4())
A_RK4_DE_ad=test_tol(ntest,rtoltest,atoltest,rep1,valrefA_DE)
ALG0=MethodOfSteps(Tsit5())
A_T5_DE_ad=test_tol(ntest,rtoltest,atoltest,rep1,valref)

Plots.plot(real.(probaA[:,6]))
Plots.plot(imag.(probaB[:,11]))
#point A
probaB=ISIM_DE(vB,(2,gmax,mult,ALG0),(1e-8,1e-8))
ALG0=MethodOfSteps(BS3())
B_BS3_DE_ad=test_tol(ntest,rtoltest,atoltest,rep1,valrefB_DE)

open("temp.txt", "w") do io
           writedlm(io, A_T5_J_dense_BF128)
       end

Plots.savefig("1.pdf")

sajt=ISIM_DE(vA,(20000,gmax,mult,ALG0),(1e-15,1e-15))
valrefn=BAS.normmax(sajt[:,end])

#valref msapmpling=20000 and reltol=abstol=1e-15
valref=precF(-0.6896644288103551 + 0.8555322147675607*im)


ALG0=MethodOfSteps(RK4())
A_RK4_DE_ad=test_tol(ntest,rtoltest,atoltest,rep1,valrefA_DE)
