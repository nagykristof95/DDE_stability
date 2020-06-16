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
using SparseArrays
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


setprecision(BigFloat,128)
setrounding(BigFloat, RoundUp)

prec=Float64
precF=ComplexF64

# prec=BigFloat
# precF=Complex{BigFloat}

omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
v=prec.([omega kappa delta epsilon b tau0])
vA=v
#USUAL PARAMETERS
abstol0=1e-8
ALG0=MethodOfSteps(RK4())
ALG0=BAS.BRK4
gmax=12
n=100
mult=12

reltol0=1e-5
abstol0=1e-5
ISIM(v,(n,gmax,mult,ALG0))
sajt=ISIM_DE(vA,(n,gmax,mult,ALG0),(reltol0,abstol0))
tn=size(sajt.t)[1]
Plots.scatter(sajt.t,zeros(Float64,tn))




valrefBn=BAS.normmax(SDM_full(vB,(n,gmax,mult,"Full"))[:,end])
valrefB=real(valrefBn)+im*abs(imag(valrefBn))

@time BAS.normmax(SDM_full(vB,(n,gmax,mult,"Full"))[:,end])
@time BAS.normmax(SDM_full_app(v,(n,gmax,mult))[:,end])
@time BAS.normmax(SDM_full_app(vB,(n,gmax,mult))[:,end])

#END OF USUAL PARAMETERS

function test(gmaxtest1,ntest1,multtest1,rep,valref)
    gmaxn=size(gmaxtest1)[1]
    ntestn=size(ntest1)[1]
    multtestn=size(multtest1)[1]
    table=zeros(prec,multtestn*ntestn*gmaxn,6)

    tabtempcomp=zeros(prec,rep)
    tabtemperr=zeros(precF,rep)
    tabtempeval=zeros(prec,rep)
    reset_timer!(BAS.to)
    reset_timer!(SI.toISIM)
    reset_timer!(SILMRK.toSILMRK)
    reset_timer!(SIDE.toSIDE)
    reset_timer!(SD.toSDM)
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
                            valtemp=@timed ISIM_LMS_RK(v,(ntest1[j3],gmaxtest1[j2],multtest1[j1],ALG0))
                            tempeval=TimerOutputs.ncalls(SILMRK.toSILMRK["fmult_LM"])
                            reset_timer!(SILMRK.toSILMRK)
                        elseif meth=="DE"
                            valtemp=@timed ISIM_DE(v,(ntest1[j3],gmaxtest1[j2],multtest1[j1],ALG0),(reltol0,abstol0))
                            tempeval=TimerOutputs.ncalls(SIDE.toSIDE["fmult_DE"])
                            reset_timer!(SIDE.toSIDE)
                        elseif meth=="DEA"
                            valtemp=@timed ISIM_DE(v,(n,gmaxtest1[j2],multtest1[j1],ALG0),(ntest1[j3],ntest1[j3]))
                            tempeval=TimerOutputs.ncalls(SIDE.toSIDE["fmult_DE"])
                            reset_timer!(SIDE.toSIDE)
                        elseif meth=="SD"
                            valtemp=@timed SDM(v,ntest1[j3])
                            tempeval=TimerOutputs.ncalls(SD.toSDM["f_SDM"])
                            reset_timer!(SD.toSDM)
                        elseif meth=="SDsp"
                            valtemp=@timed SDM_sp(v,ntest1[j3])
                            tempeval=TimerOutputs.ncalls(SD.toSDM["f_SDM"])
                            reset_timer!(SD.toSDM)
                        elseif meth=="SDfull"
                            valtemp=@timed SDM_full(v,(ntest1[j3],gmaxtest1[j2],multtest1[j1],ALG0))
                            tempeval=TimerOutputs.ncalls(SD.toSDM["f_SDM"])
                            reset_timer!(SD.toSDM)
                end
                    tabtemperr[j4]=BAS.normmax((valtemp[1])[:,end])
                    tabtempcomp[j4]=valtemp[2]
                    tabtempeval[j4]=tempeval
                end
                vecttempabs=real.(tabtemperr)+im*abs.(imag.(tabtemperr))
                tableerr=real(abs(mean(vecttempabs-valref*ones(precF,rep,1))))
                tablecomp=mean(tabtempcomp)
                tableeval=mean(tabtempeval)
                table[(j1-1)*ntestn*gmaxn+(j2-1)*ntestn+j3,:]=[convert(prec,multtest1[j1]) convert(prec,gmaxtest1[j2]) convert(prec,ntest1[j3]) tablecomp tableerr tableeval]
            end
        end
    end
    return(table)
end

#graph parameters
gmaxtest=collect(12:1:12)
ntest=collect(100:100:1500)
multtest=collect(12:1:12)
#ntest=collect(100:30:300)

rep1=2

meth="SD"
SDtest=test(gmaxtest,ntest,multtest,rep1,valrefA)
meth="SDsp"
SDsptest=test(gmaxtest,ntest,multtest,rep1,valrefA)
meth="SDfull"
SDfulltest=test(gmaxtest,ntest,multtest,rep1,valrefA)

meth="SDfull"
ALG0="Full"
SDtest=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0="ISI"
SDISI_1=test(collect(5:1:5),ntest,collect(4:1:4),rep1,valrefA)
SDISI_2=test(collect(10:1:10),ntest,collect(8:1:8),rep1,valrefA)

#Calculating reference values
ALG0=BAS.BRK4
gmax=12
mult=12
nlim=20000

#point A
omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
vA=[omega kappa delta epsilon b tau0]

ALG0=BAS.BLM4
valrefAn=BAS.normmax(ISIM_LMS_RK(v,(nlim,gmax,mult,ALG0))[:,end])

valrefAn=BAS.normmax(ISIM(v,(nlim,gmax,mult,ALG0))[:,end])
valrefA=real(valrefAn)+im*abs(imag.(valrefAn))

valrefA_lin=-0.645047283999984+0.884500092450041*im
valrefA_quad=-0.645047338561817+0.884500265852046*im
valrefA=-0.645047283999984+0.884500092450041*im
valrefAad=-0.6450472841286765+0.8845000929183954*im
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
A_RK4=test(gmaxtest,ntest,multtest,rep1,valref)
ALG0=BAS.BRK3
A_RK3=test(gmaxtest,ntest,multtest,rep1,valref)
ALG0=BAS.BRK2
A_RK2=test(gmaxtest,ntest,multtest,rep1,valref)
ALG0=BAS.BEE
A_RK1=test(gmaxtest,ntest,multtest,rep1,valref)

Plots.plot(A_RK4[:,3],hcat(A_RK4[:,5],A_RK3[:,5],A_RK2[:,5],A_RK1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4[:,6],A_RK3[:,6],A_RK2[:,6],A_RK1[:,6]),hcat(A_RK4[:,5],A_RK3[:,5],A_RK2[:,5],A_RK1[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#Linear-Multistep methods

meth="LMS"
#point A
ALG0=BAS.BLM4
valrefALMn=BAS.normmax((ISIM_LMS_RK(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefALM=real(valrefALMn)+im*abs(imag.(valrefALMn))
ALG0=BAS.BLM4
A_LM4=test(gmaxtest,ntest,multtest,rep1,valref)
ALG0=BAS.BLM3
A_LM3=test(gmaxtest,ntest,multtest,rep1,valref)
ALG0=BAS.BLM2
A_LM2=test(gmaxtest,ntest,multtest,rep1,valref)
ALG0=BAS.BLM1
A_LM1=test(gmaxtest,ntest,multtest,rep1,valref)

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_LM4[:,6],A_LM3[:,6],A_LM2[:,6],A_LM1[:,6]),hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#DifferentialEquations.jl limited step
meth="RK"

#point A
v=vA
ALG0=MethodOfSteps(RK4())
valrefARK4n=BAS.normmax((ISIM(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefARK4=real(valrefARK4n)+im*abs(imag.(valrefARK4n))

ALG0=MethodOfSteps(Euler())
A_EE_J=test(gmaxtest,ntest,multtest,rep1,valrefARK4)
ALG0=MethodOfSteps(BS3())
A_BS3=test(gmaxtest,ntest,multtest,rep1,valrefARK4)
ALG0=MethodOfSteps(RK4())
A_RK4_J=test(gmaxtest,ntest,multtest,rep1,valrefARK4)
ALG0=MethodOfSteps(Tsit5())
A_T5=test(gmaxtest,ntest,multtest,rep1,valrefARK4)


Plots.plot(A_EE_J[:,3],hcat(A_EE_J[:,5],A_BS3[:,5],A_RK4_J[:,5],A_T5[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["EE DE" "BS3 DE" "RK4 DE" "T5 DE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_EE_J[:,6],A_BS3[:,6],A_RK4_J[:,6],A_T5[:,6]),hcat(A_EE_J[:,5],A_BS3[:,5],A_RK4_J[:,5],A_T5[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["EE DE" "BS3 DE" "RK4 DE" "T5 DE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#DifferentialEquations.jl adaptive step
meth="DE"

#point A
v=vA
ALG0=MethodOfSteps(RK4())
valrefARK4Dn=BAS.normmax((ISIM_DE(vA,(nlim,gmax,mult,ALG0),(reltol0,abstol0)))[:,end])
valrefARK4D=real(valrefARK4Dn)+im*abs(imag(valrefARK4Dn))
valrefARK4D=valrefARK4
ALG0=MethodOfSteps(Euler())
A_EE_J=test(gmaxtest,ntest,multtest,rep1,valref)
ALG0=MethodOfSteps(BS3())
A_BS3=test(gmaxtest,ntest,multtest,rep1,valref)
ALG0=MethodOfSteps(RK4())
A_RK4_J=test(gmaxtest,ntest,multtest,rep1,valref)
ALG0=MethodOfSteps(Tsit5())
A_T5=test(gmaxtest,ntest,multtest,rep1,valrefARK4D)

Plots.plot(A_EE_J_ad[:,3],hcat(A_EE_J_ad[:,5],A_BS3_ad[:,5],A_RK4_J_ad[:,5],A_T5_ad[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["EE DE" "BS3 DE" "RK4 DE" "T5 DE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_EE_J_ad[:,6],A_BS3_ad[:,6],A_RK4_J_ad[:,6],A_T5_ad[:,6]),hcat(A_EE_J_ad[:,5],A_BS3_ad[:,5],A_RK4_J_ad[:,5],A_T5_ad[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["EE DE" "BS3 DE" "RK4 DE" "T5 DE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")



#RK3 & LM4 & BS3 & T5 compare

Plots.plot(A_RK3[:,3],hcat(A_RK3[:,4],A_LM4[:,4],A_BS3[:,4],A_T5[:,4]),xscale=:log10, yscale=:log10,title="Computation time A",label=["RK3" "LM4" "BS3" "T5"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_RK3[:,3],hcat(A_RK3[:,5],A_LM4[:,5],A_BS3[:,5],A_T5[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK3" "LM4" "BS3" "T5"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK3[:,4],A_LM4[:,4],A_BS3[:,4],A_T5[:,4]),hcat(A_RK3[:,5],A_LM4[:,5],A_BS3[:,5],A_T5[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["RK3" "LM4" "BS3" "T5"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#ISI numsim comparison
#reference RK4 cub, mref=20000
nlim=20000
vA=v
valrefAISIn=BAS.normmax((ISIM(vA,(nlim,gmax,mult,BAS.BRK4)))[:,end])
valrefAISI=real(valrefAISIn)+im*abs(imag(valrefAISIn))

#ranges
gmaxtest=collect(12:1:12)
ntest=collect(100:100:1500)
multtest=collect(12:1:12)
rep1=2


#RK4 lin
meth="RK"
ALG0=BAS.BRK4
A_RK4_l=test(gmaxtest,ntest,multtest,rep1,valrefAISI)

#RK4 cub
meth="RK"
ALG0=BAS.BRK4
A_RK4_c=test(gmaxtest,ntest,multtest,rep1,valrefAISI)

#AB4
#RK4 cub
meth="LMS"
ALG0=BAS.BLM4
A_AB4=test(gmaxtest,ntest,multtest,rep1,valrefAISI)

#RK4*
meth="DE"
ALG0=MethodOfSteps(RK4())
A_RK4_dense=test(gmaxtest,ntest,multtest,rep1,valref)
rep1=5
#RK4* adaptive
meth="DEA"
n=100
toltest=[10^(-4.5); 10^(-5); 10^(-5.5); 10^(-6); 10^(-6.5); 10^(-7); 10^(-7.5); 10^(-8);  10^(-8.5); 10^(-9); 10^(-9.5); 10^(-10);  10^(-10.5);  10^(-11); 10^(-11.5)]
ALG0=MethodOfSteps(Tsit5())
A_T5_dense_ad=test(gmaxtest,rtoltest,multtest,rep1,valref)

valrefn=BAS.normmax(ISIM_DE(v,(nlim,gmax,mult,ALG0),(reltol0,abstol0))[:,end])
valref=real(valrefn)+im*abs(imag(valrefn))
ALG0=MethodOfSteps(RK4())
reltol0=1e-4
abstol0=1e-4

meth="DE"
ALG0=MethodOfSteps(RK4())
A_RK4_dense_fix=test(gmaxtest,ntest,multtest,rep1,valref)


###### MM ######
# #SYSTEM PARAMETERS
mx=0.01986;my=0.02008;kx=1.60312;ky=1.155697;sx=408866;sy=413445; #structural parameters
kt=644*10^6;kr=0.368; #cutting parameters
d=0.008;ae=0.0004;z=1;strat=1.0; #tool/strategy parameters
#
#
omega=4000 #spindle speed
w=0.001 #axial immersion

# #v=[1 mx  2 my  3 kx  4 ky  5 sx  6 sy  7 kt  8 kr  9 d  10 ae  11 z  12 strat  13 omega  14 w ]
v=[mx my kx ky sx sy kt kr d ae z strat omega w]

n=100
gmax=12
mult=12

reltol0=1e-14
abstol0=1e-14
ALG0=MethodOfSteps(RK4())
sajt=ISIM_DE(v,(n,gmax,mult,ALG0),(reltol0,abstol0))

BAS.normmax(sajt[:,end])

#z1
#valrefZ1=0.6144450209445577 + 0.7019749771732842*im
valrefZ1=0.6144530338830984 + 0.7019781358042942*im #adpat 1e-14 (EZT!)
#valrefZ1=0.6144530339251405 + 0.7019781357588984*im

valrefZ3=0.5026384325865356 + 0.7478515394410706*im

meth="DE"

gmaxtest=collect(12:1:12)
ntest=collect(100:100:1500)
multtest=collect(12:1:12)

rep1=5
Z1_RK4_dense=test(gmaxtest,ntest,multtest,rep1,valrefZ1)
Z3_RK4_dense=test(gmaxtest,ntest,multtest,rep1,valrefZ3)

meth="DEA"
Z1_RK4_dense_ad=test(gmaxtest,rtoltest,multtest,rep1,valrefZ1)
Z3_RK4_dense_ad=test(gmaxtest,rtoltest,multtest,rep1,valrefZ3)
