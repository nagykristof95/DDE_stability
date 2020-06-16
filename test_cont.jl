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
includet("cont_integration_pw.jl")
includet("ISIM_CIN.jl")
using Main.ISCPW
using Main.ISCIN



setprecision(BigFloat,128)
setrounding(BigFloat, RoundUp)

prec=Float64
precF=ComplexF64

omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
v=prec.([omega kappa delta epsilon b tau0])
vA=v
#USUAL PARAMETERS

abstol0=1e-12
reltol0=1e-12
abstol0int=1e-5
reltol0int=1e-5

ALG0=MethodOfSteps(RK4())
gmax=12
n=100
mult=12

refCPW=ISIM_CPW(v,(n,gmax,mult,ALG0),(reltol0,abstol0),(reltol0int,abstol0int))

refCIN=ISIM_CIN(v,(n,gmax,mult,ALG0),(reltol0,abstol0))
valrefCINn=BAS.normmax(refCIN[:,end])
valrefCIN=real(valrefCINn)+abs(imag(valrefCINn))*im
valrefCIN_S=-0.6450472841323589 + 0.8845000929129383*im
valrefCIN_NS=-0.6896644288104126 + 0.8555322147674023*im

function test_cont(ntest1,tolabstest1,tolreltest1,rep,valref)
    ntestN=size(ntest1)[1]
    tolabstestN=size(tolabstest1)[1]
    tolreltestN=size(tolreltest1)[1]
    table=zeros(prec,ntestN*tolabstestN*tolreltestN,6)

    tabtempcomp=zeros(prec,rep)
    tabtemperr=zeros(precF,rep)
    tabtempeval=zeros(prec,rep)
    reset_timer!(ISCPW.toISCPW)
    reset_timer!(ISCIN.toISCIN)
    for j1=1:ntestN
        for j2=1:tolabstestN
            for j3=1:tolreltestN
                for j4=1:rep
                    if meth=="G"
                            #valtemp=@timed ISIM_CPW(v,(ntest1[j1],gmax,mult,ALG0),(reltol0,abstol0),(tolreltest1[j3],tolabstest1[j2]))
                            #valtemp=@timed ISIM_CPW(v,(ntest1[j1],gmax,mult,ALG0),(reltol0,abstol0),(tolreltest1[j3],tolreltest1[j3]))
                            valtemp=@timed ISIM_CPW(v,(ntest1[j1],gmax,mult,ALG0),(tolreltest1[j3],tolreltest1[j3]),(tolabstest1[j2],tolabstest1[j2]))
                            tempeval=TimerOutputs.ncalls(ISCPW.toISCPW["fmult_PW"])
                            reset_timer!(ISCPW.toISCPW)
                    elseif meth=="CIN"
                            #valtemp=@timed ISIM_CIN(v,(ntest1[j1],gmax,mult,ALG0),(tolreltest1[j3],tolabstest1[j2]))
                            valtemp=@timed ISIM_CIN(v,(ntest1[j1],gmax,mult,ALG0),(tolreltest1[j3],tolreltest1[j3]))
                            tempeval=TimerOutputs.ncalls(ISCIN.toISCIN["fmult_CIN"])
                            reset_timer!(ISCIN.toISCIN)
                end
                    tabtemperr[j4]=BAS.normmax((valtemp[1])[:,end])
                    tabtempcomp[j4]=valtemp[2]
                    tabtempeval[j4]=tempeval
                end
                vecttempabs=real.(tabtemperr)+im*abs.(imag.(tabtemperr))
                tableerr=real(abs(mean(vecttempabs-valref*ones(precF,rep,1))))
                tablecomp=mean(tabtempcomp)
                tableeval=mean(tabtempeval)
                table[(j1-1)*tolreltestN*tolabstestN+(j2-1)*tolabstestN+j3,:]=[convert(prec,ntest1[j1]) convert(prec,tolabstest1[j2]) convert(prec,tolreltest1[j3]) tablecomp tableerr tableeval]
            end
        end
    end
    return(table)
end

ntest=collect(100:100:1500)
ntest=[100]
atoltest=[1e-5]
rtoltest=[1e-5]
atoltest=[1e-2; 1e-3; 1e-4; 1e-5; 1e-6; 1e-7; 1e-8; 1e-9; 1e-10; 1e-11; 1e-12]
rtoltest=[1e-1; 1e-2; 1e-3; 1e-4; 1e-5; 1e-6; 1e-7; 1e-8; 1e-9; 1e-10; 1e-11; 1e-12]
rtoltest=[1e-3; 1e-4; 1e-5; 1e-6; 1e-7; 1e-8; 1e-9; 1e-10]
rep1=2

rtoltest=[10^(-4.5); 10^(-5); 10^(-5.5); 10^(-6); 10^(-6.5); 10^(-7); 10^(-7.5); 10^(-8);  10^(-8.5); 10^(-9); 10^(-9.5); 10^(-10);  10^(-10.5);  10^(-11); 10^(-11.5)]

meth="G"
ALG0=MethodOfSteps(Tsit5())
A_T5_G=test_cont(ntest,atoltest,rtoltest,rep1,valref)

A_RK4_S=test_cont(ntest,atoltest,rtoltest,rep1,valrefA_DE)
gr()
Plots.plot(log10.(A_RK4[:,2]),log10.(A_RK4[:,3]),log10.(A_RK4[:,5]),st=:surface,xlabel="INT tol",ylabel="NumSim tol",camera=(-0,90))
Plots.plot(log10.(A_RK4_S[:,2]),log10.(A_RK4_S[:,3]),log10.(A_RK4_S[:,5]),st=:surface,xlabel="INT tol",ylabel="NumSim tol",camera=(-0,90))

meth="CIN"
ALG0=MethodOfSteps(BS3())
A_RK4_m=test_cont(ntest,atoltest,rtoltest,rep1,valref)
A_RK4_T=test_cont(ntest,atoltest,rtoltest,rep1,valref)
A_RK4_TA=test_cont(ntest,atoltest,rtoltest,rep1,valref)

A_RK4_mNS=test_cont(ntest,atoltest,rtoltest,rep1,valrefCIN_NS)
A_RK4_TNS=test_cont(ntest,atoltest,rtoltest,rep1,valrefCIN_NS)
A_RK4_TANS=test_cont(ntest,atoltest,rtoltest,rep1,valrefCIN_NS)


ALG0=MethodOfSteps(RK4())
A_RK4_S0=test_cont(ntest,atoltest,rtoltest,rep1,valref)

valref=precF(-0.6896644288103551 + 0.8555322147675607*im)


ALG0=MethodOfSteps(RK4())
A_RK4_S6=test_cont(ntest,atoltest,rtoltest,rep1,valrefA_DE)
A_RK4_S7=test_cont(ntest,atoltest,rtoltest,rep1,valrefA_DE)
A_RK4_S8=test_cont(ntest,atoltest,rtoltest,rep1,valrefA_DE)

A_RK4_NS6=test_cont(ntest,atoltest,rtoltest,rep1,valref)
A_RK4_NS7=test_cont(ntest,atoltest,rtoltest,rep1,valref)
A_RK4_NS8=test_cont(ntest,atoltest,rtoltest,rep1,valref)

A_RK4_S6=test_cont(ntest,atoltest,rtoltest,rep1,valrefA_DE)

Plots.savefig("1.pdf")

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

n=200
gmax=12
mult=12

reltol0=0.0
abstol0=1e-6

reltol0int=1e-10
abstol0int=1e-10

ALG0=MethodOfSteps(RK4())
sajt=ISIM_CPW(v,(n,gmax,mult,ALG0),(reltol0,abstol0),(reltol0int,abstol0int))

sajt=ISIM_CIN(v,(n,gmax,mult,ALG0),(reltol0,abstol0))

BAS.normmax(sajt[:,end])

SYS.T(v)-SYS.tau(v)
meth="G"
ALG0=MethodOfSteps(RK4())
atoltest=[1e-5]
rtoltest=[1e-2; 1e-3; 1e-4; 1e-5; 1e-6; 1e-7]
ntest=[100]
rep1=2

Z1_RK4_G=test_cont(ntest,atoltest,rtoltest,rep1,valrefZ1)
println(Z1_RK4_G)
meth="CIN"

Z3_RK4_C=test_cont(ntest,atoltest,rtoltest,rep1,valrefZ3)
Z1_RK4_C=test_cont(ntest,atoltest,rtoltest,rep1,valrefZ1)
