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

function test(gmaxtest1,ntest1,multtest1,rep,point)
    valref=point
    gmaxn=size(gmaxtest1)[1]
    ntestn=size(ntest1)[1]
    multtestn=size(multtest1)[1]
    table=zeros(Float64,multtestn*ntestn*gmaxn,5)

    tabtempcomp=zeros(Float64,rep)
    tabtemperr=zeros(ComplexF64,rep)

    for j1=1:multtestn
        global mult=multtest1[j1]
        for j2=1:gmaxn
            global gmax=gmaxtest1[j2]
            for j3=1:ntestn
                global n=ntest1[j3]
                for j4=1:rep
                    if meth=="RK"
                    valtemp=@timed ISIM(v)
                elseif meth=="LMS"
                    valtemp=@timed ISIM_LMS_RK(v)
                elseif meth=="DE"
                    valtemp=@timed ISIM_DE(v)
                end
                    tabtemperr[j4]=normmax((valtemp[1])[:,end])
                    tabtempcomp[j4]=valtemp[2]
                end
                vecttempabs=real.(tabtemperr)+im*abs.(imag.(tabtemperr))
                #tableerr=real(mean(abs.(vecttempabs-valref*ones(ComplexF64,rep,1))))
                tableerr=real(abs(mean(vecttempabs-valref*ones(ComplexF64,rep,1))))
                tablecomp=mean(tabtempcomp)
                table[(j1-1)*ntestn*gmaxn+(j2-1)*ntestn+j3,:]=[convert(Float64,multtest1[j1]) convert(Float64,gmaxtest1[j2]) convert(Float64,ntest1[j3]) tablecomp tableerr]
            end
        end
    end
    return(table)
end

#Calculating reference values
method="RK"
BR=BRK4
gmax=20
mult=20
n=20000

#point A
omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
vA=[omega kappa delta epsilon b tau0]
valrefAn=normmax((ISIM(v))[:,end])
valrefA=real(valrefAn)+im*abs(imag.(valrefAn))
valrefA_lin=-0.645047283999984+0.884500092450041*im
valrefA_quad=-0.645047338561817+0.884500265852046*im

#point B
omega=1; kappa=0.2; delta=3.246; epsilon=1.0; b=-0.81; tau0=2*pi;
vB=[omega kappa delta epsilon b tau0]
valrefBn=normmax((ISIM(v))[:,end])
valrefB=real(valrefBn)+im*abs(imag.(valrefBn))

#graph parameters
gmaxtest=collect(20:1:20)
ntest=vcat([10,30,50,70,100,150,250],collect(500:1000:8000))
ntest=collect(100:200:1500)
multtest=collect(20:1:20)
rep1=2




################################################################################
#Runge-Kutta methods
meth="RK"
method="RK"
#pointA LIN
v=vA
BR=BRK4
A_RK4_l=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)
BR=BRK3
A_RK3_l=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)
BR=BRK2
A_RK2_l=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)
BR=BEE
A_EE_l=test(gmaxtest,ntest,multtest,rep1,valrefA_lin)

Plots.plot(A_RK4_l[:,3],hcat(A_RK4_l[:,4],A_RK3_l[:,4],A_RK2_l[:,4],A_EE_l[:,4]),xscale=:log10, yscale=:log10,title="Computation time",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_RK4_l[:,3],hcat(A_RK4_l[:,5],A_RK3_l[:,5],A_RK2_l[:,5],A_EE_l[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4_l[:,4],A_RK3_l[:,4],A_RK2_l[:,4],A_EE_l[:,4]),hcat(A_RK4_l[:,5],A_RK3_l[:,5],A_RK2_l[:,5],A_EE_l[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#pointA QUAD
#Runge-Kutta methods
meth="RK"
method="RK"
#pointA LIN
v=vA
BR=BRK4
A_RK4_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)
BR=BRK3
A_RK3_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)
BR=BRK2
A_RK2_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)
BR=BEE
A_EE_q=test(gmaxtest,ntest,multtest,rep1,valrefA_quad)

Plots.plot(A_RK4_q[:,3],hcat(A_RK4_q[:,4],A_RK3_q[:,4],A_RK2_q[:,4],A_EE_q[:,4]),xscale=:log10, yscale=:log10,title="Computation time",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_RK4_q[:,3],hcat(A_RK4_q[:,5],A_RK3_q[:,5],A_RK2_q[:,5],A_EE_q[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4_q[:,4],A_RK3_q[:,4],A_RK2_q[:,4],A_EE_q[:,4]),hcat(A_RK4_q[:,5],A_RK3_q[:,5],A_RK2_q[:,5],A_EE_q[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#point B
v=vB
BR=BRK4
B_RK4=test(gmaxtest,ntest,multtest,rep1,valrefB)
BR=BRK3
B_RK3=test(gmaxtest,ntest,multtest,rep1,valrefB)
BR=BRK2
B_RK2=test(gmaxtest,ntest,multtest,rep1,valrefB)
BR=BEE
B_EE=test(gmaxtest,ntest,multtest,rep1,valrefB)


Plots.plot(B_RK4[:,3],hcat(B_RK4[:,4],B_RK3[:,4],B_RK2[:,4],B_EE[:,4]),xscale=:log10, yscale=:log10,title="Computation time",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(B_RK4[:,3],hcat(B_RK4[:,5],B_RK3[:,5],B_RK2[:,5],B_EE[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(B_RK4[:,4],B_RK3[:,4],B_RK2[:,4],B_EE[:,4]),hcat(A_RK4[:,5],B_RK3[:,5],B_RK2[:,5],B_EE[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#AB compare graph

Plots.plot(A_RK4[:,3],hcat(A_RK4[:,4],A_RK3[:,4],B_RK4[:,4],B_RK3[:,4]),xscale=:log10, yscale=:log10,title="Computation time",label=["RK4 (A)" "RK3 (A)" "RK4 (B)" "RK3(B)"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_RK4[:,3],hcat(A_RK4[:,5],A_RK3[:,5],B_RK4[:,5],B_RK3[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error",label=["RK4 (A)" "RK3 (A)" "RK4 (B)" "RK3(B)"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4[:,4],A_RK3[:,4],B_RK4[:,4],B_RK3[:,4]),hcat(A_RK4[:,5],A_RK3[:,5],B_RK4[:,5],B_RK3[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods",label=["RK4 (A)" "RK3 (A)" "RK4 (B)" "RK3(B)"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#Linear-Multistep methods

meth="LMS"
#point A
v=vA
BLM=BLM4
A_LM4=test(gmaxtest,ntest,multtest,rep1,valrefA)
BLM=BLM3
A_LM3=test(gmaxtest,ntest,multtest,rep1,valrefA)
BLM=BLM2
A_LM2=test(gmaxtest,ntest,multtest,rep1,valrefA)
BLM=BLM1
A_LM1=test(gmaxtest,ntest,multtest,rep1,valrefA)

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,4],A_LM3[:,4],A_LM2[:,4],A_LM1[:,4]),xscale=:log10, yscale=:log10,title="Computation time A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_LM4[:,4],A_LM3[:,4],A_LM2[:,4],A_LM1[:,4]),hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")


#point B
v=vB
BLM=BLM4
B_LM4=test(gmaxtest,ntest,multtest,rep1,valrefB)
BLM=BLM3
B_LM3=test(gmaxtest,ntest,multtest,rep1,valrefB)
BLM=BLM2
B_LM2=test(gmaxtest,ntest,multtest,rep1,valrefB)
BLM=BLM1
B_LM1=test(gmaxtest,ntest,multtest,rep1,valrefB)

Plots.plot(B_LM4[:,3],hcat(B_LM4[:,4],B_LM3[:,4],B_LM2[:,4],B_LM1[:,4]),xscale=:log10, yscale=:log10,title="Computation time B",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(B_LM4[:,3],hcat(B_LM4[:,5],B_LM3[:,5],B_LM2[:,5],B_LM1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error B",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(B_LM4[:,4],B_LM3[:,4],B_LM2[:,4],B_LM1[:,4]),hcat(B_LM4[:,5],B_LM3[:,5],B_LM2[:,5],B_LM1[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods B",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#AB compare graph

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,4],A_LM3[:,4],B_LM4[:,4],B_LM3[:,4]),xscale=:log10, yscale=:log10,title="Computation time AB",label=["LM4 (A)" "LM3 (A)" "LM4 (B)" "LM3(B)"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,5],A_LM3[:,5],B_LM4[:,5],B_LM3[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error AB",label=["LM4 (A)" "LM3 (A)" "LM4 (B)" "LM3(B)"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_LM4[:,4],A_LM3[:,4],B_LM4[:,4],B_LM3[:,4]),hcat(A_LM4[:,5],A_LM3[:,5],B_LM4[:,5],B_LM3[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods AB",label=["LM4 (A)" "LM3 (A)" "LM4 (B)" "LM3(B)"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#DifferentialEquations.jl limited step
meth="RK"
method="Julia"

#point A
v=vA
alg1=MethodOfSteps(Euler())
A_EE_J=test(gmaxtest,ntest,multtest,rep1,valrefA)
alg1=MethodOfSteps(BS3())
A_BS3=test(gmaxtest,ntest,multtest,rep1,valrefA)
alg1=MethodOfSteps(Tsit5())
A_T5=test(gmaxtest,ntest,multtest,rep1,valrefA)
alg1=MethodOfSteps(Vern6())
A_V6=test(gmaxtest,ntest,multtest,rep1,valrefA)

alg1=MethodOfSteps(BS5())
A_BS5=test(gmaxtest,ntest,multtest,rep1,valrefA)
alg1=MethodOfSteps(RK4())
A_RK4_J=test(gmaxtest,ntest,multtest,rep1,valrefA)

Plots.plot(A_BS3[:,3],hcat(A_EE_J[:,5],A_BS3[:,5],A_RK4_J[:,5],A_BS5[:,5],A_T5[:,5],A_V6[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["EE J" "BS3" "RK4 J" "BS5" "T5" "V6"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_EE_J[:,4],A_BS3[:,4],A_RK4_J[:,4],A_BS5[:,4],A_T5[:,4],A_V6[:,4]),hcat(A_EE_J[:,5],A_BS3[:,5],A_RK4_J[:,5],A_BS5[:,5],A_T5[:,5],A_V6[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["EE J" "BS3" "RK4 J" "BS5" "T5" "V6"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(A_BS3[:,3],hcat(A_BS3[:,4],A_T5[:,4],A_RB23[:,4],A_V6[:,4]),xscale=:log10, yscale=:log10,title="Computation time A",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_BS3[:,3],hcat(A_BS3[:,5],A_T5[:,5],A_RB23[:,5],A_V6[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_BS3[:,4],A_T5[:,4],A_RB23[:,4],A_V6[:,4]),hcat(A_BS3[:,5],A_T5[:,5],A_RB23[:,5],A_V6[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#point B
v=vB
alg1=MethodOfSteps(BS3())
B_BS3=test(gmaxtest,ntest,multtest,rep1,valrefB)
alg1=MethodOfSteps(Tsit5())
B_T5=test(gmaxtest,ntest,multtest,rep1,valrefB)
alg1=MethodOfSteps(AutoTsit5(Rosenbrock23(autodiff=false)))
B_RB23=test(gmaxtest,ntest,multtest,rep1,valrefB)
alg1=MethodOfSteps(Vern6())
B_V6=test(gmaxtest,ntest,multtest,rep1,valrefB)

Plots.plot(B_BS3[:,3],hcat(B_BS3[:,4],B_T5[:,4],B_RB23[:,4],B_V6[:,4]),xscale=:log10, yscale=:log10,title="Computation time B",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(B_BS3[:,3],hcat(B_BS3[:,5],B_T5[:,5],B_RB23[:,5],B_V6[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error B",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(B_BS3[:,4],B_T5[:,4],B_RB23[:,4],B_V6[:,4]),hcat(B_BS3[:,5],B_T5[:,5],B_RB23[:,5],B_V6[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods B",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#AB compare graph

Plots.plot(A_BS3[:,3],hcat(A_BS3[:,4],B_BS3[:,4]),xscale=:log10, yscale=:log10,title="Computation time AB",label=["BS3 (A)" "BS3 (B)"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_BS3[:,3],hcat(A_BS3[:,5],B_BS3[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error AB",label=["BS3 (A)" "BS3 (B)"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_BS3[:,4],B_BS3[:,4]),hcat(A_BS3[:,5],B_BS3[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods AB",label=["BS3 (A)" "BS3 (B)"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#DifferentialEquations.jl free step
meth="RK"
method="Julia"

#point A
v=vA
alg1=MethodOfSteps(BS3())
A_BS3free=test(gmaxtest,ntest,multtest,rep1,valrefA)
alg1=MethodOfSteps(Tsit5())
A_T5free=test(gmaxtest,ntest,multtest,rep1,valrefA)

Plots.plot(A_BS3free[:,3],hcat(A_BS3free[:,4],A_T5free[:,4]),xscale=:log10, yscale=:log10,title="Computation time A",label=["BS3 free" "T5 free"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_BS3free[:,3],hcat(A_BS3free[:,5],A_T5free[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 free" "T5 free"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_BS3free[:,4],A_T5free[:,4]),hcat(A_BS3free[:,5],A_T5free[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["BS3 free" "T5 free"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#point B
v=vB
method="Julia"
alg1=MethodOfSteps(BS3())
B_BS3free=test(gmaxtest,ntest,multtest,rep1,valrefB)
alg1=MethodOfSteps(Tsit5())
B_T5free=test(gmaxtest,ntest,multtest,rep1,valrefB)

Plots.plot(B_BS3[:,3],hcat(B_BS3[:,4],B_T5[:,4],B_RB23[:,4],B_V6[:,4]),xscale=:log10, yscale=:log10,title="Computation time B",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(B_BS3[:,3],hcat(B_BS3[:,5],B_T5[:,5],B_RB23[:,5],B_V6[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error B",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(B_BS3[:,4],B_T5[:,4],B_RB23[:,4],B_V6[:,4]),hcat(B_BS3[:,5],B_T5[:,5],B_RB23[:,5],B_V6[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods B",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#AB compare graph

Plots.plot(A_BS3[:,3],hcat(A_BS3[:,4],B_BS3[:,4]),xscale=:log10, yscale=:log10,title="Computation time AB",label=["BS3 (A)" "BS3 (B)"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_BS3[:,3],hcat(A_BS3[:,5],B_BS3[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error AB",label=["BS3 (A)" "BS3 (B)"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_BS3[:,4],B_BS3[:,4]),hcat(A_BS3[:,5],B_BS3[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods AB",label=["BS3 (A)" "BS3 (B)"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#RK3 & LM4 & BS3 & T5 compare

Plots.plot(A_RK3[:,3],hcat(A_RK3[:,4],A_LM4[:,4],A_BS3[:,4],A_T5[:,4]),xscale=:log10, yscale=:log10,title="Computation time A",label=["RK3" "LM4" "BS3" "T5"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_RK3[:,3],hcat(A_RK3[:,5],A_LM4[:,5],A_BS3[:,5],A_T5[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK3" "LM4" "BS3" "T5"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK3[:,4],A_LM4[:,4],A_BS3[:,4],A_T5[:,4]),hcat(A_RK3[:,5],A_LM4[:,5],A_BS3[:,5],A_T5[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["RK3" "LM4" "BS3" "T5"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

arrn=vcat(collect(100:100:900),collect(1000:5000:25000))
arrn1=collect(30000:10000:60000)
method="RK"
BR=BRK4

function test_lq(arrn0)
    sn=size(arrn0)[1]
    table=zeros(ComplexF64,sn,2)
    for j=1:sn
        global n=arrn0[j]
        table[j,1]=arrn0[j]
        temp=normmax((ISIM(v))[:,end])
        table[j,2]=real(temp)+im*abs(imag.(temp))
    end
    return(table)
end

inttyp0=BSpline(Linear())
lin=test_lq(arrn)
inttyp0=BSpline(Quadratic(Line(OnGrid())))
quad=test_lq(arrn)

inttyp0=BSpline(Linear())
lin1=test_lq(arrn1)
inttyp0=BSpline(Quadratic(Line(OnGrid())))
quad1=test_lq(arrn1)

open("temp.txt", "w") do io
            writedlm(io,imag(quad))
end

#DifferentialEquations.jl free step
meth="DE"

#point A
v=vA
alg1=MethodOfSteps(Euler())
A_EE_DE=test(gmaxtest,ntest,multtest,rep1,valrefA)
alg1=MethodOfSteps(BS3())
A_BS3_DE=test(gmaxtest,ntest,multtest,rep1,valrefA)
alg1=MethodOfSteps(Tsit5())
A_T5_DE=test(gmaxtest,ntest,multtest,rep1,valrefA)

Plots.plot(A_BS3free[:,3],hcat(A_BS3free[:,4],A_T5free[:,4]),xscale=:log10, yscale=:log10,title="Computation time A",label=["BS3 free" "T5 free"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_BS3free[:,3],hcat(A_BS3free[:,5],A_T5free[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 free" "T5 free"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_BS3free[:,4],A_T5free[:,4]),hcat(A_BS3free[:,5],A_T5free[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["BS3 free" "T5 free"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")
