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

function test(gmaxtest1,ntest1,multtest1,rep)
    valref=valrefA
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
                    valtemp=@timed ISIM_LSM(v)
                end
                    tabtemperr[j4]=normmax((valtemp[1])[:,end])
                    tabtempcomp[j4]=valtemp[2]
                end
                vecttempabs=real.(tabtemperr)+im*abs.(imag.(tabtemperr))
                tableerr=real(mean(abs.(vecttempabs-valref*ones(ComplexF64,rep,1))))
                tablecomp=mean(tabtempcomp)
                table[(j1-1)*ntestn*gmaxn+(j2-1)*ntestn+j3,:]=[convert(Float64,multtest1[j1]) convert(Float64,gmaxtest1[j2]) convert(Float64,ntest1[j3]) tablecomp tableerr]
            end
        end
    end
    return(table)
end

#DifferentialEquations.jl tolerance
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
ntest=collect(1000:250:5000)
multtest=collect(20:1:20)

method="Julia"
alg1=MethodOfSteps(BS3())
abstol0=1e-12
A_BS3_tol12=test(gmaxtest,ntest,multtest,2)
abstol0=1e-8
A_BS3_tol8=test(gmaxtest,ntest,multtest,2)
abstol0=1e-4
A_BS3_tol4=test(gmaxtest,ntest,multtest,2)


Plots.plot(A_BS3_tol12[:,3],hcat(A_BS3_tol12[:,4],A_BS3_tol8[:,4],A_BS3_tol4[:,4]),xscale=:log10, yscale=:log10,title="Computation time A",label=["BS3 tol12" "BS3 tol8" "BS3 tol4"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_BS3_tol12[:,3],hcat(A_BS3_tol12[:,5],A_BS3_tol8[:,5],A_BS3_tol4[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 tol12" "BS3 tol8" "BS3 tol4"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_BS3_tol12[:,4],A_BS3_tol8[:,4],A_BS3_tol4[:,4]),hcat(A_BS3_tol12[:,5],A_BS3_tol8[:,5],A_BS3_tol4[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["BS3 tol12" "BS3 tol8" "BS3 tol4"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")


#Runge-Kutta methods
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
ntest=collect(1000:250:5000)
multtest=collect(20:1:20)

BR=BRK4
A_RK4=test(gmaxtest,ntest,multtest,5)
BR=BRK3
A_RK3=test(gmaxtest,ntest,multtest,5)
BR=BRK2
A_RK2=test(gmaxtest,ntest,multtest,5)
BR=BEE
A_EE=test(gmaxtest,ntest,multtest,5)

Plots.plot(A_RK4[:,3],hcat(A_RK4[:,4],A_RK3[:,4],A_RK2[:,4],A_EE[:,4]),xscale=:log10, yscale=:log10,title="Computation time",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_RK4[:,3],hcat(A_RK4[:,5],A_RK3[:,5],A_RK2[:,5],A_EE[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4[:,4],A_RK3[:,4],A_RK2[:,4],A_EE[:,4]),hcat(A_RK4[:,5],A_RK3[:,5],A_RK2[:,5],A_EE[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#point B
omega=1; kappa=0.2; delta=3.246; epsilon=1.0; b=-0.81; tau0=2*pi;
v=[omega kappa delta epsilon b tau0]
gmax=20
mult=20
n=20000
method="RK"
BR=BRK4
valrefAn=normmax((ISIM(v))[:,end])
valrefA=real(valrefAn)+im*abs(imag.(valrefAn))

gmaxtest=collect(20:1:20)
ntest=collect(1000:250:5000)
multtest=collect(20:1:20)

BR=BRK4
B_RK4=test(gmaxtest,ntest,multtest,5)
BR=BRK3
B_RK3=test(gmaxtest,ntest,multtest,5)
BR=BRK2
B_RK2=test(gmaxtest,ntest,multtest,5)
BR=BEE
B_EE=test(gmaxtest,ntest,multtest,5)

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
ntest=collect(1000:250:5000)
multtest=collect(20:1:20)

BLM=BLM4
A_LM4=test(gmaxtest,ntest,multtest,5)
BLM=BLM3
A_LM3=test(gmaxtest,ntest,multtest,5)
BLM=BLM2
A_LM2=test(gmaxtest,ntest,multtest,5)
BLM=BLM1
A_LM1=test(gmaxtest,ntest,multtest,5)

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,4],A_LM3[:,4],A_LM2[:,4],A_LM1[:,4]),xscale=:log10, yscale=:log10,title="Computation time A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_LM4[:,4],A_LM3[:,4],A_LM2[:,4],A_LM1[:,4]),hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")


#point B
omega=1; kappa=0.2; delta=3.246; epsilon=1.0; b=-0.81; tau0=2*pi;
v=[omega kappa delta epsilon b tau0]
gmax=20
mult=20
n=20000
method="RK"
BR=BRK4
valrefAn=normmax((ISIM(v))[:,end])
valrefA=real(valrefAn)+im*abs(imag.(valrefAn))

gmaxtest=collect(20:1:20)
ntest=collect(1000:250:5000)
multtest=collect(20:1:20)

BLM=BLM4
B_LM4=test(gmaxtest,ntest,multtest,5)
BLM=BLM3
B_LM3=test(gmaxtest,ntest,multtest,5)
BLM=BLM2
B_LM2=test(gmaxtest,ntest,multtest,5)
BLM=BLM1
B_LM1=test(gmaxtest,ntest,multtest,5)

Plots.plot(B_LM4[:,3],hcat(B_LM4[:,4],B_LM3[:,4],B_LM2[:,4],B_LM1[:,4]),xscale=:log10, yscale=:log10,title="Computation time B",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(B_LM4[:,3],hcat(B_LM4[:,5],B_LM3[:,5],B_LM2[:,5],B_LM1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error B",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(B_LM4[:,4],B_LM3[:,4],B_LM2[:,4],B_LM1[:,4]),hcat(B_LM4[:,5],B_LM3[:,5],B_LM2[:,5],B_LM1[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods B",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#AB compare graph

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,4],A_LM3[:,4],B_LM4[:,4],B_LM3[:,4]),xscale=:log10, yscale=:log10,title="Computation time AB",label=["LM4 (A)" "LM3 (A)" "LM4 (B)" "LM3(B)"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,5],A_LM3[:,5],B_LM4[:,5],B_LM3[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error AB",label=["LM4 (A)" "LM3 (A)" "LM4 (B)" "LM3(B)"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_LM4[:,4],A_LM3[:,4],B_LM4[:,4],B_LM3[:,4]),hcat(A_LM4[:,5],A_LM3[:,5],B_LM4[:,5],B_LM3[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods AB",label=["LM4 (A)" "LM3 (A)" "LM4 (B)" "LM3(B)"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")



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
ntest=collect(1000:250:5000)
multtest=collect(20:1:20)

method="Julia"
alg1=MethodOfSteps(BS3())
A_BS3=test(gmaxtest,ntest,multtest,5)
alg1=MethodOfSteps(Tsit5())
A_T5=test(gmaxtest,ntest,multtest,5)
alg1=MethodOfSteps(AutoTsit5(Rosenbrock23(autodiff=false)))
A_RB23=test(gmaxtest,ntest,multtest,5)
alg1=MethodOfSteps(Vern6())
A_V6=test(gmaxtest,ntest,multtest,5)

Plots.plot(A_BS3[:,3],hcat(A_BS3[:,4],A_T5[:,4],A_RB23[:,4],A_V6[:,4]),xscale=:log10, yscale=:log10,title="Computation time A",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_BS3[:,3],hcat(A_BS3[:,5],A_T5[:,5],A_RB23[:,5],A_V6[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_BS3[:,4],A_T5[:,4],A_RB23[:,4],A_V6[:,4]),hcat(A_BS3[:,5],A_T5[:,5],A_RB23[:,5],A_V6[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")


meth="RK"
#point B
omega=1; kappa=0.2; delta=3.246; epsilon=1.0; b=-0.81; tau0=2*pi;
v=[omega kappa delta epsilon b tau0]
gmax=20
mult=20
n=20000
method="RK"
BR=BRK4
valrefAn=normmax((ISIM(v))[:,end])
valrefA=real(valrefAn)+im*abs(imag.(valrefAn))

gmaxtest=collect(20:1:20)
ntest=collect(1000:250:5000)
multtest=collect(20:1:20)

method="Julia"
alg1=MethodOfSteps(BS3())
B_BS3=test(gmaxtest,ntest,multtest,5)
alg1=MethodOfSteps(Tsit5())
B_T5=test(gmaxtest,ntest,multtest,5)
alg1=MethodOfSteps(AutoTsit5(Rosenbrock23(autodiff=false)))
B_RB23=test(gmaxtest,ntest,multtest,5)
alg1=MethodOfSteps(Vern6())
B_V6=test(gmaxtest,ntest,multtest,5)

Plots.plot(B_BS3[:,3],hcat(B_BS3[:,4],B_T5[:,4],B_RB23[:,4],B_V6[:,4]),xscale=:log10, yscale=:log10,title="Computation time B",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(B_BS3[:,3],hcat(B_BS3[:,5],B_T5[:,5],B_RB23[:,5],B_V6[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error B",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(B_BS3[:,4],B_T5[:,4],B_RB23[:,4],B_V6[:,4]),hcat(B_BS3[:,5],B_T5[:,5],B_RB23[:,5],B_V6[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods B",label=["BS3 lim" "T5 lim" "RB23 lim" "V6 lim"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#AB compare graph

Plots.plot(A_BS3[:,3],hcat(A_BS3[:,4],B_BS3[:,4]),xscale=:log10, yscale=:log10,title="Computation time AB",label=["BS3 (A)" "BS3 (B)"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_BS3[:,3],hcat(A_BS3[:,5],B_BS3[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error AB",label=["BS3 (A)" "BS3 (B)"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_BS3[:,4],B_BS3[:,4]),hcat(A_BS3[:,5],B_BS3[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods AB",label=["BS3 (A)" "BS3 (B)"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

#DifferentialEquations.jl free step
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
ntest=collect(1000:250:5000)
multtest=collect(20:1:20)

method="Julia"
alg1=MethodOfSteps(BS3())
A_BS3free=test(gmaxtest,ntest,multtest,5)
alg1=MethodOfSteps(Tsit5())
A_T5free=test(gmaxtest,ntest,multtest,5)

Plots.plot(A_BS3free[:,3],hcat(A_BS3free[:,4],A_T5free[:,4]),xscale=:log10, yscale=:log10,title="Computation time A",label=["BS3 free" "T5 free"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_BS3free[:,3],hcat(A_BS3free[:,5],A_T5free[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 free" "T5 free"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_BS3free[:,4],A_T5free[:,4]),hcat(A_BS3free[:,5],A_T5free[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["BS3 free" "T5 free"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")


meth="RK"
#point B
omega=1; kappa=0.2; delta=3.246; epsilon=1.0; b=-0.81; tau0=2*pi;
v=[omega kappa delta epsilon b tau0]
gmax=20
mult=20
n=20000
method="RK"
BR=BRK4
valrefAn=normmax((ISIM(v))[:,end])
valrefA=real(valrefAn)+im*abs(imag.(valrefAn))

gmaxtest=collect(20:1:20)
ntest=collect(1000:250:5000)
multtest=collect(20:1:20)

method="Julia"
alg1=MethodOfSteps(BS3())
B_BS3free=test(gmaxtest,ntest,multtest,5)
alg1=MethodOfSteps(Tsit5())
B_T5free=test(gmaxtest,ntest,multtest,5)

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


open("temp.txt", "w") do io
            writedlm(io,B_V6)
end
