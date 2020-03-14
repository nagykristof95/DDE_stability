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
includet("ISIM_DE_sol.jl")
includet("ISIM_LMS_RK.jl")
includet("ISIM.jl")
using Main.SIDE
using Main.SILMRK
using Main.SI

omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
v=[omega kappa delta epsilon b tau0]
vA=v
#USUAL PARAMETERS
abstol0=1e-8
ALG0=MethodOfSteps(BS3())
gmax=20
n=20000
mult=20
#END OF USUAL PARAMETERS

ALG0=MethodOfSteps(Tsit5())
sajt=ISIM_DE(v,(n,gmax,mult,ALG0))

sajt=ISIM_LMS_RK(v,(n,gmax,mult,BAS.BLM4))

sajt=ISIM(v,(n,gmax,mult,MethodOfSteps(Tsit5())))
print(BAS.normmax(sajt[:,end]))
ISIM(v,(n,gmax,mult,BAS.BRK4))

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
    reset_timer!(SILMRK.toSILMRK)
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
                            valtemp=@timed ISIM_LMS_RK(v,(ntest1[j3],gmaxtest1[j2],multtest1[j1],ALG0))
                            tempeval=TimerOutputs.ncalls(SILMRK.toSILMRK["fmult_LM"])
                            reset_timer!(SILMRK.toSILMRK)
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
ntest=vcat([10,30,50,70,100,150,250],collect(500:1000:5000))
multtest=collect(20:1:20)
rep1=1

#Calculating reference values
ALG0=BAS.BRK4
gmax=20
mult=20
nlim=20000

#point A
omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
vA=[omega kappa delta epsilon b tau0]
valrefAn=BAS.normmax(ISIM(v,(n,gmax,mult,ALG0))[:,end])
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
valrefARK4n=BAS.normmax((ISIM(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefARK4=real(valrefARK4n)+im*abs(imag.(valrefARK4n))
A_RK4=test(gmaxtest,ntest,multtest,rep1,valrefA)
ALG0=BAS.BRK3
valrefARK3n=BAS.normmax((ISIM(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefARK3=real(valrefARK3n)+im*abs(imag.(valrefARK3n))
A_RK3=test(gmaxtest,ntest,multtest,rep1,valrefARK3)
ALG0=BAS.BRK2
valrefARK2n=BAS.normmax((ISIM(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefARK2=real(valrefARK2n)+im*abs(imag.(valrefARK2n))
A_RK2=test(gmaxtest,ntest,multtest,rep1,valrefARK2)
ALG0=BAS.BEE
valrefARK1n=BAS.normmax((ISIM(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefARK1=real(valrefARK1n)+im*abs(imag.(valrefARK1n))
A_RK1=test(gmaxtest,ntest,multtest,rep1,valrefARK1)

Plots.plot(A_RK4[:,3],hcat(A_RK4[:,5],A_RK3[:,5],A_RK2[:,5],A_RK1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4[:,6],A_RK3[:,6],A_RK2[:,6],A_EE[:,6]),hcat(A_RK4[:,5],A_RK3[:,5],A_RK2[:,5],A_EE[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#Linear-Multistep methods

meth="LMS"
#point A
v=vA
ALG0=BAS.BLM4
valrefALM4n=BAS.normmax((ISIM_LMS_RK(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefALM4=real(valrefALM4n)+im*abs(imag.(valrefALM4n))
A_LM4=test(gmaxtest,ntest,multtest,rep1,valrefALM4)
ALG0=BAS.BLM3
valrefALM3n=BAS.normmax((ISIM_LMS_RK(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefALM3=real(valrefALM3n)+im*abs(imag.(valrefALM3n))
A_LM3=test(gmaxtest,ntest,multtest,rep1,valrefALM3)
ALG0=BAS.BLM2
valrefALM2n=BAS.normmax((ISIM_LMS_RK(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefALM2=real(valrefALM2n)+im*abs(imag.(valrefALM2n))
A_LM2=test(gmaxtest,ntest,multtest,rep1,valrefALM2)
ALG0=BAS.BLM1
valrefALM1n=BAS.normmax((ISIM_LMS_RK(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefALM1=real(valrefALM1n)+im*abs(imag.(valrefALM1n))
A_LM1=test(gmaxtest,ntest,multtest,rep1,valrefALM1)

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_LM4[:,6],A_LM3[:,6],A_LM2[:,6],A_LM1[:,6]),hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#DifferentialEquations.jl limited step
meth="DE"

#point A
v=vA
ALG0=MethodOfSteps(Euler())
valrefARK1Jn=BAS.normmax((ISIM_DE(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefARK1J=real(valrefARK1Jn)+im*abs(imag.(valrefARK1Jn))
A_EE_J=test(gmaxtest,ntest,multtest,rep1,valrefARK1J)
ALG0=MethodOfSteps(BS3())
valrefABS3n=BAS.normmax((ISIM_DE(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefABS3=real(valrefABS3n)+im*abs(imag.(valrefABS3n))
A_BS3=test(gmaxtest,ntest,multtest,rep1,valrefABS3)
ALG0=MethodOfSteps(RK4())
valrefARK4Jn=BAS.normmax((ISIM_DE(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefARK4J=real(valrefARK4Jn)+im*abs(imag.(valrefARK4Jn))
A_RK4_J=test(gmaxtest,ntest,multtest,rep1,valrefARK4J)
ALG0=MethodOfSteps(Tsit5())
valrefAT5n=BAS.normmax((ISIM_DE(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefAT5=real(valrefAT5n)+im*abs(imag.(valrefAT5n))
A_T5=test(gmaxtest,ntest,multtest,rep1,valrefAT5)


Plots.plot(A_EE_J[:,3],hcat(A_EE_J[:,5],A_BS3[:,5],A_RK4_J[:,5],A_T5[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["EE DE" "BS3 DE" "RK4 DE" "T5 DE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_EE_J[:,6],A_BS3[:,6],A_RK4_J[:,6],A_T5[:,6]),hcat(A_EE_J[:,5],A_BS3[:,5],A_RK4_J[:,5],A_T5[:,5],A_V6[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["EE DE" "BS3 DE" "RK4 DE" "T5 DE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#DifferentialEquations.jl adaptive step
meth="DE"

#point A
v=vA
ALG0=MethodOfSteps(BS3())
valrefABS3n_ad=BAS.normmax((ISIM_DE(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefABS3_ad1=real(valrefABS3n_ad)+im*abs(imag.(valrefABS3n_ad))
valrefABS3n_ad=BAS.normmax((ISIM_DE(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefABS3_ad2=real(valrefABS3n_ad)+im*abs(imag.(valrefABS3n_ad))
valrefABS3n_ad=BAS.normmax((ISIM_DE(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefABS3_ad3=real(valrefABS3n_ad)+im*abs(imag.(valrefABS3n_ad))
valrefABS3_ad=mean([valrefABS3_ad1,valrefABS3_ad2,valrefABS3_ad3])
A_BS3_ad=test(gmaxtest,ntest,multtest,rep1,valrefABS3_ad)
ALG0=MethodOfSteps(RK4())
valrefARK4Jn_ad=BAS.normmax((ISIM_DE(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefARK4J_ad=real(valrefARK4Jn_ad)+im*abs(imag.(valrefARK4Jn_ad))
A_RK4_J_ad=test(gmaxtest,ntest,multtest,rep1,valrefARK4J_ad)
ALG0=MethodOfSteps(Tsit5())
valrefAT5n_ad=BAS.normmax((ISIM_DE(vA,(nlim,gmax,mult,ALG0)))[:,end])
valrefAT5_ad=real(valrefAT5n_ad)+im*abs(imag.(valrefAT5n_ad))
A_T5_ad=test(gmaxtest,ntest,multtest,rep1,valrefAT5_ad)


Plots.plot(A_EE_J_ad[:,3],hcat(A_EE_J_ad[:,5],A_BS3_ad[:,5],A_RK4_J_ad[:,5],A_T5_ad[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["EE DE" "BS3 DE" "RK4 DE" "T5 DE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_EE_J_ad[:,6],A_BS3_ad[:,6],A_RK4_J_ad[:,6],A_T5_ad[:,6]),hcat(A_EE_J_ad[:,5],A_BS3_ad[:,5],A_RK4_J_ad[:,5],A_T5_ad[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["EE DE" "BS3 DE" "RK4 DE" "T5 DE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")



#RK3 & LM4 & BS3 & T5 compare

Plots.plot(A_RK3[:,3],hcat(A_RK3[:,4],A_LM4[:,4],A_BS3[:,4],A_T5[:,4]),xscale=:log10, yscale=:log10,title="Computation time A",label=["RK3" "LM4" "BS3" "T5"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(A_RK3[:,3],hcat(A_RK3[:,5],A_LM4[:,5],A_BS3[:,5],A_T5[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK3" "LM4" "BS3" "T5"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK3[:,4],A_LM4[:,4],A_BS3[:,4],A_T5[:,4]),hcat(A_RK3[:,5],A_LM4[:,5],A_BS3[:,5],A_T5[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["RK3" "LM4" "BS3" "T5"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")



function tangent((a1,b1),(a2,b2))
    (log10(b2)-log10(b1))/(log10(a2)-log10(a1))
end

tangent((A_RK4[1,3],A_RK4[1,5]),(A_RK4[end,3],A_RK4[end,5]))
tangent((A_RK2[1,3],A_RK2[1,5]),(A_RK2[end,3],A_RK2[end,5]))
tangent((A_RK1[1,3],A_RK1[1,5]),(A_RK1[end,3],A_RK1[end,5]))

Plots.plot(A_RK4[:,3],hcat(A_RK4[:,5],A_RK2[:,5],A_RK1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK4" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

RK4_cont=[20.0 20.0 10.0 0.0911723 0.15013512671960227 720.0; 20.0 20.0 30.0 0.2609634 0.006400808495078936 2320.0; 20.0 20.0 50.0 0.2997562 0.0020762038148035037 3920.0; 20.0 20.0 70.0 0.4448522 0.0010280452129538655 5520.0; 20.0 20.0 100.0 0.596373199 0.000495019520860554 7920.0; 20.0 20.0 150.0 1.004914299 0.00021759303695555093 11920.0; 20.0 20.0 250.0 1.7056869 7.774851426473388e-5 19920.0; 20.0 20.0 500.0 3.446358101 1.933465486168876e-5 39920.0; 20.0 20.0 1500.0 10.348499399 2.1313702158966124e-6 119920.0; 20.0 20.0 2500.0 17.5302291 7.591613638786559e-7 199920.0; 20.0 20.0 3500.0 29.0276514 3.813374791559021e-7 279920.0; 20.0 20.0 4500.0 31.297900801 2.2589780201182932e-7 359920.0]

RK4_disc=[20.0 20.0 10.0 0.1030186 0.187733310680753 720.0; 20.0 20.0 30.0 0.2897011 0.020785371749734275 2320.0; 20.0 20.0 50.0 0.344457799 0.010816993953030419 3920.0; 20.0 20.0 70.0 0.476365001 0.007346604958992352 5520.0; 20.0 20.0 100.0 0.698001501 0.004342301130358067 7920.0; 20.0 20.0 150.0 1.070955 0.003237492550089733 11920.0; 20.0 20.0 250.0 1.754919101 0.001913613871875905 19920.0; 20.0 20.0 500.0 3.5986862 0.0008862282950386961 39920.0; 20.0 20.0 1500.0 11.1406305 0.0002821569080972849 119920.0; 20.0 20.0 2500.0 19.7989574 0.00016035999251738477 199920.0; 20.0 20.0 3500.0 25.566053199 0.000108060263682328 279920.0; 20.0 20.0 4500.0 32.6547996 7.8978602785724e-5 359920.0]
println(A_RK4)
