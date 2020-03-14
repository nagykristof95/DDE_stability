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

Ac=3.0
t_end=1.0
C=1.0
to=TimerOutput()

function ft(t,y)
    @timeit to "f" Ac*y
end

function y_ex(t,Ac0,C0)
    [C0*exp(Ac0*t)]
end

f_t(u,p,t) = ft(t,u)

dimt=1

BaRK5=[0.0 0 0 0 0 0; 0.16666666666666666667 0 0 0 0 0; -0.21627570527696895373 0.54960903861030228706 0 0 0 0; 0.08482881411262012706 0.04162653285051884260 0.37354465303686103035 0 0 0; -0.08651098424575942561 0.37955562705964599292 0.01753570971622337002 0.35608631413655672933 0 0; -0.12499755969423778621 0.72695084642093284094 -0.38363171852137430626 0.29492374551818501854 0.32008801960982756632 0]
BbRK5=[0.07892564703041163884, 0.15537176484794180580, 0.08925647030411638840, 0.51074352969588361160, -0.30537176484794180580, 0.47107435296958836116]
BcRK5=[0.0, 1/6, 1/3, 1/2, 2/3, 5/6]
BRK5=(BaRK5,BbRK5,BcRK5)

BaRK4=[0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0]
BbRK4=[1/6, 1/3, 1/3, 1/6]
BcRK4=[0, 0.5, 0.5, 1]
BRK4=(BaRK4,BbRK4,BcRK4)

BaRK3=[0 0 0; 0.5 0 0;-1.0 2.0 0;]
BbRK3=[1/6, 2/3, 1/6]
BcRK3=[0, 0.5, 1.0]
BRK3=(BaRK3,BbRK3,BcRK3)

BaRK2=[0 0; 0.5 0]
BbRK2=[0, 1.0]
BcRK2=[0, 0.5]
BRK2=(BaRK2,BbRK2,BcRK2)

BaEE=[0.0]
BbEE=[1.0]
BcEE=[0.0]
BEE=(BaEE,BbEE,BcEE)

BLM1=([-1.0,0.0],[1.0,0.0])
BLM2=([0.0,-1.0,0.0],[-1/2,3/2,0.0])
BLM3=([0.0,0.0,-1.0,0.0],[5/12,-16/12,23/12,0.0])
BLM4=([0.0,0.0,0.0,-1.0,0.0],[-9/24,37/24,-59/24,55/24,0.0])

function value(t_1,algorithm,n0)
    dt1=t_1/(n0-1)
    tvec=collect(0:dt1:t_1)
    soly=zeros(Float64,n0,dimt)
    sol=hcat(tvec,soly)
    sol[1,2:end]=[C]
    if algorithm=="RK4" || algorithm=="RK3" || algorithm=="RK2" || algorithm=="RK1"
        if algorithm=="RK4"
            B0=BRK4
        elseif algorithm=="RK3"
            B0=BRK3
        elseif algorithm=="RK2"
            B0=BRK2
        elseif algorithm=="RK1"
            B0=BEE
        end
            for j=1:n0-1
                sol[j+1,2:end]=butcher_t(sol[j,1],sol[j,2:end],dt1,B0)
            end
    elseif algorithm=="LM4" || algorithm=="LM3" || algorithm=="LM2" || algorithm=="LM1"
        if algorithm=="LM4"
            L0=BLM4
        elseif algorithm=="LM3"
            L0=BLM3
        elseif algorithm=="LM2"
            L0=BLM2
        elseif algorithm=="LM1"
            L0=BLM1
        end
        s=size(L0[1])[1]-1
        ni=1
        dti=dt1/ni
        for j=1:s-1
            initsim=zeros(Float64,ni+1,dimt)
            initsim[1,:]=sol[1,2:end]
            for jj=1:ni
                initsim[jj+1,:]=butcher_t(sol[jj,1],sol[j,2:end],dti,BRK5)
            end
            sol[1+j,2:end]=initsim[end,:]
        end
        fmemory=zeros(Float64,size(tvec)[1],dimt)
        for j=1:s-1
            fmemory[j,:]=ft(sol[j,1],sol[j,2:end])
        end
        for j=1+(s-1):n0-1
            Y=LM_RK_t(sol,j,dt1,L0,fmemory)
            sol[j+1,2:end]=Y[1]
            fmemory[j,:]=Y[2]
        end
    elseif algorithm=="V6" || algorithm=="BS3"
        if algorithm=="V6"
            alg=Vern6()
        elseif algorithm=="BS3"
            alg=BS3()
        end
        u0_t = sol[1,2:end]
        tspan_t = (0.0,t_1)
        prob_t = ODEProblem(f_t,u0_t,tspan_t)
        sol_J = solve(prob_t,alg,adaptive=false,dt=dt1,progress=true)
        sol[end,2:end]=sol_J(t_1)
    end
    return(sol[end,2:end])
end

function butcher_t(t,y,dt1,(Ba1,Bb1,Bc1)) #one step by Butcher table (explicit only!)
    s=size(Bb1)[1]
    kvec=zeros(Float64,dimt,s)
    Svec=zeros(Float64,dimt,s)
    for j=1:s
        for jj=1:s
            Svec[:,j]=Svec[:,j]+Ba1[j,jj]*kvec[:,jj]
        end
        kvec[:,j]=dt1*ft(t+Bc1[j]*dt1,y+Svec[:,j])
    end
    yn=y
    for j=1:s
        yn=yn+Bb1[j]*kvec[:,j]
    end
    return(yn)
end

function LM_RK_t(soltable,jstart,dt0,BLM00,fmultarr0)
    (aco,bco)=BLM00
    s=size(aco)[1]-1
    Yy=zeros(Float64,dimt)
    Yf=zeros(Float64,dimt)
    for j=1:s
        Yy=Yy+aco[j]*soltable[(jstart-s)+j,2:end]
    end
    for j=1:s-1
        Yf=Yf+bco[j]*fmultarr0[(jstart-s)+j,:]
    end
    fnext=ft(soltable[jstart,1],soltable[jstart,2:end])
    Yf=Yf+bco[s]*fnext
    return(-Yy+dt0*Yf,fnext)
end

arr_n0=vcat([10,30,50,70],collect(100:500:10000))

rep=5
function test_t(arr_n00,method0,rep0)
    valexact0=y_ex(t_end,Ac,C)
    ns=size(arr_n00)[1]
    table=zeros(Float64,ns,4)
    table_rep=zeros(Float64,rep0,3)
    reset_timer!(to)
        for jj=1:ns
            for j=1:rep0
            table[jj,1]=arr_n0[jj]
            n00=arr_n0[jj]
            valtemp=@timed abs((value(t_end,method0,n00)-valexact0)[1])
            tempeval=TimerOutputs.ncalls(to["f"])
            reset_timer!(to)
            table_rep[j,1]=valtemp[1]
            table_rep[j,2]=valtemp[2]
            table_rep[j,3]=tempeval
            end
            table[jj,2]=mean(table_rep[:,1])
            table[jj,3]=mean(table_rep[:,2])
            table[jj,4]=mean(table_rep[:,3])
        end
    return(table)
end

RK4_t=test_t(arr_n0,"RK4",rep)
RK3_t=test_t(arr_n0,"RK3",rep)
RK2_t=test_t(arr_n0,"RK2",rep)
RK1_t=test_t(arr_n0,"RK1",rep)
LM4_t=test_t(arr_n0,"LM4",rep)
LM3_t=test_t(arr_n0,"LM3",rep)
LM2_t=test_t(arr_n0,"LM2",rep)
LM1_t=test_t(arr_n0,"LM1",rep)
V6_t=test_t(arr_n0,"V6",rep)
BS3_t=test_t(arr_n0,"BS3",rep)


Plots.plot(RK4_t[:,1],hcat(RK4_t[:,3],RK3_t[:,3],RK2_t[:,3],RK1_t[:,3],LM4_t[:,3],LM3_t[:,3],LM2_t[:,3],LM1_t[:,3],V6_t[:,3],BS3_t[:,3]),xscale=:log10, yscale=:log10,title="Computation time ODE TEST",label=["RK4" "RK3" "RK2" "RK1" "LM4" "LM3" "LM2" "LM1" "V6" "BS3"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(RK4_t[:,1],hcat(RK4_t[:,2],RK3_t[:,2],RK2_t[:,2],RK1_t[:,2],LM4_t[:,2],LM3_t[:,2],LM2_t[:,2],LM1_t[:,2],V6_t[:,2],BS3_t[:,2]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error ODE TEST",label=["RK4" "RK3" "RK2" "RK1" "LM4" "LM3" "LM2" "LM1" "V6" "BS3"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(RK4_t[:,3],RK3_t[:,3],RK2_t[:,3],RK1_t[:,3],LM4_t[:,3],LM3_t[:,3],LM2_t[:,3],LM1_t[:,3],V6_t[:,3],BS3_t[:,3]),hcat(RK4_t[:,2],RK3_t[:,2],RK2_t[:,2],RK1_t[:,2],LM4_t[:,2],LM3_t[:,2],LM2_t[:,2],LM1_t[:,2],V6_t[:,2],BS3_t[:,2]),xscale=:log10, yscale=:log10,title="Comparaison of methods ODE TEST",label=["RK4" "RK3" "RK2" "RK1" "LM4" "LM3" "LM2" "LM1" "V6" "BS3"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")


Plots.plot(RK4_t[:,1],hcat(RK4_t[:,3],RK3_t[:,3],RK2_t[:,3],RK1_t[:,3]),xscale=:log10, yscale=:log10,title="Computation time ODE TEST",label=["RK4" "RK3" "RK2" "RK1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(RK4_t[:,1],hcat(RK4_t[:,2],RK3_t[:,2],RK2_t[:,2],RK1_t[:,2]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error ODE TEST",label=["RK4" "RK3" "RK2" "RK1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(RK4_t[:,3],RK3_t[:,3],RK2_t[:,3],RK1_t[:,3]),hcat(RK4_t[:,2],RK3_t[:,2],RK2_t[:,2],RK1_t[:,2]),xscale=:log10, yscale=:log10,title="Comparaison of methods ODE TEST",label=["RK4" "RK3" "RK2" "RK1"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

Plots.plot(LM4_t[:,1],hcat(LM4_t[:,3],LM3_t[:,3],LM2_t[:,3],LM1_t[:,3]),xscale=:log10, yscale=:log10,title="Computation time ODE TEST",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(LM4_t[:,1],hcat(LM4_t[:,2],LM3_t[:,2],LM2_t[:,2],LM1_t[:,2]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error ODE TEST",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(LM4_t[:,3],LM3_t[:,3],LM2_t[:,3],LM1_t[:,3]),hcat(LM4_t[:,2],LM3_t[:,2],LM2_t[:,2],LM1_t[:,2]),xscale=:log10, yscale=:log10,title="Comparaison of methods ODE TEST",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")

Plots.plot(BS3_t[:,1],hcat(BS3_t[:,3],V6_t[:,3]),xscale=:log10, yscale=:log10,title="Computation time ODE TEST",label=["BS3" "V6"],markershapes= [:circle],xlabel = "number of steps",ylabel = "computation time")

Plots.plot(BS3_t[:,1],hcat(BS3_t[:,2],V6_t[:,2]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error ODE TEST",label=["BS3" "V6"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(BS3_t[:,3],V6_t[:,3]),hcat(BS3_t[:,2],V6_t[:,2]),xscale=:log10, yscale=:log10,title="Comparaison of methods ODE TEST",label=["BS3" "V6"],markershapes= [:circle],xlabel = "calculation time",ylabel = "eigenvalue error")
