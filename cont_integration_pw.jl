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
using QuadGK
rng = MersenneTwister(1234)

grids=20
mult=2
gmax=6
n=100

#inittest=randn!(rng, zeros(Float64,grids,mult*dim))

function ISIMcontfirst(v1)
    Hval0=zeros(ComplexF64,mult)
    #inital interpolation array
    dt=tau(v1)/(grids-1)
    tvec=collect(-tau(v1):dt:0)
    init=randn!(rng, zeros(ComplexF64,grids,mult*dim))
    #init=inittest
    int=it(hcat(tvec,init))

    alg =  MethodOfSteps(BS3())
    #alg =  MethodOfSteps(AutoTsit5(Rosenbrock23(autodiff=false)))
    #alg=MethodOfSteps(RK4())
    lags = [tau(v1)]
    p = (tau(v1),Ai(v1),Bi(v1),mult)
    global out1=zeros(ComplexF64,dim*mult)
    h(out1,p,t)=(out1.=sub(int,t))
    u0 =sub(int,0)
    tspan = (0.0,T(v1))
    prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
    sol1=solve(prob,alg,reltol=1e-4,dtmax=T(v1)/(n),save_everystep=true)

    SS=quadgk(t -> transpose(res1(sub(int,t-tau(v1))))*res1(sub(int,t-tau(v1))),0,tau(v1),rtol=1e-5)[1]
    SV=quadgk(t -> transpose(res1(sub(int,t-tau(v1))))*res1(sol1(t+(T(v1)-tau(v1)))),0,tau(v1),rtol=1e-5)[1]


    H=inv(SS)*SV
    eigH=eigen(H)
    Hvec=eigH.vectors
    Hval=eigH.values #eigenvalue calculation
        return(sol1,Hvec,Hval)
end

function ISIMcontloop(v1,initarr)
        sol11=initarr[1]
        Hvec=initarr[2]
        Hval0=initarr[3]
        p = (tau(v1),Ai(v1),Bi(v1),mult)
        tspan = (0.0,T(v1))
        lags = [tau(v1)]
        alg=MethodOfSteps(BS3())
    for j=1:gmax-1
        global out1=zeros(ComplexF64,dim*mult)
        h(out1,p,t)=(out1.=res2(res1(sol11(t+T(v1)))*Hvec))
        u0 = res2(res1(sol11(T(v1)))*Hvec)

        prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
        sol2=solve(prob,alg,reltol=1e-4,dtmax=tau(v1)/n,save_everystep=true)

        SS=quadgk(t -> transpose(res1(sol11(t+(T(v1)-tau(v1))))*Hvec)*(res1(sol11(t+(T(v1)-tau(v1))))*Hvec),0,tau(v1),rtol=1e-5)[1]
        SV=quadgk(t -> transpose(res1(sol11(t+(T(v1)-tau(v1))))*Hvec)*res1(sol2(t+(T(v1)-tau(v1)))),0,tau(v1),rtol=1e-5)[1]

        H=inv(SS)*SV

        eigH=eigen(H)

        Hval=eigH.values #eigenvalue calculation

        Hval0=hcat(Hval0,Hval)
        Hvec=eigH.vectors

        sol11=sol2
    end
    return(Hval0)
end

function ISIM_cont_pw(v1)
    if abs(v1[5])<1e-10
        global mult=2
    else
        global mult=6
    end
    first=ISIMcontfirst(v1)
    return(ISIMcontloop(v1,first))
end

ISIMcont(v)
