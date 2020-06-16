module ISCI

using Printf
using Statistics
using LinearAlgebra
using SparseArrays
using Interpolations
using Plots
pyplot()
using PyPlot
pygui(true);
using DelimitedFiles
using Random
using DifferentialEquations
using TimerOutputs
using Revise
includet("system_definition.jl")
includet("basic_functions.jl")
using Main.SYS
using Main.BAS

export ISIM_CI,toISCI

grids=50
rng = MersenneTwister(1234)
toISCI = TimerOutput()

function bc_model_int(du,u,h,p,t) #DE.jl problem definiton
    tau,AA,BB,mult1,T = p
    hist=h( out1, p, t-tau)
    histT=h( out2, p, t-T)
    du=(zeros(ComplexF64,mult1*dim),zeros(ComplexF64,mult1,mult1),zeros(ComplexF64,mult1,mult1))
    @timeit toISCI "fmult_CI" for i1=0:dim:(mult1-1)*dim
        for j=1:dim
            for jj=1:dim
            du[1][j+i1] = du[1][j+i1]+AA(t)[j,jj]*u[1][jj+i1]+BB(t)[j,jj]*out1[1][jj+i1]
            end
        end
    end
    du[2]=transpose(out2[2])*out2[2] #constructing SS
    du[3]=transpose(out2[2])*res1(u[1])  #constructing SV
end

function fsol_int_first(v1,tolvar,mult1,algvar) #DE.jl solution defintion from 0 to tend
    dt1=tau(v1)/(grids-1) #pseudo-timestep

    tvec=collect(-tau(v1):dt1:1e-14*dt1)
    gridst=size(tvec)[1]
    sol00=randn!(rng, zeros(ComplexF64,gridst,mult1*dim))

    sol0=hcat(tvec,vcat(sol00))

    lags = [tau(v1),T(v1)]
    p = (tau(v1),BAS.Ai(v1),BAS.Bi(v1),mult1,T(v1))
    interp=it(sol0)
    global out1=(zeros(ComplexF64,mult1*dim),zeros(ComplexF64,mult1,mult1),zeros(ComplexF64,mult1,mult1))
    global out2=(zeros(ComplexF64,mult1*dim),zeros(ComplexF64,mult1,mult1),zeros(ComplexF64,mult1,mult1))  #define a cache variable
    h(out1,p,t)=(out1.=(itpdelay(interp,mult1,t,tau(v1)),zeros(ComplexF64,mult1,mult1),zeros(ComplexF64,mult1,mult1)))
    h(out2,p,t)=(out2.=(itpdelay(interp,mult1,t,tau(v1)),zeros(ComplexF64,mult1,mult1),zeros(ComplexF64,mult1,mult1)))
    tspan = (0.0,T(v1))
    u0 = (itpdelay(interp,mult1,0.0,tau(v1)),zeros(ComplexF64,mult1,mult1),zeros(ComplexF64,mult1,mult1))
    prob = DDEProblem(bc_model_int,u0,h,tspan,p; constant_lags=lags)
    sol1i=solve(prob,algvar,abstol=tolvar*1e2,reltol=tolvar*1e2,adaptive=true,save_everystep=true,progress=true)
    END=sol1i(T(v1))
    return(sol1i(T(v1)))
    SS=zeros(ComplexF64,mult1,mult1) #collecting SS and SV
    SV=zeros(ComplexF64,mult1,mult1)
    for j=1:mult1
        for jj=1:mult1
            SS[j,jj]=END[mult1*dim+(j-1)*mult1+jj]
            SV[j,jj]=END[mult1*dim+mult1*mult1+(j-1)*mult1+jj]
        end
    end
    H=inv(SS)*SV
    eigH=eigen(H)
    Hvec=eigH.vectors
    Hval=eigH.values #eigenvalue calculation

    normmult=solnormalize(sol1i,Hvec,mult1)
    return(sol1i,Hvec,Hval,normmult)
end

function fsol_int_loop(v1,tolvar,mult1,algvar,initarr,gmaxvar) #DE.jl solution defintion from 0 to tend
    sol1=initarr[1]
    Hvec=initarr[2]
    Hval0=initarr[3]
    normmult=initarr[4]
    lags = [tau(v1),T(v1)]
    p = (tau(v1),BAS.Ai(v1),BAS.Bi(v1),mult1,T(v1))
    for j=2:gmaxvar
        global out1=zeros(ComplexF64,dim*mult1+mult1*mult1*2)
        global out2=zeros(ComplexF64,dim*mult1+mult1*mult1*2)   #define a cache variable
        h(out1,p,t)=(out1.=vcat(res2((res1(sol1(t+T(v1))[1:mult1*dim])*Hvec)*normmult),zeros(ComplexF64,2*mult1*mult1)))
        h(out2,p,t)=(out2.=vcat(res2((res1(sol1(t+T(v1))[1:mult1*dim])*Hvec)*normmult),zeros(ComplexF64,2*mult1*mult1)))
        u0=vcat(res2((res1(sol1(T(v1))[1:mult1*dim])*Hvec)*normmult),zeros(ComplexF64,2*mult1*mult1))
        tspan = (0.0,T(v1))
        prob = DDEProblem(bc_model_int,u0,h,tspan,p; constant_lags=lags)
        sol2=solve(prob,algvar,abstol=tolvar,reltol=tolvar,adaptive=true,save_everystep=true,progress=true)
        println(size(sol2.u))
        END=sol2(T(v1))

        SS=zeros(ComplexF64,mult1,mult1)
        SV=zeros(ComplexF64,mult1,mult1)
        for j=1:mult1 #collecting SS and SV
            for jj=1:mult1
                SS[j,jj]=END[mult1*dim+(j-1)*mult1+jj]
                SV[j,jj]=END[mult1*dim+mult1*mult1+(j-1)*mult1+jj]
            end
        end
        H=inv(SS)*SV
        eigH=eigen(H)
        Hvec=eigH.vectors
        normmult=solnormalize(sol2,Hvec,mult1)
        Hval=eigH.values #eigenvalue calculation
        Hval0=hcat(Hval0,Hval)
        sol1=sol2
    end
    print(gmaxvar)
    return(Hval0)
end

#piecewise inital function defintion - first
function itpdelay(interp0,mult1,t,tau1)
    if t<-tau1
        return(zeros(ComplexF64,mult1*dim))
    else
        return(sub(interp0,t))
    end
end

#piecewise inital function defintion - loop
function itpdelayDE(sol1v,mult1,t,tau1,T1,Hvec0,normmult0)
    if t<-tau1
        return(zeros(ComplexF64,mult1*dim))
    else
        return(vcat(res2((res1(sol1v(t+T1)[1:mult1*dim])*Hvec0)*normmult0)))
    end
end

function ISIM_CI(v1,(tolvar,gmaxvar,multvar,ALG0))
            first=fsol_int_first(v1,tolvar,multvar,ALG0)
            retun(first)
            return(fsol_int_loop(v1,tolvar,multvar,ALG0,first,gmaxvar))
end

function solnormalize(solarray,Hvec0,multvar)
    n1=size(solarray.u)[1]
    arr=zeros(ComplexF64,n1,multvar*dim)
    normvec=zeros(Float64,multvar)
    for j=1:n1
        arr[j,:]=(solarray[j])[1:multvar*dim]
    end
    arr=res1full(arr)*Hvec0
    for jj=1:multvar
        normvec[jj]=1/norm(arr[1:n1,jj])
    end
    return(diagm(vec(normvec)))
end

end #module
