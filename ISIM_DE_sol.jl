module SIDE

using TimerOutputs
using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using Random
using DifferentialEquations
using Revise
includet("system_definition.jl")
includet("basic_functions.jl")
using Main.SYS
using Main.BAS

export ISIM_DE,toSIDE

toSIDE=TimerOutput()

#USUAL PARAMETERS
abstol0=1e-8
#END OF USUAL PARAMETERS

rng = MersenneTwister(1234)

function bc_model(du,u,h,p,t) #DE.jl problem definiton
    tau,AA,BB,mult1 = p
    h( out1, p, t-tau)
    dutemp=zeros(ComplexF64,mult1*dim)
    @timeit toSIDE "fmult_DE" for i1=0:dim:(mult1-1)*dim
        for j=1:dim
            for jj=1:dim
            dutemp[j+i1] = dutemp[j+i1]+AA(t)[j,jj]*u[jj+i1]+BB(t)[j,jj]*out1[jj+i1]
            end
        end
    end
    for j=1:mult1*dim
    du[j] = dutemp[j]
    end
end

function iter_DE(S1,V1)
    n1=size(S1)[1]
    mult1=trunc(Int,size(S1)[2]/dim)
    S=zeros(ComplexF64,n1*dim,mult1)
    V=zeros(ComplexF64,n1*dim,mult1)
    for p=1:n1
           for s=1:mult1
               for q=1:dim
                   S[1+(p-1)*dim+(q-1),s]=S1[p,1+(s-1)*dim+(q-1)]
                   V[1+(p-1)*dim+(q-1),s]=V1[p,1+(s-1)*dim+(q-1)]
               end
           end
       end
       H=pinv(S)*V #pseudo-inverse calculation
       eigH=eigen(H)
       return((eigH.values,eigH.vectors))
end


function subs(solarr::ODESolution,t)
    solarr(t)
end

function ISIM_DE_first(v1,(nvar,gmaxvar,multvar,alg1var))
    dtn=tau(v1)/(nvar-1) #timestep
    nmax=floor(Int,round((T(v1)/dtn)))
    kint=floor(Int,(T(v1)/tau(v1)))
    nrest=nmax-kint*(nvar-1)
    tvec=collect(-tau(v1):dtn:(nmax*dtn)+1e-10*dtn)
    sol00=randn!(rng, zeros(ComplexF64,nvar,multvar*dim))
    sol=zeros(ComplexF64,nmax,multvar*dim)  #empty solution matrix

    sol0=hcat(tvec,vcat(sol00,sol))
    int=it(sol0)

    Hval0=zeros(ComplexF64,multvar)

    alg =  alg1var
    lags = [tau(v1)]
    p = (tau(v1),BAS.Ai(v1),BAS.Bi(v1),multvar)
    global out1=zeros(ComplexF64,dim*multvar)
    h(out1,p,t)=(out1.=sub(int,t))
    u0 =sub(int,0)
    tspan = (0.0,T(v1))
    prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
    sol1=solve(prob,alg,adaptive=false,dt=dtn,progress=true)
    S=sol0[1:nvar,2:end]
    V=transpose(subs(sol1,collect(T(v1)-tau(v1):dtn:T(v1))))
    IT=iter_DE(S,V)
    Hval=IT[1]
    Hval0=hcat(Hval0,Hval)
    Hvec=IT[2]
    #eigenvalue calculation
    return(sol1,Hvec,Hval)
end

function ISIM_DE_loop(v1,initarr,(nvar,gmaxvar,multvar,alg1var))
        dtn=tau(v1)/(nvar-1) #timestep
        sol11=initarr[1]
        Hvec=initarr[2]
        Hval0=initarr[3]
        p = (tau(v1),BAS.Ai(v1),BAS.Bi(v1),multvar)
        tspan = (0.0,T(v1))
        lags = [tau(v1)]
        alg=alg1var
    for g=1:gmaxvar-1
        global out1=zeros(ComplexF64,dim*multvar)
        h(out1,p,t)=(out1.=res2(res1(sol11(t+T(v1)))*Hvec))
        u0 = res2(res1(sol11(T(v1)))*Hvec)

        prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
        sol2=solve(prob,alg,adaptive=false,dt=dtn,progress=true)

        S=zeros(ComplexF64,nvar,multvar*dim)
        V=zeros(ComplexF64,nvar,multvar*dim)
        for j=1:nvar
            S[j,:]=res1(sol11((j-1)*dtn+(T(v1)-tau(v1))))*Hvec
            V[j,:]=res1(sol2((j-1)*dtn+(T(v1)-tau(v1))))
        end

        IT=iter_DE(S,V)
        Hval=IT[1]
        Hval0=hcat(Hval0,Hval)
        Hvec=IT[2]
        sol11=sol2
    end
    print(gmaxvar)
    return(Hval0)
end

function ISIM_DE(v1,(nvar,gmaxvar,multvar,alg1var))
    first=ISIM_DE_first(v1,(nvar,gmaxvar,multvar,alg1var))
    return(ISIM_DE_loop(v1,first,(nvar,gmaxvar,multvar,alg1var)))
end

end #module
