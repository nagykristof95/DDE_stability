module SIDEA

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

export ISIM_DE_append,toSIDEA

toSIDEA=TimerOutput()

#USUAL PARAMETERS
abstol0=1e-7
reltol0=1e-7
grids=50
#END OF USUAL PARAMETERS
mult_add=4
gmax_add=4

rng = MersenneTwister(1234)

function bc_model(du,u,h,p,t) #DE.jl problem definiton
    tau,AA,BB,mult1 = p
    h( out1, p, t-tau)
    dutemp=zeros(ComplexF64,mult1*dim)
    @timeit toSIDEA "fmult_DE_append" for i1=0:dim:(mult1-1)*dim
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

function ISIM_DE_first(v1,(nvar,gmaxvar,multvar,mult0,alg1var))
    dtgrids=tau(v1)/(grids-1) #timestep grid
    dtn=tau(v1)/(nvar-1) #timestep
    nmax=floor(Int,round((T(v1)/dtgrids)))
    tvec=collect(-tau(v1):dtgrids:(nmax*dtgrids)+1e-10*dtgrids)
    sol00=randn!(rng, zeros(ComplexF64,grids,multvar*dim))
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
    #sol1=solve(prob,alg,adaptive=false,dt=dtn,save_everystep=true,progress=true)
    sol1=solve(prob,alg,abstol=abstol0*1e2,reltol=reltol0*1e2,adaptive=true,save_everystep=true,progress=true)
    S=zeros(Complex,nvar,multvar*dim)
    for j=1:nvar
        S[j,:]=sub(int,-tau(v1)+(j-1)*dtn)
    end
    V=transpose(subs(sol1,collect(T(v1)-tau(v1):dtn:T(v1)+1e-10*dtn)))
    IT=iter_DE(S,V)
    Hval=IT[1]
    Hval0=hcat(Hval0,BAS.sorteigvals(Hval)[1:mult0])
    Hvec=IT[2]
    #eigenvalue calculation
    normmult=solnormalize(sol1,Hvec)
    return(sol1,Hvec,Hval0,normmult)
end

function ISIM_DE_loop(v1,initarr,(nvar,gmaxvar,multvar,mult0,alg1var))
        dtn=tau(v1)/(nvar-1) #timestep
        sol11=initarr[1]
        Hvec=initarr[2]
        Hval0=initarr[3]
        normmult=initarr[4]
        p = (tau(v1),BAS.Ai(v1),BAS.Bi(v1),multvar)
        tspan = (0.0,T(v1))
        lags = [tau(v1)]
        alg=alg1var
    for g=2:gmaxvar
        global out1=zeros(ComplexF64,dim*multvar)
        println(size(res2((res1(sol11(0.0+T(v1)))*Hvec)*normmult)))
        h(out1,p,t)=(out1.=res2((res1(sol11(t+T(v1)))*Hvec)*normmult))
        u0 = res2((res1(sol11(T(v1)))*Hvec)*normmult)
        #h(out1,p,t)=(out1.=sol11(t+T(v1)))
        #u0 = sol11(T(v1))
        prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
        #sol2=solve(prob,alg,adaptive=false,dt=dtn,save_everystep=true,progress=true)
        sol2=solve(prob,alg,abstol=abstol0,reltol=reltol0,adaptive=true,save_everystep=true,progress=true)
        println(size(sol2.u)[1])
        S=zeros(ComplexF64,nvar,multvar*dim)
        V=zeros(ComplexF64,nvar,multvar*dim)
        for j=1:nvar
            S[j,:]=(res1(sol11((j-1)*dtn+(T(v1)-tau(v1))))*Hvec)*normmult
            V[j,:]=res1(sol2((j-1)*dtn+(T(v1)-tau(v1))))
        end

        IT=iter_DE(S,V)
        Hval=IT[1]
        Hval0=hcat(Hval0,BAS.sorteigvals(Hval)[1:mult0])
        Hvec=IT[2]
        normmult=solnormalize(sol2,Hvec)
        sol11=sol2
        #sol11=sol2

    end
    print(gmaxvar)
    return(sol11,Hvec,Hval0,normmult)
end

function ISIM_DE_append_first(v1,initarr,(nvar,gmaxvar,mult0,mult_add0,alg1var))
        dtn=tau(v1)/(nvar-1) #timestep

        dtgrids=tau(v1)/(grids-1) #timestep grids
        ngrids=floor(Int,round((T(v1)/dtgrids)))
        tvec=collect(-tau(v1):dtgrids:1e-10*dtgrids)
        sol00=randn!(rng, zeros(ComplexF64,grids,mult_add0*dim))

        sol0=hcat(tvec,sol00)
        int=it(sol0)
        multvar=mult0+mult_add0
        sol11=initarr[1]
        Hvec=initarr[2]
        Hval0=initarr[3]
        normmult=initarr[4]
        p = (tau(v1),BAS.Ai(v1),BAS.Bi(v1),multvar)
        tspan = (0.0,T(v1))
        lags = [tau(v1)]
        alg=alg1var
        global out1=zeros(ComplexF64,dim*multvar)
        println(typeof(res2((res1(sol11(0.0+T(v1)))*Hvec)*normmult)))
        println(typeof(sub(int,0.0)))
        #println(size(res2((res1(sol11(0.0+T(v1)))*Hvec)*normmult)))
        h(out1,p,t)=(out1.=[res2((res1(sol11(t+T(v1)))*Hvec)*normmult); sub(int,t)])
        #println(size(sub(int,0.0)))
        u0 = [res2((res1(sol11(T(v1)))*Hvec)*normmult); sub(int,0.0)]
        prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
        #sol2=solve(prob,alg,adaptive=false,dt=dtn,save_everystep=true,progress=true)
        sol2=solve(prob,alg,abstol=abstol0,reltol=reltol0,adaptive=true,save_everystep=true,progress=true)
        #println(size(sol2.u)[1])

        S=zeros(ComplexF64,nvar,multvar*dim)
        V=zeros(ComplexF64,nvar,multvar*dim)
        for j=1:nvar
            S[j,:]=[(res1(sol11((j-1)*dtn+(T(v1)-tau(v1))))*Hvec)*normmult res1(sub(int,-tau(v1)+(j-1)*dtn))]
            V[j,:]=res1(sol2((j-1)*dtn+(T(v1)-tau(v1))))
        end

        IT=iter_DE(S,V)
        Hval=IT[1]
        Hval0=hcat(Hval0,BAS.sorteigvals(Hval)[1:mult0])
        Hvec=IT[2]
        normmult=solnormalize(sol2,Hvec)
        sol11=sol2
        return(sol11,Hvec,Hval0,normmult)
end

function ISIM_DE_append(v1,(nvar,gmaxvar,multvar,alg1var))
    first=ISIM_DE_first(v1,(nvar,gmaxvar,multvar,multvar,alg1var))
    loop=ISIM_DE_loop(v1,first,(nvar,gmaxvar,multvar,multvar,alg1var))
    if convreq(loop[3][:,end],3,0.1)
        append_first=ISIM_DE_append_first(v1,loop,(nvar,gmaxvar,multvar,mult_add,alg1var))
        gmaxvar=gmax_add
        append_loop=ISIM_DE_loop(v1,append_first,(nvar,gmaxvar,multvar+mult_add,multvar,alg1var))
        return(append_loop[3][:,2:end])
    else
        return(loop[3][:,2:end])
    end
end

function solnormalize(solarray,Hvec0)
    n1=size(solarray.u)[1]
    mult1=trunc(Int,size(solarray.u[1])[1]/dim)
    arr=zeros(ComplexF64,n1,mult1*dim)
    normvec=zeros(Float64,mult1)
    for j=1:n1
        arr[j,:]=solarray.u[j]
    end
    arr=res1full(arr)*Hvec0
    for jj=1:mult1
        normvec[jj]=1/norm(arr[1:n1,jj])
    end
    return(diagm(vec(normvec)))
end

end #module
