module ISCPW

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
using GenericSchur
using Random
using DifferentialEquations
using QuadGK
using Revise
using TimerOutputs
includet("system_definition.jl")
includet("basic_functions.jl")
using Main.SYS
using Main.BAS

export ISIM_CPW,toISCPW

grids=50
toISCPW = TimerOutput()
rng = MersenneTwister(1234)

setprecision(BigFloat,128)
setrounding(BigFloat, RoundUp)

prec=Float64
precF=ComplexF64
# prec=BigFloat
# precF=Complex{BigFloat}

function bc_model(du,u,h,p,t) #DE.jl problem definiton
    tau,AA,BB,mult1 = p
    h( out1, p, t-tau)
    dutemp=zeros(precF,mult1*dim)
    @timeit toISCPW "fmult_PW" for i1=0:dim:(mult1-1)*dim
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


function ISIMcontfirst(v1,(n,gmaxvar,multvar,ALG),(reltol0,abstol0),(reltol0int,abstol0int))
        #inital interpolation array
        dtgrids=tau(v1)/(grids-1)
        dtn=tau(v1)/(n-1)
        tvec=collect(-tau(v1):dtgrids:0+1e-25*dtgrids)
        init=precF.(randn!(rng, zeros(Float64,grids,multvar*dim)))
        int=it(hcat(tvec,init))

        Hval0=zeros(precF,multvar)

        lags = [tau(v1)]
        p = (tau(v1),BAS.Ai(v1),BAS.Bi(v1),multvar)
        global out1=zeros(precF,dim*multvar)
        h(out1,p,t)=(out1.=sub(int,t))
        u0 =sub(int,0.0)
        tspan = (0.0,T(v1))
        prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
        sol1=solve(prob,ALG,abstol=abstol0*1e3,reltol=reltol0*1e3,save_everystep=true,progress=true)
        #sol1=solve(prob,ALG,adaptive=false,dt=dtn,save_everystep=true,progress=true)
        println(size(sol1.u))
        SSX=quadgk(t -> transpose(res1(sub(int,t-tau(v1))))*res1(sub(int,t-tau(v1))),0.0,tau(v1),rtol=reltol0int,atol=abstol0int,maxevals=1e7)
        SS=SSX[1]
        println(SSX[2])
        SVX=quadgk(t -> transpose(res1(sub(int,t-tau(v1))))*res1(sol1(t+(T(v1)-tau(v1)))),0.0,tau(v1),rtol=reltol0int,atol=abstol0int,maxevals=1e7)
        SV=SVX[1]
        println(SVX[2])

        H=inv(SS)*SV
        eigH=eigen(H)
        Hvec=eigH.vectors
        Hval=eigH.values #eigenvalue calculation
        normmult=solnormalize(sol1,Hvec)

    return(sol1,Hvec,Hval,normmult)
end

function ISIMcontloop(v1,initarr,(n,gmaxvar,multvar,ALG),(reltol0,abstol0),(reltol0int,abstol0int))
        dtn=tau(v1)/(n-1)
        sol11=initarr[1]
        Hvec=initarr[2]
        Hval0=initarr[3]
        normmult=initarr[4]
        p = (tau(v1),BAS.Ai(v1),BAS.Bi(v1),multvar)
        tspan = (0.0,T(v1))
        lags = [tau(v1)]
            for j=2:gmaxvar
                global out1=zeros(precF,dim*multvar)
                h(out1,p,t)=(out1.=res2((res1(sol11(t+T(v1)))*Hvec)*normmult))
                u0 = res2((res1(sol11(T(v1)))*Hvec)*normmult)

                prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
                sol2=solve(prob,ALG,abstol=abstol0,reltol=reltol0,save_everystep=true,progress=true)
                #sol2=solve(prob,ALG,adaptive=false,dt=dtn,save_everystep=true,progress=true)
                println(size(sol2.u))
                SSX=quadgk(t -> transpose((res1(sol11(t+(T(v1)-tau(v1))))*Hvec)*normmult)*(res1(sol11(t+(T(v1)-tau(v1))))*Hvec)*normmult,0.0,tau(v1),rtol=reltol0int,atol=abstol0int,maxevals=1e7)
                SS=SSX[1]
                println(SSX[2])
                SVX=quadgk(t -> transpose((res1(sol11(t+(T(v1)-tau(v1))))*Hvec)*normmult)*res1(sol2(t+(T(v1)-tau(v1)))),0.0,tau(v1),rtol=reltol0int,atol=abstol0int,maxevals=1e7)
                SV=SVX[1]
                println(SVX[2])

                H=inv(SS)*SV

                eigH=eigen(H)

                Hval=eigH.values #eigenvalue calculation
                Hval0=hcat(Hval0,Hval)
                Hvec=eigH.vectors
                normmult=solnormalize(sol2,Hvec)
                sol11=sol2
        end
    print(gmaxvar)
    return(Hval0)
end

function ISIM_CPW(v1,(nvar,gmaxvar,multvar,ALG),(reltol0,abstol0),(reltol0int,abstol0int))
    if abs(v1[5])<1e-12
        multvar=2
    else
        multvar=multvar
    end
    first=ISIMcontfirst(v1,(nvar,gmaxvar,multvar,ALG),(reltol0,abstol0),(reltol0int,abstol0int))

    return(ISIMcontloop(v1,first,(nvar,gmaxvar,multvar,ALG),(reltol0,abstol0),(reltol0int,abstol0int)))
end

function solnormalize(solarray,Hvec0)
    n1=size(solarray.u)[1]
    #n1half=trunc(Int,n1/2)
    mult1=trunc(Int,size(solarray.u[1])[1]/dim)
    arr=zeros(precF,n1,mult1*dim)
    #narr=zeros(precF,n1*dim,mult1)
    normvec=zeros(prec,mult1)
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
