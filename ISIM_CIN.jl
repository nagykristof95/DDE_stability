module ISCIN

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

export ISIM_CIN,toISCIN

function internalnorm(u,t)
    nxF=norm(u.x[1])
    nx1=norm(u.x[2])
    nx2=norm(u.x[3])
    # if t>0.0010471975511965978
         return(maximum([nxF;nx1;nx2]))
    # else
    #     return(nxF)
    # end

end
internalnorm(u::Number,t) = norm(u)

toISCIN = TimerOutput()
rng = MersenneTwister(1234)
grids=50

function bc_model_int(du,u,h,p,t) #DE.jl problem definiton
    tau,AA,BB,mult1,T = p
    hist=h( out1, p, t-tau)
    histT=h( out2, p, t-T)
    dutemp1=zeros(ComplexF64,mult1*dim)
    dutemp2=zeros(ComplexF64,mult1,mult1)
    dutemp3=zeros(ComplexF64,mult1,mult1)
    @timeit toISCIN "fmult_CIN" for i1=0:dim:(mult1-1)*dim
        for j=1:dim
            for jj=1:dim
            dutemp1[j+i1] = dutemp1[j+i1] + AA(t)[j,jj]*u.x[1][jj+i1] + BB(t)[j,jj]*out1.x[1][jj+i1]
            end
        end
    end
    dutemp2=transpose(res1(out2.x[1]))*res1(out2.x[1]) #constructing SS
    dutemp3=transpose(res1(out2.x[1]))*res1(u.x[1])  #constructing SV
    du[1:end]=ArrayPartition(dutemp1,dutemp2,dutemp3)
end

function fsol_int_first(v1,(n1,gmax1,mult1,ALG0),(reltol1,abstol1))
    dt1=tau(v1)/(grids-1) #pseudo-timestep
    tvec=collect(-tau(v1):dt1:1e-14*dt1)
    gridst=size(tvec)[1]
    sol00=randn!(rng, zeros(ComplexF64,gridst,mult1*dim))

    sol0=hcat(tvec,vcat(sol00))

    lags = [tau(v1),T(v1)]
    p = (tau(v1),BAS.Ai(v1),BAS.Bi(v1),mult1,T(v1))
    interpX=it(sol0)

    global out1=ArrayPartition(zeros(ComplexF64,mult1*dim),zeros(ComplexF64,mult1,mult1),zeros(ComplexF64,mult1,mult1))
    global out2=ArrayPartition(zeros(ComplexF64,mult1*dim),zeros(ComplexF64,mult1,mult1),zeros(ComplexF64,mult1,mult1))

    h(out1,p,t)=(out1.=ArrayPartition(itpdelay(interpX,mult1,t,tau(v1)),zeros(ComplexF64,mult1,mult1),zeros(ComplexF64,mult1,mult1)))
    h(out2,p,t)=(out2.=ArrayPartition(itpdelay(interpX,mult1,t,tau(v1)),zeros(ComplexF64,mult1,mult1),zeros(ComplexF64,mult1,mult1)))

    abstolARR=ArrayPartition(ones(Float64,mult1*dim),abstol1*ones(Float64,mult1,mult1),abstol1*ones(Float64,mult1,mult1))
    reltolARR=ArrayPartition(ones(Float64,mult1*dim),reltol1*ones(Float64,mult1,mult1),reltol1*ones(Float64,mult1,mult1))

    tspan = (0.0,T(v1))
    u0 = ArrayPartition(itpdelay(interpX,mult1,0.0,tau(v1)),zeros(ComplexF64,mult1,mult1),zeros(ComplexF64,mult1,mult1))
    prob = DDEProblem(bc_model_int,u0,h,tspan,p; constant_lags=lags)
    #SOL=solve(prob,ALG0,dt=T(v1)/n1,adaptive=false,save_everystep=true,progress=true)
    SOL=solve(prob,ALG0,abstol=abstol1*1e3,reltol=reltol1*1e3,adaptive=true,save_everystep=true,progress=true,internalnorm=internalnorm)
    #SOL=solve(prob,ALG0,abstol=1e3*abstolARR,reltol=1e3*reltolARR,adaptive=true,save_everystep=true,progress=true)
    println(size(SOL.u))

    SS=SOL(T(v1)).x[2]
    SV=SOL(T(v1)).x[3]

    H=inv(SS)*SV

    eigH=eigen(H)

    Hvec=eigH.vectors
    Hval=eigH.values #eigenvalue calculation
    normmult=solnormalize(SOL,Hvec,mult1,v1)
    return(SOL,Hvec,Hval,normmult)
end

function fsol_int_loop(v1,(nvar,gmaxvar,mult1,ALG0),(reltol1,abstol1),initarr) #DE.jl solution defintion from 0 to tend
    sol1=initarr[1]
    Hvec=initarr[2]
    Hval0=initarr[3]
    normmult=sparse(initarr[4])
    lags = [tau(v1),T(v1)]
    tspan = (0.0,T(v1))
    p = (tau(v1),BAS.Ai(v1),BAS.Bi(v1),mult1,T(v1))

    h(out1,p,t)=(out1.=ArrayPartition(itpdelayDE(sol1,mult1,t,tau(v1),T(v1),Hvec,normmult),zeros(ComplexF64,mult1,mult1),zeros(ComplexF64,mult1,mult1)))
    h(out2,p,t)=(out2.=ArrayPartition(itpdelayDE(sol1,mult1,t,tau(v1),T(v1),Hvec,normmult),zeros(ComplexF64,mult1,mult1),zeros(ComplexF64,mult1,mult1)))

    abstolARR=ArrayPartition(ones(Float64,mult1*dim),abstol1*ones(Float64,mult1,mult1),abstol1*ones(Float64,mult1,mult1))
    reltolARR=ArrayPartition(ones(Float64,mult1*dim),reltol1*ones(Float64,mult1,mult1),reltol1*ones(Float64,mult1,mult1))

    for j=2:gmaxvar
  #define a cache variable
        u0 = ArrayPartition(itpdelayDE(sol1,mult1,0.0,tau(v1),T(v1),Hvec,normmult),zeros(ComplexF64,mult1,mult1),zeros(ComplexF64,mult1,mult1))
        prob = DDEProblem(bc_model_int,u0,h,tspan,p; constant_lags=lags)
        #SOLloop=solve(prob,ALG0,dt=T(v1)/nvar,adaptive=false,save_everystep=true,progress=true)
        SOLloop=solve(prob,ALG0,abstol=abstol1,reltol=reltol1,adaptive=true,save_everystep=true,progress=true,internalnorm=internalnorm)
        #SOLloop=solve(prob,ALG0,abstol=abstolARR,reltol=reltolARR,adaptive=true,save_everystep=true,progress=true)
        println(size(SOLloop.u))

        SS=SOLloop(T(v1)).x[2]
        SV=SOLloop(T(v1)).x[3]

        H=inv(SS)*SV
        eigH=eigen(H)
        Hvec=eigH.vectors
        normmult=solnormalize(SOLloop,Hvec,mult1,v1)
        Hval=eigH.values #eigenvalue calculation
        Hval0=hcat(Hval0,Hval)
        sol1=SOLloop
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
        return(res2((res1(sol1v(t+T1).x[1])*Hvec0)*sparse(normmult0)))
    end
end

function solnormalize(solarray,Hvec0,multvar,v1) #lazy normalize
    n1=size(solarray.u)[1]
    n0=Int(ceil((1-tau(v1)/T(v1))*n1+1))
    arr=zeros(ComplexF64,n1-n0+1,multvar*dim)
    normvec=zeros(Float64,multvar)
    for j=n0:n1
        arr[j-n0+1,:]=(solarray[j].x[1])
    end
    arr=res1full(arr)*Hvec0
    for jj=1:multvar
        normvec[jj]=1/norm(arr[:,jj])
    end
    return(diagm(vec(normvec)))
end

function ISIM_CIN(v1,(n1,gmax1,mult1,ALG0),(reltol1,abstol1))
            first=fsol_int_first(v1,(n1,gmax1,mult1,ALG0),(reltol1,abstol1))
            return(fsol_int_loop(v1,(n1,gmax1,mult1,ALG0),(reltol1,abstol1),first))
end



end #module
