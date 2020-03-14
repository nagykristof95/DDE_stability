module ISCI

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
    dutemp=zeros(ComplexF64,mult1*dim)
    @timeit toISCI "fmult_CI" for i1=0:dim:(mult1-1)*dim
        for j=1:dim
            for jj=1:dim
            dutemp[j+i1] = dutemp[j+i1]+AA(t)[j,jj]*u[jj+i1]+BB(t)[j,jj]*out1[jj+i1]
            end
        end
    end
    for j=1:mult1*dim
    du[j] = dutemp[j]
    end
    utemp=zeros(ComplexF64,mult1*mult1)
    utemp=transpose(res1(out2[1:dim*mult1]))*res1(out2[1:dim*mult1]) #constructing SS
    for j=1:mult1
        for jj=1:mult1
            du[mult1*dim+(j-1)*mult1+jj]=utemp[j,jj]
        end
    end
    utemp=transpose(res1(out2[1:dim*mult1]))*res1(u[1:dim*mult1]) #constructing SV
    for j=1:mult1
        for jj=1:mult1
            du[mult1*dim+mult1*mult1+(j-1)*mult1+jj]=utemp[j,jj]
        end
    end
end

function fsol_int_first(v1,tolvar,mult1,algvar) #DE.jl solution defintion from 0 to tend
    dt1=tau(v1)/(grids-1) #pseudo-timestep
    gridsmax=floor(Int,round((T(v1)/dt1)))

    tvec=collect(-tau(v1):dt1:(gridsmax*dt1)+1e-10*dt1)
    sol00=randn!(rng, zeros(ComplexF64,grids,mult1*dim))
    sol=zeros(ComplexF64,gridsmax,mult1*dim)  #empty solution matrix
    solint=zeros(ComplexF64,size(tvec)[1],mult1*mult1*2)

    sol0=hcat(tvec,vcat(sol00,sol),solint)

    lags = [tau(v1),T(v1)]
    p = (tau(v1),BAS.Ai(v1),BAS.Bi(v1),mult1,T(v1))
    interp=it(sol0)
    global out1=zeros(ComplexF64,dim*mult1+mult1*mult1*2)
    global out2=zeros(ComplexF64,dim*mult1+mult1*mult1*2)   #define a cache variable
    h(out1,p,t)=(out1.=sub(interp,t))
    tspan = (0.0,T(v1))
    u0 = sub(interp,0.0)
    prob = DDEProblem(bc_model_int,u0,h,tspan,p; constant_lags=lags)
    sol1i=solve(prob,algvar,abstol=tolvar,reltol=tolvar,adaptive=true,save_everystep=true,progress=true)
    END=sol1i(T(v1))

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
        #h(out1,p,t)=(out1.=sol1(t+T(v1)))
        h(out1,p,t)=(out1.=vcat(res2((res1(sol1(t+T(v1))[1:mult1*dim])*Hvec)*normmult),zeros(ComplexF64,2*mult1*mult1)))
        u0=vcat(res2((res1(sol1(T(v1))[1:mult1*dim])*Hvec)*normmult),zeros(ComplexF64,2*mult1*mult1))
        #u0 = sol1(T(v1))
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
        # for j=1:size(sol2.t)[1]
        # sol2[1:mult1*dim,j]=res2(res1(sol2[1:mult1*dim,j])*Hvec)
        # end
        # sol2[mult1*dim+1:end,:].=0.0
        sol1=sol2
    end
    print(gmaxvar)
    return(Hval0)
end


function ISIM_CI(v1,(tolvar,gmaxvar,multvar,ALG0))
            first=fsol_int_first(v1,tolvar,multvar,ALG0)
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
