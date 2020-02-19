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

function bc_model_int(du,u,h,p,t) #DE.jl problem definiton
    tau,AA,BB,mult1,T = p
    hist = h( out1, p, t-tau)
    histT = h( out2, p, t-T)
    dutemp=zeros(ComplexF64,mult*dim)
    for i1=0:dim:(mult1-1)*dim
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
    utemp=transpose(res1(out2[1:dim*mult1]))*res1(out2[1:dim*mult1])
    for j=1:mult1
        for jj=1:mult1
            du[mult1*dim+(j-1)*mult1+jj]=utemp[j,jj]
        end
    end
    utemp=transpose(res1(out2[1:dim*mult1]))*res1(u[1:dim*mult1])
    for j=1:mult1
        for jj=1:mult1
            du[mult1*dim+mult1*mult1+(j-1)*mult1+jj]=utemp[j,jj]
        end
    end
end

function fsol_int_ini(tau,IF,AA,BB,T,tend,mult1,dt1) #DE.jl solution defintion from 0 to tend
    alg=MethodOfSteps(RK4())
    lags = [tau,T]
    p = (tau,AA,BB,mult1,T)
    interp=it(IF)
    global out1=zeros(ComplexF64,dim*mult1+mult1*mult1*2)
    global out2=zeros(ComplexF64,dim*mult1+mult1*mult1*2)   #define a cache variable
    h(out1,p,t)=(out1.=sub(interp,t))
    tspan = (0.0,tend)
    u0 = sub(interp,0.0)
    prob = DDEProblem(bc_model_int,u0,h,tspan,p; constant_lags=lags)
    sol1i=(solve(prob,alg,abstol=abstol0,reltol=1e3*abstol0,adaptive=true))

    END=sol1i(tend)

    SS=zeros(ComplexF64,mult1,mult1)
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
    for j=1:size(sol1i.t)[1]
    sol1i[1:mult1*dim,j]=res2(res1(sol1i[1:mult1*dim,j])*Hvec)
    end
    sol1i[mult1*dim+1:end,:].=0.0
        return(sol1i,Hval)
end

function fsol_int(tau,solprev,AA,BB,T,tend,mult1,dt1) #DE.jl solution defintion from 0 to tend
    alg=MethodOfSteps(RK4())
    lags = [tau,T]
    p = (tau,AA,BB,mult1,T)
    global out1=zeros(ComplexF64,dim*mult1+mult1*mult1*2)
    global out2=zeros(ComplexF64,dim*mult1+mult1*mult1*2)   #define a cache variable
    h(out1,p,t)=(out1.=solprev(t+T))
    tspan = (0.0,tend)
    u0 = solprev(T)
    prob = DDEProblem(bc_model_int,u0,h,tspan,p; constant_lags=lags)

    sol1=(solve(prob,alg,abstol=abstol0,reltol=1e3*abstol0,adaptive=true))

    END=sol1(tend)

    SS=zeros(ComplexF64,mult1,mult1)
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
    for j=1:size(sol1.t)[1]
    sol1[1:mult1*dim,j]=res2(res1(sol1[1:mult1*dim,j])*Hvec)
    end
    sol1[mult1*dim+1:end,:].=0.0
        return(sol1,Hval)
end


function ISIM_cont(v1)
        dt=tau(v1)/(grids-1) #timestep
        gridsmax=floor(Int,round((T(v1)/dt)))

        tvec=collect(-tau(v1):dt:(gridsmax*dt)+1e-10*dt)
        sol00=randn!(rng, zeros(ComplexF64,grids,mult*dim))
        sol=zeros(ComplexF64,gridsmax,mult*dim)  #empty solution matrix
        solint=zeros(ComplexF64,size(tvec)[1],mult*mult*2)

        sol00=hcat(tvec,vcat(sol00,sol),solint)
        Hval0=zeros(ComplexF64,mult)

        solp=fsol_int_ini(tau(v1),sol00,Ai(v1),Bi(v1),T(v1),T(v1),mult,dt)
        solnext=solp[1]
        Hval0=solp[2]


        for g=2:gmax

            solploop=fsol_int(tau(v1),solnext,Ai(v1),Bi(v1),T(v1),T(v1),mult,dt)

            solnext=solploop[1]
            Hval0=hcat(Hval0,solploop[2])
        end
        print(gmax)
        return(Hval0[:,2:end])
end

abstol0=1e-6
grids=80
mult=8
gmax=5

@time isimc=ISIM_cont(v)
