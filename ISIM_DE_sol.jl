module SIDE

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

export ISIM_DE

#USUAL PARAMETERS
abstol0=1e-8
alg1= MethodOfSteps(RK4())
gmax=10
n=100
mult=4
#END OF USUAL PARAMETERS

rng = MersenneTwister(1234)

function bc_model_DE(du,u,h,p,t) #DE.jl problem definiton
    tau,AA,BB,mult1 = p
    h( out1, p, t-tau)
    dutemp=zeros(ComplexF64,mult1*dim)
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
end

function fsol_DE_ini(tau,IF,AA,BB,tend,mult1,dt1,v1) #DE.jl solution defintion from 0 to tend
    alg=alg1
    lags = [tau]
    p = (tau,AA,BB,mult1)
    interp=it(IF)
    global out1=zeros(ComplexF64,dim*mult1)  #define a cache variable
    h(out1,p,t)=(out1.=sub(interp,t))
    tspan = (0.0,tend)
    u0 = sub(interp,0.0)
    prob = DDEProblem(bc_model_DE,u0,h,tspan,p; constant_lags=lags)
    return(solve(prob,alg,adaptive=false,dt=dt1,progress=true))
end

function fsol_DE(tau,solprev,AA,BB,tend,mult1,dt1,v1) #DE.jl solution defintion from 0 to tend
    alg=alg1
    lags = [tau]
    p = (tau,AA,BB,mult1)
    history(t)=subs(solprev,t+T(v1))
    global out1=zeros(ComplexF64,dim*mult1)  #define a cache variable
    h(out1,p,t)=(out1.=history(t))
    tspan = (0.0,tend)
    u0 = history(0.0)
    prob = DDEProblem(bc_model_DE,u0,h,tspan,p; constant_lags=lags)
    return(solve(prob,alg,adaptive=false,dt=dt1,progress=true))
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



function weighting(solarr,Hvec00)
    n1=size(solarr)[1]
    mult1=trunc(Int,size(solarr[1])[1]/dim)
    solstruct=zeros(ComplexF64,n1*dim,mult1)
    for p=1:n1
           for s=1:mult1
               for q=1:dim
                   solstruct[1+(p-1)*dim+(q-1),s]=solarr[p][1+(s-1)*dim+(q-1)]
               end
           end
    end
    Vj0=solstruct*Hvec00 #calculating of new set of eigenvectors
    Vj=zeros(ComplexF64,n1*dim,mult1)
    for j=1:mult1
        Vj[:,j]=normalize(Vj0[:,j])
    end
    for p=1:n1
        for s=1:mult1
            for q=1:dim
                     solarr[p][1+(s-1)*dim+(q-1)]=Vj[1+dim*(p-1)+(q-1),s]
            end
        end
    end
    return(solarr)
end

function giteriation(solprev,v1,dt1)
    S01=transpose(subs(solprev,collect((T(v1)-tau(v1)):dt1:T(v1))))
    solDEret=fsol_DE(tau(v1),solprev,BAS.Ai(v1),BAS.Bi(v1),T(v1),mult,dt1,v1)
    V01=transpose(subs(solDEret,collect(T(v1)-tau(v1):dt1:T(v1))))
    IT=iter_DE(S01,V01)
    Hval=IT[1]
    Hvec=IT[2]
    solDEret.u[:]=weighting(solDEret.u[:],Hvec)
    return(solDEret,Hval)
end

function subs(solarr,t)
    solarr(t)
end


############
function ISIM_DE(v1)
        dtn=tau(v1)/(n-1) #timestep
        nmax=floor(Int,round((T(v1)/dtn)))
        kint=floor(Int,(T(v1)/tau(v1)))
        nrest=nmax-kint*(n-1)
        tvec=collect(-tau(v1):dtn:(nmax*dtn)+1e-10*dtn)
        sol00=randn!(rng, zeros(ComplexF64,n,mult*dim))
        sol=zeros(ComplexF64,nmax,mult*dim)  #empty solution matrix

        Hval0=zeros(ComplexF64,mult)

        sol0=hcat(tvec,vcat(sol00,sol))
        solDE=fsol_DE_ini(tau(v1),sol0,BAS.Ai(v1),BAS.Bi(v1),T(v1),mult,dtn,v1)
        S=sol0[1:n,2:end]
        V=transpose(subs(solDE,collect(T(v1)-tau(v1):dtn:T(v1))))
        IT=iter_DE(S,V)
        Hval=IT[1]
        Hval0=hcat(Hval0,Hval)
        Hvec=IT[2]
        solDE.u[:]=weighting(solDE.u[:],Hvec)[:]
        for g=1:gmax-1
            it0=giteriation(solDE,v1,dtn)
            solDE=it0[1]
            Hval0=hcat(Hval0,it0[2])
        end
        print(gmax)
        #eturn(solDE)
        return(Hval0[:,2:end])
end

end #module


# Juno.@profiler ISIM_DE(v)
# Juno.@enter ISIM_DE(v)
#
# @time ISIM_DE(v)
