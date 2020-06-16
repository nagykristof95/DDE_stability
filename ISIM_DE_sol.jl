module SIDE

using TimerOutputs
using Printf
using Statistics
using LinearAlgebra
using GenericSchur
using Interpolations
using DelimitedFiles
using Random
using SparseArrays
using DifferentialEquations
using Revise
includet("system_definition.jl")
includet("basic_functions.jl")
using Main.SYS
using Main.BAS

export ISIM_DE,toSIDE

toSIDE=TimerOutput()

#USUAL PARAMETERS
abstol00=1e-7
reltol00=1e-7
grids=50
#END OF USUAL PARAMETERS

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
    dutemp=sparse(zeros(precF,mult1*dim))
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
    S=zeros(precF,n1*dim,mult1)
    V=zeros(precF,n1*dim,mult1)
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
      #  EIGVEC=eigH.vectors
      #  EIGVALS=eigH.values
      #
      # println(maximum(abs.(H-EIGVEC*diagm(EIGVALS)*inv(EIGVEC))))
       return((eigH.values,eigH.vectors))
end


function subs(solarr::ODESolution,t)
    solarr(t)
end

function ISIM_DE_first(v1,(nvar,gmaxvar,multvar,alg1var),(reltol0,abstol0))
    dtgrids=tau(v1)/(grids-1) #timestep grid
    dtn=tau(v1)/(nvar-1) #timestep
    nmax=floor(Int,round(T(v1)/dtgrids))
    tvec=precF.(collect(-tau(v1):dtgrids:nmax*dtgrids+1e-14*dtgrids))
    sol00=precF.(randn!(rng, zeros(ComplexF64,grids,multvar*dim)))
    sol=zeros(precF,nmax,multvar*dim)  #empty solution matrix

    sol0=hcat(tvec,vcat(sol00,sol))
    int=it(sol0)

    Hval0=zeros(precF,multvar)

    alg =  alg1var
    lags = [tau(v1)]
    p = (tau(v1),BAS.Ai(v1),BAS.Bi(v1),multvar)
    global out1=zeros(precF,dim*multvar)
    h(out1,p,t)=(out1.=sub(int,t))
    u0 =sub(int,0)
    tspan = (0.0,T(v1))
    prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
    #sol1=solve(prob,alg,adaptive=false,dt=dtn,save_everystep=true,progress=true)
    sol1=solve(prob,alg,abstol=abstol0*1e3,reltol=reltol0*1e3,adaptive=true,save_everystep=true,progress=true)
    S=zeros(precF,nvar,multvar*dim)
    for j=1:nvar
        S[j,:]=sub(int,-tau(v1)+(j-1)*dtn)
    end
    V=transpose(subs(sol1,collect(T(v1)-tau(v1):dtn:T(v1)+1e-12*dtn)))
    IT=iter_DE(S,V)
    Hval=IT[1]
    Hval0=hcat(Hval0,Hval)
    Hvec=IT[2]
    normmult=solnormalize(sol1,Hvec)
    return(sol1,Hvec,Hval,normmult)
end

function ISIM_DE_loop(v1,initarr,(nvar,gmaxvar,multvar,alg1var),(reltol0,abstol0))
        dtn=tau(v1)/(nvar-1) #timestep
        sol11=initarr[1]
        Hvec=initarr[2]
        Hval0=initarr[3]
        normmult=initarr[4]
        p = (tau(v1),BAS.Ai(v1),BAS.Bi(v1),multvar)
        tspan = (0.0,T(v1))
        lags = [tau(v1)]
        alg=alg1var
        S=zeros(precF,nvar,multvar*dim)
        V=zeros(precF,nvar,multvar*dim)
    for g=2:gmaxvar
        global out1=zeros(precF,dim*multvar)
        h(out1,p,t)=(out1.=res2((res1(sol11(t+T(v1)))*Hvec)*normmult))
        u0 = res2((res1(sol11(T(v1)))*Hvec)*normmult)
        prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
        #sol2=solve(prob,alg,adaptive=false,dt=dtn,save_everystep=true,progress=true)
        sol2=solve(prob,alg,abstol=abstol0,reltol=reltol0,adaptive=true,save_everystep=true,progress=true)
        #println(size(sol2.u)[1])
        for j=1:nvar
            S[j,:]=(res1(sol11((j-1)*dtn+(T(v1)-tau(v1))))*Hvec)*normmult
            V[j,:]=res1(sol2((j-1)*dtn+(T(v1)-tau(v1))))
        end
        IT=iter_DE(S,V)
        Hval=IT[1]
        Hval0=hcat(Hval0,Hval)
        Hvec=IT[2]
        normmult=solnormalize(sol2,Hvec)
        sol11=sol2
        #sol11=sol2

    end
    print(gmaxvar)
    #return(sol11)
    return(Hval0)
end

function ISIM_DE(v1,(nvar,gmaxvar,multvar,alg1var),(reltol0,abstol0))
    first=ISIM_DE_first(v1,(nvar,gmaxvar,multvar,alg1var),(reltol0,abstol0))
    return(ISIM_DE_loop(v1,first,(nvar,gmaxvar,multvar,alg1var),(reltol0,abstol0)))
end

function solnormalize(solarray,Hvec0)
    n1=size(solarray.u)[1]
    mult1=trunc(Int,size(solarray.u[1])[1]/dim)
    arr=zeros(precF,n1,mult1*dim)
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
