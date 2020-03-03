module SILMS

using TimerOutputs
using Printf
using Statistics
using LinearAlgebra
using GenericSchur
using Interpolations
using DelimitedFiles
using Random
using DifferentialEquations
using Revise
includet("system_definition.jl")
includet("basic_functions.jl")
using Main.SYS
using Main.BAS

export ISIM_LMS, toSILMS

rng = MersenneTwister(1234)
setprecision(BigFloat,128)
setrounding(BigFloat, RoundUp)

toSILMS = TimerOutput()

function LM_init(soltable,jstart,dt0,BLM00,v1,mult0,n0)
    (aco,bco)=BLM00
    s=size(aco)[1]-1
    fmultarr=zeros(Complex{BigFloat},(s-1),mult0*dim)
    for j=1:(s-1)
    @timeit toSILMS "fmult_LMS" fmultarr[j,:]=BAS.fmult(soltable[(jstart-s)+j,1],soltable[(jstart-s)+j,2:end],soltable[(jstart-s)+j-(n0-1),2:end],v1,mult0)
    end
    return(fmultarr)
end

function LM(soltable,jstart,dt0,BLM00,v1,mult0,fmultarr0,n0)
    (aco,bco)=BLM00
    s=size(aco)[1]-1
    Yy=zeros(Complex{BigFloat},mult0*dim)
    Yf=zeros(Complex{BigFloat},mult0*dim)
    for j=1:s
        Yy=Yy+aco[j]*soltable[(jstart-s)+j,2:end]
    end
    for j=1:s-1
        Yf=Yf+bco[j]*fmultarr0[(jstart-s)+j,:]
    end
    @timeit toSILMS "fmult_LMS" fnext=BAS.fmult(soltable[jstart,1],soltable[jstart,2:end],soltable[jstart-(n0-1),2:end],v1,mult0)
    Yf=Yf+bco[s]*fnext
    return(-Yy+dt0*Yf,fnext)
end

function iter_LSM(S10,V10,s0)
    n1=size(S10)[1]-(s0-1)
    mult1=trunc(Int,size(S10)[2]/dim)

    S1=zeros(Complex{BigFloat},(n1+(s0-1))*dim,mult1)
    V1=zeros(Complex{BigFloat},(n1+(s0-1))*dim,mult1)
    for p=1:(n1+(s0-1))
           for s=1:mult1
               for q=1:dim
                   S1[1+(p-1)*dim+(q-1),s]=S10[p,1+(s-1)*dim+(q-1)]
                   V1[1+(p-1)*dim+(q-1),s]=V10[p,1+(s-1)*dim+(q-1)]
               end
           end
       end

       S=S1[dim*(s0-1)+1:end,:]
       V=V1[dim*(s0-1)+1:end,:]
       H=pinv(S)*V #pseudo-inverse calculation
       eigH=eigen(H)
       Hval=eigH.values #eigenvalue calculation
       Vj0=V1*eigH.vectors #calculating of new set of eigenvectors
       Vj=zeros(Complex{BigFloat},(n1+(s0-1))*dim,mult1)
       for j=1:mult1
           Vj[:,j]=normalize(Vj0[:,j])
       end
       Sj=zeros(Complex{BigFloat},n1+(s0-1),mult1*dim) #creating new initial solution array
       for p=1:n1+(s0-1)
           for s=1:mult1
               for q=1:dim
                        Sj[p,1+(s-1)*dim+(q-1)]=Vj[1+dim*(p-1)+(q-1),s]
               end
           end
       end
       return((Sj,Hval))
end





###################### ISIM with LSM #######################

function ISIM_LMS(v1,(nvar,gmaxvar,multvar,BLM))
        s=size(BLM[1])[1]-1
        dt=tau(v1)/(nvar-1) #timestep
        nmax=floor(Int,round((T(v1)/dt)))
        tvec=collect(-((s-1)*dt+tau(v1)):dt:(nmax*dt)+1e-10*dt)
        sol00=Complex{BigFloat}.(randn!(rng, zeros(ComplexF64,nvar+(s-1),multvar*dim)))
        sol=zeros(Complex{BigFloat},nmax,multvar*dim)  #empty solution matrix
        sol0m=zeros(Complex{BigFloat},nvar+(s-1),multvar*dim)
        sol0=zeros(Complex{BigFloat},size(tvec)[1],multvar*dim+1)
        Hval0=zeros(Complex{BigFloat},multvar)

        for g=1:gmaxvar
            sol0=hcat(tvec,vcat(sol00,sol))

            fmemory=zeros(Complex{BigFloat},nvar+nmax+(s-1),multvar*dim)
            fmemory[nvar+1+(s-1)-(s-1):nvar+1+(s-1)-1,:]=LM_init(sol0,nvar+(s-1),dt,BLM,v1,multvar,nvar)

            for j=0:nmax-1
                Y=LM(sol0,j+nvar+(s-1),dt,BLM,v1,multvar,fmemory,nvar)

                sol0[j+nvar+1+(s-1),2:end]=Y[1]
                fmemory[j+nvar+(s-1),:]=Y[2]
            end
            sol0m=sol0[end-(nvar+s-2):end,2:end]
            resit=iter_LSM(sol0[1:nvar+(s-1),2:end],sol0m,s)
            sol00=resit[1]
            Hval0=hcat(Hval0,resit[2])
        end
        print(gmaxvar)
        return(Hval0[:,2:end])
end

end #module end
