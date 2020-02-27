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

BLM1=([-1.0,0.0],[1.0,0.0])
BLM2=([0.0,-1.0,0.0],[-1/2,3/2,0.0])
BLM3=([0.0,0.0,-1.0,0.0],[5/12,-16/12,23/12,0.0])
BLM4=([0.0,0.0,0.0,-1.0,0.0],[-9/24,37/24,-59/24,55/24,0.0])


function LM_init(soltable,jstart,dt0,BLM00,v1,mult0)
    (aco,bco)=BLM00
    s=size(aco)[1]-1
    fmultarr=zeros(ComplexF64,(s-1),mult0*dim)
    for j=1:(s-1)
        fmultarr[j,:]=fmult(soltable[(jstart-s)+j,1],soltable[(jstart-s)+j,2:end],soltable[(jstart-s)+j-(n-1),2:end],v1,mult0)
    end
    return(fmultarr)
end

function LM(soltable,jstart,dt0,BLM00,v1,mult0,fmultarr0)
    (aco,bco)=BLM00
    s=size(aco)[1]-1
    Yy=zeros(ComplexF64,mult0*dim)
    Yf=zeros(ComplexF64,mult0*dim)
    for j=1:s
        Yy=Yy+aco[j]*soltable[(jstart-s)+j,2:end]
    end
    for j=1:s-1
        Yf=Yf+bco[j]*fmultarr0[(jstart-s)+j,:]
    end
    fnext=fmult(soltable[jstart,1],soltable[jstart,2:end],soltable[jstart-(n-1),2:end],v1,mult0)
    Yf=Yf+bco[s]*fnext
    return(-Yy+dt0*Yf,fnext)
end

function iter_LSM(S10,V10,s0)
    n1=size(S10)[1]-(s0-1)
    mult1=trunc(Int,size(S10)[2]/dim)

    S1=zeros(ComplexF64,(n1+(s0-1))*dim,mult1)
    V1=zeros(ComplexF64,(n1+(s0-1))*dim,mult1)
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
       Vj=zeros(ComplexF64,(n1+(s0-1))*dim,mult1)
       for j=1:mult1
           Vj[:,j]=normalize(Vj0[:,j])
       end
       Sj=zeros(ComplexF64,n1+(s0-1),mult1*dim) #creating new initial solution array
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

function ISIM_LSM(v1)
        s=size(BLM[1])[1]-1
        dt=tau(v1)/(n-1) #timestep
        nmax=floor(Int,round((T(v1)/dt)))
        tvec=collect(-((s-1)*dt+tau(v1)):dt:(nmax*dt)+1e-10*dt)
        sol00=randn!(rng, zeros(ComplexF64,n+(s-1),mult*dim))
        sol=zeros(ComplexF64,nmax,mult*dim)  #empty solution matrix
        sol0m=zeros(ComplexF64,n+(s-1),mult*dim)
        sol0=zeros(ComplexF64,size(tvec)[1],mult*dim+1)
        Hval0=zeros(ComplexF64,mult)

        for g=1:gmax
            sol0=hcat(tvec,vcat(sol00,sol))

            fmemory=zeros(ComplexF64,n+nmax+(s-1),mult*dim)
            fmemory[n+1+(s-1)-(s-1):n+1+(s-1)-1,:]=LM_init(sol0,n+(s-1),dt,BLM,v1,mult)

            for j=0:nmax-1
                Y=LM(sol0,j+n+(s-1),dt,BLM,v1,mult,fmemory)

                sol0[j+n+1+(s-1),2:end]=Y[1]
                fmemory[j+n+(s-1),:]=Y[2]
            end
            sol0m=sol0[end-(n+s-2):end,2:end]
            resit=iter_LSM(sol0[1:n+(s-1),2:end],sol0m,s)
            sol00=resit[1]
            Hval0=hcat(Hval0,resit[2])
        end
        print(gmax)
        return(Hval0[:,2:end])
end

n=50
mult=4
gmax=5
BLM=BLM2

sajt=ISIM_LSM(v)
