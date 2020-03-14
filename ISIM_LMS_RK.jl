module SILMRK

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
using TimerOutputs
using DifferentialEquations
using Revise
includet("system_definition.jl")
includet("basic_functions.jl")
using Main.SYS
using Main.BAS

export ISIM_LMS_RK, toSILMRK
toSILMRK = TimerOutput()
rng = MersenneTwister(1234)


BLM1=([-1.0,0.0],[1.0,0.0])
BLM2=([0.0,-1.0,0.0],[-1/2,3/2,0.0])
BLM3=([0.0,0.0,-1.0,0.0],[5/12,-16/12,23/12,0.0])
BLM4=([0.0,0.0,0.0,-1.0,0.0],[-9/24,37/24,-59/24,55/24,0.0])
BLM5=([0.0,0.0,0.0,0.0,-1.0,0.0],[251/720,-1274/720,2616/720,-2774/720,1901/720,0.0])


function LM_RK(soltable,jstart,dt0,BLM00,v1,mult0,fmultarr0,nvar0)
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
    @timeit toSILMRK "fmult_LM" fnext=BAS.fmult(soltable[jstart,1],soltable[jstart,2:end],soltable[jstart-(nvar0-1),2:end],v1,mult0)
    Yf=Yf+bco[s]*fnext
    return(-Yy+dt0*Yf,fnext)
end

function iter(S1,V1)
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
       Hval=eigH.values #eigenvalue calculation
       # EIGVEC=eigH.vectors
       # EIGVALS=diagm(eigH.values)
       #
       # println(maximum(abs.(H-EIGVEC*EIGVALS*inv(EIGVEC))))
       Vj0=V*eigH.vectors #calculating of new set of eigenvectors
       Vj=zeros(ComplexF64,n1*dim,mult1)
       for j=1:mult1
           Vj[:,j]=normalize(Vj0[:,j])
       end
       Sj=zeros(ComplexF64,n1,mult1*dim) #creating new initial solution array
       for p=1:n1
           for s=1:mult1
               for q=1:dim
                        Sj[p,1+(s-1)*dim+(q-1)]=Vj[1+dim*(p-1)+(q-1),s]
               end
           end
       end
       return((Sj,Hval))
end

###################### ISIM with LSM #######################

function ISIM_LMS_RK(v1,(nvar,gmaxvar,multvar,algvar))
        s=size(algvar[1])[1]-1
        #s=5
        dt=tau(v1)/(nvar-1) #timestep
        nmax=floor(Int,round((T(v1)/dt)))
        tvec=collect(-(tau(v1)):dt:(nmax*dt)+1e-10*dt)
        sol00=randn!(rng, zeros(ComplexF64,nvar,multvar*dim))
        sol=zeros(ComplexF64,nmax,multvar*dim)  #empty solution matrix
        sol0m=zeros(ComplexF64,nvar,multvar*dim)
        sol0=zeros(ComplexF64,size(tvec)[1],multvar*dim+1)
        Hval0=zeros(ComplexF64,multvar)
        ni=1
        dti=dt/ni
        for g=1:gmaxvar
            sol0=hcat(tvec,vcat(sol00,sol))

            fmemory=zeros(ComplexF64,size(tvec)[1],multvar*dim)
            int=it(sol0)
            for j=1:s-1
                initsim=zeros(ComplexF64,ni+1,multvar*dim)
                initsim[1,:]=sol0[nvar+(j-1),2:end]
                for jj=1:ni
                    initsim[jj+1,:]=butcher(real(sol0[nvar+(j-1),1])+dti*(jj-1)+dt*1e-10, int, initsim[jj,:], dti, BAS.BRK5 , v1, tau(v1), multvar)
                end
                sol0[nvar+j,2:end]=initsim[end,:]
            end


            for j=1:s-1
                fmemory[nvar+(j-1),:]=BAS.fmult(real(sol0[nvar+(j-1),1]),sol0[nvar+(j-1),2:end],sol0[j,2:end],v1,multvar)
            end


            for j=nvar+(s-1):nvar+nmax-1
                Y=LM_RK(sol0,j,dt,algvar,v1,multvar,fmemory,nvar)

                sol0[j+1,2:end]=Y[1]
                fmemory[j,:]=Y[2]
            end
            sol0m=sol0[end-(nvar-1):end,2:end]
            resit=iter(sol0[1:nvar,2:end],sol0m)
            sol00=resit[1]
            Hval0=hcat(Hval0,resit[2])
        end
        print(gmaxvar)
        return(Hval0[:,2:end])
end

end #module
