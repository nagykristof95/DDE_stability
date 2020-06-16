module SDapp

using TimerOutputs
using Printf
using Statistics
using LinearAlgebra
using Interpolations
using SparseArrays
using Random
using Revise
includet("system_definition.jl")
includet("basic_functions.jl")
using Main.SYS
using Main.BAS

export SDM_full_app

rng = MersenneTwister(1234)
toSDapp = TimerOutput()
###################### SDM ########################
mult_add=2
gmax_add=6

function SDM_full_Z(v1,(nvar,gmaxvar,multvar))
        dt=tau(v1)/(nvar)
        nmax=floor(Int,round((T(v1)/dt)))
        m2=trunc(Int,nvar/2)
        dim2=trunc(Int,dim/2)
        dimg=dim*(m2+1)
        P=Array{Array{Float64,2}}(undef,nmax) #construction of Pi matrices
        for j=1:nmax
            @timeit toSDapp "f_SDMapp" P[j]=exp(A((j-1)*dt,v1)*dt)
        end
        R=Array{Array{Float64,2}}(undef,nmax) #construction of Ri matrices
        for j=1:nmax
            R[j]=0.5*((exp(A((j-1)*dt,v1)*dt)-Matrix{Float64}(I, dim, dim))*inv(A((j-1)*dt,v1))*B((j-1)*dt,v1))
        end
        Zw=Array{Array{Float64,2}}(undef,nmax) #construction of Zi matrices
        for j=1:nmax
            P1=P[j]
            R1=R[j]
            Zw[j]=vcat(hcat(P1,zeros(Float64,dim,(nmax-2)*dim),R1,R1),hcat(Matrix{Float64}(I,(nmax)*dim,(nmax)*dim),zeros(Float64,(nmax)*dim,dim)))
        end
        #construction of final Z matrix
        Z=Zw[end]*Zw[end-1]
        for j=1:nmax-2
            Z=Z*Zw[nmax-1-j]
        end
        #evaluating stability based on largest eigenvalue
        return(Z)
end

function iter_app(Phi,multvar,mult0,gmaxvar,Sprev)
    n1=size(Phi)[1]
    if multvar==mult0
        S=randn!(rng, zeros(ComplexF64,n1,multvar))
    else
        S=hcat(Sprev,randn!(rng, zeros(ComplexF64,n1,multvar)))
    end
    V=zeros(ComplexF64,n1,multvar)
    H=zeros(ComplexF64,multvar,multvar)
    Hval0=zeros(ComplexF64,mult0,gmaxvar)
    for g=1:gmaxvar
        V=Phi*S
        H=pinv(S)*V #pseudo-inverse calculation
        eigH=eigen(H)
        Hval0[:,g]=BAS.sorteigvals(eigH.values)[1:mult0] #eigenvalue calculation
        Vj0=V*eigH.vectors #calculating of new set of eigenvectors
        Vj=zeros(ComplexF64,n1,multvar)
       for j=1:multvar
           Vj[:,j]=normalize(Vj0[:,j])
       end
       S=Vj #creating new initial solution array
   end
       return(Hval0,S)
end

function SDM_full_app(v1,(nvar,gmaxvar,multvar))
    Z=SDM_full_Z(v1,(nvar,gmaxvar,multvar))
    loop=iter_app(Z,multvar,multvar,gmaxvar,Nothing)
    if convreq(loop[1][:,end],3,0.2)
        append=iter_app(Z,multvar+mult_add,multvar,gmax_add,loop[2])
        print(gmaxvar+gmax_add)
        return(hcat(loop[1],append[1]))
    else
        print(gmaxvar)
        return(loop[1])
    end
end

end #module
