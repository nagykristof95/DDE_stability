module SD

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

export SDM, toSDM, SDM_RK4, SDM_sp, SDM_full

rng = MersenneTwister(1234)
toSDM = TimerOutput()
###################### SDM ########################
function SDM(v1,nvar)
        dt=tau(v1)/(nvar)
        nmax=floor(Int,round((T(v1)/dt)))
        m2=trunc(Int,nvar/2)
        dim2=trunc(Int,dim/2)
        dimg=dim*(m2+1)
        P=Array{Array{Float64,2}}(undef,nmax) #construction of Pi matrices
        for j=1:nmax
            @timeit toSDM "f_SDM" P[j]=exp(A((j-1)*dt,v1)*dt)
        end
        R=Array{Array{Float64,2}}(undef,nmax) #construction of Ri matrices
        for j=1:nmax
            R[j]=0.5*((exp(A((j-1)*dt,v1)*dt)-Matrix{Float64}(I, dim, dim))*inv(A((j-1)*dt,v1))*B((j-1)*dt,v1))
        end
        Zw=Array{Array{Float64,2}}(undef,nmax) #construction of Zi matrices
        for j=1:nmax
            P1=P[j]
            R1=R[j]
            Zw[j]=vcat(hcat(P1,zeros(dim,dim*(m2-1)),R1[:,dim2+1:dim],R1[:,dim2+1:dim]),hcat(zeros(m2*dim,dim2),Matrix{Float64}(I, m2*dim,m2*dim),zeros(m2*dim,dim2)))
        end
        #construction of final Z matrix
        Z=Zw[end]*Zw[end-1]
        for j=1:nmax-2
            Z=Z*Zw[nmax-1-j]
        end
        #evaluating stability based on largest eigenvalue
        return(eigvals(Z))
        #return(P)
end

function SDM_sp(v1,nvar)
        dt=tau(v1)/(nvar)
        nmax=floor(Int,round((T(v1)/dt)))
        m2=trunc(Int,nvar/2)
        dim2=trunc(Int,dim/2)
        dimg=dim*(m2+1)
        P=Array{Array{Float64,2}}(undef,nmax) #construction of Pi matrices
        for j=1:nmax
            @timeit toSDM "f_SDM" P[j]=exp(A((j-1)*dt,v1)*dt)
        end
        R=Array{Array{Float64,2}}(undef,nmax) #construction of Ri matrices
        for j=1:nmax
            R[j]=0.5*((exp(A((j-1)*dt,v1)*dt)-Matrix{Float64}(I, dim, dim))*inv(A((j-1)*dt,v1))*B((j-1)*dt,v1))
        end
        Zw=Array{SparseMatrixCSC{Float64,Int64}}(undef,nmax) #construction of Zi matrices
        for j=1:nmax
            P1=P[j]
            R1=R[j]
            Zw[j]=sparse(vcat(hcat(P1,zeros(dim,dim*(m2-1)),R1[:,dim2+1:dim],R1[:,dim2+1:dim]),hcat(zeros(m2*dim,dim2),Matrix{Float64}(I, m2*dim,m2*dim),zeros(m2*dim,dim2))))
        end
        #construction of final Z matrix
        Z=Zw[end]*Zw[end-1]
        for j=1:nmax-2
            Z=Z*Zw[nmax-1-j]
        end
        #evaluating stability based on largest eigenvalue
        return(eigvals(Array(Z)))
        #return(P)
end

function SDM_full(v1,(nvar,gmaxvar,multvar,ALG0))
        dt=tau(v1)/(nvar)
        nmax=floor(Int,round((T(v1)/dt)))
        P=Array{Array{Float64,2}}(undef,nmax) #construction of Pi matrices
        for j=1:nmax
            @timeit toSDM "f_SDM" P[j]=exp(A((j-1)*dt,v1)*dt)
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
        if ALG0=="Full"
            print("F")

        elseif ALG0=="ISI"
            print(gmaxvar)
            return(iter(Z,multvar,gmaxvar))
            end
        return(eigvals(Z))
end

function iter(Phi,multvar,gmaxvar)
    n1=size(Phi)[1]
    S=randn!(rng, zeros(ComplexF64,n1,multvar))
    V=zeros(ComplexF64,n1,multvar)
    H=zeros(ComplexF64,multvar,multvar)
    Hval0=zeros(ComplexF64,multvar,gmaxvar)
    for g=1:gmaxvar
        V=Phi*S
        H=pinv(S)*V #pseudo-inverse calculation
        eigH=eigen(H)
        Hval0[:,g]=eigH.values #eigenvalue calculation
        Vj0=V*eigH.vectors #calculating of new set of eigenvectors
        Vj=zeros(ComplexF64,n1,multvar)
       for j=1:multvar
           Vj[:,j]=normalize(Vj0[:,j])
       end
       S=Vj #creating new initial solution array
   end
       return(Hval0)
end

end #module
