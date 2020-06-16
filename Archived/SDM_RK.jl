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

function SDM_RK4(v1,nvar)
        dt=tau(v1)/(nvar)
        h=dt
        nmax=floor(Int,round((T(v1)/dt)))+1
        P=Array{Array{Float64,2}}(undef,nmax)
        Q=Array{Array{Float64,2}}(undef,nmax)
        R=Array{Array{Float64,2}}(undef,nmax)
        for j=1:nmax
            #defining sub-matrices
            Imat=Matrix{Float64}(I,dim,dim)
            Ai=A((j-1)*dt,v1)
            Ai2=A((j-1/2)*dt,v1)
            Ai1=A((j)*dt,v1)
            Bi=B((j-1)*dt,v1)
            Bi2=B((j-1/2)*dt,v1)
            Bi1=B((j)*dt,v1)
            P[j]=Imat + (Ai*h)/6 + (2*Ai2*h)/3 + 1/6*h^2*Ai1*Ai2 + (h*Ai1)/6 + 1/6*h^2*Ai2*Ai + 1/6*h^2*Ai2*Ai2 + 1/12*h^3*Ai1*Ai2*Ai2 + 1/12*h^3*Ai2*Ai2*Ai + 1/24*h^4*Ai1*Ai2*Ai2*Ai
            Q[j]=(Bi1*h)/6 + (Bi2*h)/3 + 1/12*h^2*Ai1*Bi2 + 1/12*h^2*Ai2*Bi2 + 1/24*h^3*Ai1*Ai2*Bi2
            R[j]=(Bi*h)/6 + (Bi2*h)/3 + 1/12*h^2*Ai1*Bi2 + 1/6*h^2*Ai2*Bi + 1/12*h^2*Ai2*Bi2 + 1/24*h^3*Ai1*Ai2*Bi2 + 1/12*h^3*Ai2*Ai2*Bi + 1/24*h^4*Ai1*Ai2*Ai2*Bi
        end
        Zw=Array{Array{Float64,2}}(undef,nmax) #construction of Zi matrices
        for j=1:nmax
            P1=P[j]
            Q1=Q[j]
            R1=R[j]
            Zw[j]=vcat(hcat(P1,zeros(Float64,dim,(nmax-3)*dim),Q1,R1),hcat(Matrix{Float64}(I,(nmax-1)*dim,(nmax-1)*dim),zeros(Float64,(nmax-1)*dim,dim)))
        end
        #construction of final Z matrix
        Z=Zw[end]*Zw[end-1]
        for j=1:nmax-2
            Z=Z*Zw[nmax-1-j]
        end
        #evaluating stability based on largest eigenvalue
        return(eigs(Z))
        #return(P)
end


end #module
