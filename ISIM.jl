module SI

using SparseArrays
using TimerOutputs
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

export ISIM, toISIM

#USUAL PARAMETERS
#END OF USUAL PARAMETERS

toISIM = TimerOutput()
rng = MersenneTwister(1234)
gmax_add=4
mult_add=4

function bc_model(du,u,h,p,t) #DE.jl problem definiton
    tau,AA,BB,mult1 = p
    h( out1, p, t-tau)
    dutemp=sparse(zeros(ComplexF64,mult1*dim))
    @timeit toISIM "fmult_Julia" for i1=0:dim:(mult1-1)*dim
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

function fsol(tau,IF,AA,BB,tend,mult1,dt1,alg0) #DE.jl solution defintion from 0 to tend
    lags = [tau]
    p = (tau,AA,BB,mult1)
    interp=it(IF)
    global out1=zeros(ComplexF64,dim*mult1)  #define a cache variable
    h(out1,p,t)=(out1.=sub(interp,t))
    tspan = (0.0,tend)
    u0 = sub(interp,0.0)
    prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
    return(solve(prob,alg0,adaptive=false,dt=dt1,progress=true))
end



###################### ISIM #######################

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

function ISIM(v1,(nvar,gmaxvar,multvar,ALG))
        if typeof(ALG) <: MethodOfSteps
            method="Julia"
        elseif typeof(ALG) <: Tuple
            method="RK"
        end
        dt=tau(v1)/(nvar-1) #timestep
        nmax=floor(Int,round((T(v1)/dt)))
        kint=floor(Int,(T(v1)/tau(v1)))
        nrest=nmax-kint*(nvar-1)
        tvec=collect(-tau(v1):dt:(nmax*dt)+1e-10*dt)
        sol00=randn!(rng, zeros(ComplexF64,nvar,multvar*dim))
        sol=zeros(ComplexF64,nmax,multvar*dim)  #empty solution matrix
        sol0m=zeros(ComplexF64,nvar,multvar*dim)
        sol0=zeros(ComplexF64,size(tvec)[1],multvar*dim+1)
        Hval0=zeros(ComplexF64,multvar)
        precised=true
        mult0=multvar
        g=1
        while g <= gmaxvar
            sol0=hcat(tvec,vcat(sol00,sol))

            if method == "Julia"
                solarr=fsol(tau(v1),sol0[1:nvar,:],BAS.Ai(v1),BAS.Bi(v1),nmax*dt,multvar,dt,ALG)

                for tv=0:nvar-1
                    for j=1:multvar*dim
                        sol0m[tv+1,j]=solarr((nmax*dt-tau(v1))+tv*dt)[j]
                    end
                end

            elseif method == "RK"
                for k=1:kint
                interp=it(sol0[1+(k-1)*(nvar-1):nvar+(k-1)*(nvar-1),:])
                    for j=1:(nvar-1)
                        sol0[nvar+j+(k-1)*(nvar-1),2:end]=transpose(butcher(real(sol0[nvar+(j-1)+(k-1)*(nvar-1),1]),interp,sol0[nvar+(j-1)+(k-1)*(nvar-1),2:end],dt,ALG,v1,tau(v1),multvar))
                    end
                end
                if nrest>0
                    interp=it(sol0[1+(kint-1)*(nvar-1):nvar+(kint-1)*(nvar-1),:])
                    for j=1:nrest
                        sol0[nvar+j+(kint-1)*(nvar-1),2:end]=transpose(butcher(real(sol0[nvar+j+(kint-1)*(nvar-1),1]),interp,sol0[nvar+j-1+(kint-1)*(nvar-1),2:end],dt,ALG,v1,tau(v1),multvar))
                    end
                end
                sol0m=sol0[end-(nvar-1):end,2:end]

            end
            resit=iter(sol0[1:nvar,2:multvar*dim+1],sol0m)
            sol00=resit[1]
            Hval0=hcat(Hval0,BAS.sorteigvals(resit[2])[1:mult0])
            if convreq(Hval0[:,end],3,0.1) && g==gmaxvar && precised==false
                multvar=multvar+mult_add
                gmaxvar=gmaxvar+gmax_add
                sol=zeros(ComplexF64,nmax,multvar*dim)
                sol00=hcat(sol00,randn!(rng, zeros(ComplexF64,nvar,mult_add*dim)))
                precised=true
            end
            g += 1
        end
        print(gmaxvar)
        return(Hval0[:,2:end])
end

end #module


#
# SI.ISIM(v)
