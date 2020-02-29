module SI

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

export ISIM

#USUAL PARAMETERS
method="RK"
abstol0=1e-8
alg1= MethodOfSteps(BS3())
gmax=5
n=100
mult=4
BR=BAS.BRK4
#END OF USUAL PARAMETERS

rng = MersenneTwister(1234)

function bc_model(du,u,h,p,t) #DE.jl problem definiton
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

function fsol(tau,IF,AA,BB,tend,mult1,dt1) #DE.jl solution defintion from 0 to tend
    alg=alg1
    #alg =  MethodOfSteps(BS3())
    #alg=MethodOfSteps(Tsit5())
    #alg =  MethodOfSteps(AutoTsit5(Rosenbrock23(autodiff=false)))
    #alg=MethodOfSteps(RK4())
    #alg=MethodOfSteps(Vern6())
    lags = [tau]
    p = (tau,AA,BB,mult1)
    interp=it(IF)
    global out1=zeros(ComplexF64,dim*mult1)  #define a cache variable
    h(out1,p,t)=(out1.=sub(interp,t))
    tspan = (0.0,tend)
    u0 = sub(interp,0.0)
    prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
    x=solve(prob,alg,adaptive=false,dt=dt1,progress=true)
    print(size(x.t)[1])
    return(x)
    #return(solve(prob,alg,abstol=1e-6,reltol=1e-3,saveat=collect(0.0:dt1:tend),adaptive=true))
    #return(solve(prob,alg,abstol=abstol0,reltol=abstol0*1e3,saveat=collect(0.0:dt1:tend),adaptive=true))
    #return(solve(prob,alg,abstol=abstol0,reltol=abstol0*1e3,adaptive=true))
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

function ISIM(v1)
        dt=tau(v1)/(n-1) #timestep
        nmax=floor(Int,round((T(v1)/dt)))
        kint=floor(Int,(T(v1)/tau(v1)))
        nrest=nmax-kint*(n-1)
        tvec=collect(-tau(v1):dt:(nmax*dt)+1e-10*dt)
        sol00=randn!(rng, zeros(ComplexF64,n,mult*dim))
        sol=zeros(ComplexF64,nmax,mult*dim)  #empty solution matrix
        sol0m=zeros(ComplexF64,n,mult*dim)
        sol0=zeros(ComplexF64,size(tvec)[1],mult*dim+1)
        Hval0=zeros(ComplexF64,mult)

        for g=1:gmax
            sol0=hcat(tvec,vcat(sol00,sol))

            if method == "Julia"
                solarr=fsol(tau(v1),sol0[1:n+1,:],BAS.Ai(v1),BAS.Bi(v1),nmax*dt,mult,dt)

                for tv=0:n-1
                    for j=1:mult*dim
                        sol0m[tv+1,j]=solarr((nmax*dt-tau(v1))+tv*dt)[j]
                    end
                end

            elseif method == "RK"
                for k=1:kint
                interp=it(sol0[1+(k-1)*(n-1):n+1+(k-1)*(n-1),:])
                    for j=1:(n-1)
                        sol0[n+j+(k-1)*(n-1),2:end]=transpose(butcher(real(sol0[n+(j-1)+(k-1)*(n-1),1])+1e-14,interp,sol0[n+(j-1)+(k-1)*(n-1),2:end],dt,BR,v1,tau(v1),mult))
                    end
                end
                if nrest>0
                    interp=it(sol0[1+(kint-1)*(n-1):n+(kint-1)*(n-1),:])
                    for j=1:nrest
                        sol0[n+j+(kint-1)*(n-1),2:end]=transpose(butcher(real(sol0[n+j+(kint-1)*(n-1),1]),interp,sol0[n+j-1+(kint-1)*(n-1),2:end],dt,BR,v1,tau(v1),mult))
                    end
                end
                sol0m=sol0[end-(n-1):end,2:end]

            end
            resit=iter(sol0[1:n,2:mult*dim+1],sol0m)
            sol00=resit[1]
            Hval0=hcat(Hval0,resit[2])
        end
        print(gmax)
        return(Hval0[:,2:end])
end


end #module


#
# SI.ISIM(v)
