using Printf
using Statistics
using LinearAlgebra
using MDBM
using Plots
pyplot()
using PyPlot
pygui(true);
using Interpolations
using DifferentialEquations
using Random

rng = MersenneTwister(1234)

#creating iteration array
function it(A)
    inttyp=BSpline(Quadratic(Line(OnGrid())))
    #inttyp=BSpline(Constant())
    matrdim=size(A,2)-1
    step=abs(real(A[end,1]-A[1,1]))/(size(A,1)-1)
    scaleitp=real(A[1,1]):step:real(A[end,1])
    ARe=real(A)
    AIm=imag(A)
    imRe=[scale(interpolate(ARe[:,2],inttyp),scaleitp)]
    for j = 3:matrdim+1
        imRetemp=[scale(interpolate(ARe[:,j],inttyp),scaleitp)]
        imRe=vcat(imRe,imRetemp)
    end
    imIm=[scale(interpolate(AIm[:,2],inttyp),scaleitp)]
    for j = 3:matrdim+1
        imImtemp=[scale(interpolate(AIm[:,j],inttyp),scaleitp)]
        imIm=vcat(imIm,imImtemp)
    end
    return(hcat(imRe,imIm))
end

#substitution in interpolating function
function sub(it,t)
    subdim=size(it,1)
    out=zeros(ComplexF64,subdim)
    for j = 1:subdim
        out[j]=it[j,1](t)+it[j,2](t)*im
    end
    return(out)
end

dim=2

function A(t,vv)
    Amatr=zeros(Float64,dim,dim)
    #vv[1]=kappa #vv[2]=delta #vv[3]=epsilon #vv[4]=omega
    Amatr=[-vv[1] -(vv[2]+vv[3]*cos(vv[4]*t)); 1 0]
    return(Amatr)
end

function B(t,v)
    Bmatr=zeros(dim*dim)
    Bmatr=[0 v[5]; 0 0]
    return(Bmatr)
end

function bc_model(du,u,h,p,t)
    kappa,delta,epsilon,omega,b,tau,AA,BB = p
    hist = h(p, t-tau)
    dutemp=zeros(ComplexF64,dim)
    for j=1:dim
        for jj=1:dim
        dutemp[j] = dutemp[j]+AA(t)[j,jj]*u[jj]+BB(t)[j,jj]*hist[jj]
        end
    end
    for j=1:dim
    du[j] = dutemp[j]
    end
end


h(p, t) = ones(dim)*(1+im)
h(p, t)=[1.0+3.0*t 2.0]

tau=2*pi;
kappa=0.1;
delta=1.7;
epsilon=1.0;
b=0.25;
k=1;
T=k*tau;
omega=2*pi/T;

lags = [tau]

AA=Ai(v)
BB=Bi(v)

function Ai(v)
    return(t->A(t,v))
end

function Bi(v)
    return(t->B(t,v))
end

v=[kappa, delta, epsilon, 2*pi/(k*tau),b]

tspan = (0.0,5*k*tau)
p = (kappa,delta,epsilon,omega,b,tau,AA,BB)

u0 = [1.0, 2.0]
prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)

alg = MethodOfSteps(Tsit5())

sol = solve(prob,alg,abstol=1e-6,reltol=1e-6,maxiters=Int(1e6))

Plots.plot(sol)

n=30
dtt=tau/(n-1)
tvect0t=collect(-tau:dtt:tau*k+0.0000001*dtt) #discretized time vector
solt=zeros(ComplexF64,k*(n-1),dim) #empty solution vector
sol00=randn!(rng, zeros(n,dim))
sol0=hcat(tvect0t,vcat(sol00,solt))

function bc_model(du,u,h,p,t)
    kappa,delta,epsilon,omega,b,tau,AA,BB = p
    hist = h(p, t-tau)
    dutemp=zeros(ComplexF64,dim)
    for j=1:dim
        for jj=1:dim
        dutemp[j] = dutemp[j]+AA(t)[j,jj]*u[jj]+BB(t)[j,jj]*hist[jj]
        end
    end
    for j=1:dim
    du[j] = dutemp[j]
    end
end
prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)


function fsol(v,tau,k,IF)
    alg = MethodOfSteps(Tsit5())
    lags = [tau]
    p = (v[1],v[2],v[3],v[4],v[5],tau,Ai(v),Bi(v))
    interp=it(IF)
    h(p,t)=sub(interp,t)
    tspan = (0,k*tau)
    u0 = sub(interp,0)
    prob = DDEProblem(bc_model,u0,h,tspan,p; constant_lags=lags)
    return(solve(prob,alg,abstol=1e-6,reltol=1e-6,maxiters=Int(1e6)))
end

sol1=fsol(v,tau,k,sol0)
Plots.plot(sol1)
sol1(0.0:0.1:1.0)[1:2,1:end]

sol0m[n+j+1,2:2+(dim-1)]=
