using TimerOutputs
using Printf
using Statistics
using LinearAlgebra
using Interpolations
using Random
using DifferentialEquations
using QuadGK
using Plots
pyplot()
using PyPlot
pygui(true);

rng = MersenneTwister(1234)


function it(A) #creating complex iteration array
    inttyp=inttyp0
    #inttyp=BSpline(Quadratic(Line(OnGrid())))
    #inttyp=BSpline(Cubic(Line(OnGrid())))
    #inttyp=BSpline(Linear())
    matrdim=size(A,2)-1
    step=abs(real(A[end,1]-A[1,1]))/(size(A,1)-1)
    scaleitp=real(A[1,1]):step:real(A[end,1])
    ARe=real(A); AIm=imag(A);
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


function sub(it,t) #substitution in interpolating function
    subdim=size(it,1)
    # deltat=it[1].ranges[1][2]-it[1].ranges[1][1]
    # t=t-deltat/2
    out=zeros(ComplexF64,subdim)
    tend=(it[1].ranges[1])[end]
    tstart=(it[1].ranges[1])[1]
    if t>tend
        t=tend-1e-16
        println("Interpolation adjusted at end of range!")
    elseif t<tstart
        t=tstart+1e-16
        println("Interpolation adjusted at start of range!")
    end
    for j = 1:subdim
        out[j]=it[j,1](t)+it[j,2](t)*im
    end
    return(out)
end

f(t)=sin(t)
Dt=0.3
array=hcat(collect(0:Dt:2*pi),f.(collect(0:Dt:2*pi)))
T1=collect(0:Dt:2*pi)[end]
itp=it(array)
fit(t)=real(sub(itp,t))[1]

inttyp0=BSpline(Constant())
INTk=quadgk(t -> abs(fit(t)),0.0,T1,rtol=1e-10,atol=1-10)[1]
inttyp0=BSpline(Linear())
INTl=quadgk(t -> abs(fit(t)),0.0,T1,rtol=1e-10,atol=1-10)[1]
inttyp0=BSpline(Quadratic(Line(OnGrid())))
INTq=quadgk(t -> abs(fit(t)),0.0,T1,rtol=1e-10,atol=1-10)[1]
inttyp0=BSpline(Cubic(Line(OnGrid())))
INTc=quadgk(t -> abs(fit(t)),0.0,T1,rtol=1e-10,atol=1-10)[1]

#0.3
INTk
INTl
INTq
INTc

#0.1
INTk
INTl
INTq
INTc

T1=collect(0:Dt:2*pi)[end]
