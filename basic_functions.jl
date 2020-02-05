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

###################################################################################

function f(t,y,ytau,v1) #right-hand side for TNS
    A(t,v1)*y+B(t,v1)*ytau
end

function Amult(t,v1,mult1)
    Atemp=zeros(ComplexF64,mult1*dim,mult1*dim)
    for i1=0:dim:(mult1-1)*dim
        for j=1:dim
            for jj=1:dim
            Atemp[i1+j,i1+jj] =  A(t,v1)[j,jj]
            end
        end
    end
    return(Atemp)
end

function Bmult(t,v1,mult1)
    Btemp=zeros(ComplexF64,mult1*dim,mult1*dim)
    for i1=0:dim:(mult1-1)*dim
        for j=1:dim
            for jj=1:dim
            Btemp[i1+j,i1+jj] =  B(t,v1)[j,jj]
            end
        end
    end
    return(Btemp)
end

function fmult(t,y,ytau,v1,mult1) #right-hand side for TNS
    Amult(t,v1,mult1)*y+Bmult(t,v1,mult1)*ytau
end



#auxiliar functions definitions for DE.jl
function Ai(v1) #coefficient matrix A DE.jl
    return(t->A(t,v1))
end

function Bi(v1) #coefficient matrix B DE.jl
    return(t->B(t,v1))
end

function bc_model(du,u,h,p,t) #DE.jl problem definiton
    tau,AA,BB,mult1 = p
    hist = h( out1,p, t-tau)
    dutemp=zeros(ComplexF64,mult*dim)
    for i1=0:dim:(mult-1)*dim
        for j=1:dim
            for jj=1:dim
            dutemp[j+i1] = dutemp[j+i1]+AA(t)[j,jj]*u[jj+i1]+BB(t)[j,jj]*out1[jj+i1]
            end
        end
    end
    for j=1:mult*dim
    du[j] = dutemp[j]
    end
end

BaRK4=[0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0]
BbRK4=[1/6, 1/3, 1/3, 1/6]
BcRK4=[0, 0.5, 0.5, 1]
BRK4=(BaRK4,BbRK4,BcRK4)

BaRK3=[0 0 0; 0.5 0 0;-1.0 2.0 0;]
BbRK3=[1/6, 2/3, 1/6]
BcRK3=[0, 0.5, 1.0]
BRK3=(BaRK3,BbRK3,BcRK3)

BaRK2=[0 0; 0.5 0]
BbRK2=[0, 1.0]
BcRK2=[0, 0.5]
BRK2=(BaRK2,BbRK2,BcRK2)

BaEE=[0.0]
BbEE=[1.0]
BcEE=[0.0]
BEE=(BaEE,BbEE,BcEE)

################## General functions ###############
function it(A) #creating complex iteration array
    #inttyp=BSpline(Linear(Line(OnGrid())))
    inttyp=BSpline(Quadratic(Reflect(OnCell())))
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
    out=zeros(ComplexF64,subdim)
    for j = 1:subdim
        out[j]=it[j,1](t)+it[j,2](t)*im
    end
    return(out)
end

function normmax(a) #finding the greatest eigenvalue
    s=maximum(size(a))
    norma=zeros(s)
    for j = 1:s
        norma[j]=norm(a[j])
    end
    return(a[findmax(norma)[2]])
end

function normmin(a) #finding the smallest eigenvalue
    s=maximum(size(a))
    norma=zeros(s)
    for j = 1:s
        norma[j]=norm(a[j])
    end
    return(a[findmin(norma)[2]])
end

function infnancheck(a) #checking inversion feasability
    var=false
    if any(isnan,a)==true || any(isinf,a)==true
    var=true
    end
return(var)
end

function evplot(a,k) #eigenvalue plot with circle
    U=a[:,k]
    eigRe1=[real(U[n,1]) for n in 1:size(U,1)]
    eigIm1=[imag(U[n,1]) for n in 1:size(U,1)]
    plot1=Plots.plot(eigRe1,eigIm1,seriestype=:scatter)
    Re_circle=[cos(n) for n in 0:0.05:2*pi]
    Im_circle=[sin(n) for n in 0:0.05:2*pi]
    return(Plots.plot!(plot1,Re_circle,Im_circle))
end
