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
    h( out1, p, t-tau)
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

BaRK5=[0.0 0 0 0 0 0; 0.16666666666666666667 0 0 0 0 0; -0.21627570527696895373 0.54960903861030228706 0 0 0 0; 0.08482881411262012706 0.04162653285051884260 0.37354465303686103035 0 0 0; -0.08651098424575942561 0.37955562705964599292 0.01753570971622337002 0.35608631413655672933 0 0; -0.12499755969423778621 0.72695084642093284094 -0.38363171852137430626 0.29492374551818501854 0.32008801960982756632 0]
BbRK5=[0.07892564703041163884, 0.15537176484794180580, 0.08925647030411638840, 0.51074352969588361160, -0.30537176484794180580, 0.47107435296958836116]
BcRK5=[0.0, 1/6, 1/3, 1/2, 2/3, 5/6]
BRK5=(BaRK5,BbRK5,BcRK5)

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
    inttyp=inttyp0
    #inttyp=BSpline(Quadratic(Line(OnGrid())))
    #inttyp=BSpline(Cubic(Reflect(OnCell())))
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

function inversecontrol(A0)
    s=size(A0)[1]
    B0=A0-Array(I,s,s)
    return(round(sum(B0)))
end

function res1(solvar)
    s=size(solvar)[1]
    multvar=convert(Int,s/dim)
    ret=zeros(ComplexF64,dim,multvar)
    for j=1:multvar
        ret[:,j]=solvar[1+(j-1)*dim:j*dim]
    end
    return(ret)
end

function res2(solvar)
    s=size(solvar)[2]
    multvar=convert(Int,s)
    ret=zeros(ComplexF64,dim*multvar)
    for j=1:multvar
        ret[1+(j-1)*dim:j*dim]=solvar[:,j]
    end
    return(ret)
end
