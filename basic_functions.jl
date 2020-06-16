module BAS
using TimerOutputs
using Printf
using Statistics
using LinearAlgebra
using Interpolations
using Random
using DifferentialEquations
using SparseArrays
using Plots
pyplot()
using PyPlot
pygui(true);
using Revise
includet("system_definition.jl")
using Main.SYS

export bc_model,it,sub, res1, res2, evplot,butcher,out1,to,res1full,res2full,convreq,sorteigvals

to = TimerOutput()

rng = MersenneTwister(1234)
setprecision(BigFloat,128)
setrounding(BigFloat, RoundUp)

prec=Float64
precF=ComplexF64
# prec=BigFloat
# precF=Complex{BigFloat}

#inttyp0=BSpline(Constant())
#inttyp0=BSpline(Linear())
inttyp0=BSpline(Quadratic(Line(OnGrid())))
#inttyp0=BSpline(Cubic(Line(OnGrid())))

function f(t,y,ytau,v1) #right-hand side for TNS
    SYS.A(t,v1)*y+SYS.B(t,v1)*ytau
end

function Amult(t,v1,mult1)
    Atemp=sparse(zeros(precF,mult1*dim,mult1*dim))
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
    Btemp=sparse(zeros(precF,mult1*dim,mult1*dim))
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

BLM1=([-1.0,0.0],[1.0,0.0])
BLM2=([0.0,-1.0,0.0],[-1/2,3/2,0.0])
BLM3=([0.0,0.0,-1.0,0.0],[5/12,-16/12,23/12,0.0])
BLM4=([0.0,0.0,0.0,-1.0,0.0],[-9/24,37/24,-59/24,55/24,0.0])

function butcher(t,inttau,y,dt1,(Ba1,Bb1,Bc1),v1,tau1,mult1) #one step by Butcher table (explicit only!)
    s=size(Bb1)[1]
    kvec=zeros(precF,dim*mult1,s)
    Svec=zeros(precF,dim*mult1,s)
    for j=1:s
        for jj=1:s
            Svec[:,j]=Svec[:,j]+Ba1[j,jj]*kvec[:,jj]
        end
@timeit to "fmult" kvec[:,j]=dt1*fmult(t+Bc1[j]*dt1,y+Svec[:,j],sub(inttau,t+Bc1[j]*dt1-tau1),v1,mult1)
    end
    yn=y
    for j=1:s
        yn=yn+Bb1[j]*kvec[:,j]
    end
    return(yn)
end

################## General functions ###############
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
    out=zeros(precF,subdim)
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
    plot1=Plots.plot(eigRe1,eigIm1,seriestype=:scatter,aspect_ratio=true)
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
    ret=zeros(precF,dim,multvar)
    for j=1:multvar
        ret[:,j]=solvar[1+(j-1)*dim:j*dim]
    end
    return(ret)
end

function res2(solvar)
    s=size(solvar)[2]
    multvar=convert(Int,s)
    ret=zeros(precF,dim*multvar)
    for j=1:multvar
        ret[1+(j-1)*dim:j*dim]=solvar[:,j]
    end
    return(ret)
end

function res1full(solvar)
    n1=size(solvar)[1]
    multvar=convert(Int,size(solvar)[2]/dim)
    ret=zeros(precF,dim*n1,multvar)
    for j=1:multvar
        for jj=1:n1
        ret[1+(jj-1)*dim:1+(jj-1)*dim+(dim-1),j]=solvar[jj,1+(j-1)*dim:1+(j-1)*dim+(dim-1)]
        end
    end
    return(ret)
end

function res2full(solvar)
    n1=convert(Int,size(solvar)[1]/dim)
    multvar=size(solvar)[2]
    ret=zeros(precF,n1,dim*multvar)
    for j=1:multvar
        for jj=1:n1
        ret[jj,1+(j-1)*dim:1+(j-1)*dim+(dim-1)]=solvar[1+(jj-1)*dim:1+(jj-1)*dim+(dim-1),j]
        end
    end
    return(ret)
end

function convreq(eigval0,N,tol)
    mult1=size(eigval0)[1] #number of eigenvalues calcualted
    if mult1<3 return(false) end #if only one or two eigenvalues are caluclated we won't append anything
    eigval1=norm.(eigval0) #array with norms of eigenvalues
    eigvalsort=sort(eigval1,rev = true) #sorted norms
    diff=[abs(eigvalsort[j]-eigvalsort[jj]) for j=1:N, jj=1:N]
    nover0=length(diff[diff .< tol])
    nover=(nover0-N)/2
    if nover>1
        return(true)
    else
        return(false)
    end
end

function sorteigvals(eigval0)
    eigvalsort=sort(eigval0,by=norm,rev=true)
end

end #module end
