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

function test(gmaxtest1,ntest1,multtest1,rep)
    valref=valrefA
    gmaxn=size(gmaxtest1)[1]
    ntestn=size(ntest1)[1]
    multtestn=size(multtest1)[1]
    table=zeros(Float64,multtestn*ntestn*gmaxn,5)

    tabtempcomp=zeros(Float64,rep)
    tabtemperr=zeros(ComplexF64,rep)

    for j1=1:multtestn
        global mult=multtest1[j1]
        for j2=1:gmaxn
            global gmax=gmaxtest1[j2]
            for j3=1:ntestn
                global n=ntest1[j3]
                for j4=1:rep
                    valtemp=@timed ISIM(v)
                    tabtemperr[j4]=normmax((valtemp[1])[:,end])
                    tabtempcomp[j4]=valtemp[2]
                end
                vecttempabs=real.(tabtemperr)+im*abs.(imag.(tabtemperr))
                tableerr=real(mean(abs.(vecttempabs-valref*ones(ComplexF64,rep,1))))
                tablecomp=mean(tabtempcomp)
                table[(j1-1)*ntestn*gmaxn+(j2-1)*ntestn+j3,:]=[convert(Float64,multtest1[j1]) convert(Float64,gmaxtest1[j2]) convert(Float64,ntest1[j3]) tablecomp tableerr]
            end
        end
    end
    return(table)
end

gmax=16
mult=32
n=2002
method="RK4"
valrefAn=normmax((ISIM(v))[:,end])
valrefA=real(valrefAn)+im*abs(imag.(valrefAn))

gmaxtest=collect(2:1:3)
ntest=collect(100:100:1500)
multtest=collect(4:1:6)

text=test(gmaxtest,ntest,multtest,2)

open("temp.txt", "w") do io
            writedlm(io,text)
end
