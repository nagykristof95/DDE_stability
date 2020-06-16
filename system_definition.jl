module SYS

using Printf
using Statistics
using LinearAlgebra
using Interpolations
using Revise

export T, tau, dim, A, B

# ################### SYSTEM DEFINITION ###################
# #function definiton : f(t)=A(t)y(t)+B(t)y(t-tau)
dim=4 #SYSTEM DIMENSION (Post-Cauchy)


#time delay
function tau(v1)
    (2*pi)/(v1[13]*v1[11])
end

#period
function T(v1)
    (2*pi)/v1[13]
end

function Mm(t,v1) #structural matrix Mm
    return([v1[1] 0; 0 v1[2]])
end

function Km(t,v1) #structural matrix Km
    return([v1[3] 0; 0 v1[4]])
end

function Sm(t,v1) #structural matrix Sm
    return([v1[5] 0; 0 v1[6]])
end

function g(phi,phi_in1,phi_out1) #sign function (working teeth)
    g_val=0
    if phi < phi_in1
        g_val=0
    elseif (phi >= phi_in1) && (phi <= phi_out1)
        g_val=1
    else
        g_val=0
    end
    return(g_val)
end

function W(t,v1) #periodic matrix
    if v1[12]==1.0 #up-milling
        aepu=v1[10]/v1[9]; phi_in=0; phi_out=acos(1-2*aepu);
    else #down-milling
        aepd=2-v1[10]/v1[9]; phi_in=acos(1-2(aepd-1)); phi_out=pi;
    end
    zint=trunc(Int,v1[11])
    phi=zeros(zint)
    for j=0:(zint-1)
    phi[j+1]=mod(t*v1[13]+j*2*pi/v1[11],2*pi)
    end
    W_f=zeros(Float64,2,2)
    for j = 0:(zint-1)
        W_f=W_f+g(phi[j+1],phi_in,phi_out)*[-sin(2*phi[j+1])-v1[8]+v1[8]*cos(2*phi[j+1]) -1-cos(2*phi[j+1])-v1[8]*sin(2*phi[j+1]); 1-cos(2*phi[j+1])-v1[8]*sin(2*phi[j+1]) sin(2*phi[j+1])-v1[8]-v1[8]*cos(2*phi[j+1])]
    end
    return ((v1[14]*v1[7]/2)*W_f)
end

function A(t,v1) #coefficient matrix A
    vcat(hcat(-inv(Mm(t,v1))*Km(t,v1),-inv(Mm(t,v1))*(Sm(t,v1)-W(t,v1))),hcat(Matrix{Float64}(I, 2, 2),zeros(2,2)))
end

function B(t,v1) #coefficient matrix B (delay)
    vcat(hcat(zeros(2,2),-inv(Mm(t,v1))*W(t,v1)),hcat(zeros(2,2),zeros(2,2)))
end

###################################################################################

################### SYSTEM DEFINITION ###################

#dim=2 #DoF of the system (after Cauchy transcription)

#omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;

#v=[1 omega  2 kappa  3 delta  4 epsilon  5 b  6 tau0]
#v=[omega kappa delta epsilon b tau0]

#time delay
# function tau(v1)
#     v1[6]
# end

#period
# function T(v1)
#     (2*pi)/v1[1]
# end

# # matrix coefficients
# function A(t,v1)
#     [-v1[2] -(v1[3]+v1[4]*cos(v1[1]*t)); 1.0 0.0]
# end

# function A(t,v1)
#     [-v1[2] -(v1[3]+v1[4]*sign(cos(v1[1]*t))); 1.0 0.0]
# end

# function B(t,v1)
#     [0.0 v1[5]; 0.0 0.0]
# end

end #module end
