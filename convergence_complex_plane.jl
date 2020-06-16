using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using Plots
using TimerOutputs
pyplot()
using PyPlot
pygui(true);
using Revise
includet("ISIM_DE_sol.jl")
includet("ISIM_LMS_RK.jl")
includet("ISIM.jl")
using Main.SIDE
using Main.SILMRK
using Main.SI

includet("cont_integration_pw.jl")
includet("ISIM_CIN.jl")
using Main.ISCPW
using Main.ISCIN

omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
v=[omega kappa delta epsilon b tau0]


#USUAL PARAMETERS
abstol1=1e-4
reltol1=1e-4
ALG0=MethodOfSteps(RK4())
gmax=12
n=1000
mult=12

function convergence(ntest0,rep)
    ns=size(ntest0)[1]
    table=zeros(ComplexF64,ns,2)
    tabletemp=zeros(ComplexF64,rep)
    for j=1:ns
        for jj=1:rep
            if meth=="RK"
                valtemp=BAS.normmax(ISIM(v,(ntest0[j],gmax,mult,ALG0))[:,end])
            elseif meth=="LMS"
                valtemp=BAS.normmax(ISIM_LMS_RK(v,(ntest0[j],gmax,mult,ALG0))[:,end])
            elseif meth=="DE"
                valtemp=BAS.normmax(ISIM_DE(v,(ntest0[j],gmax,mult,ALG0),(reltol1,abstol1))[:,end])
                #adaptive
                #valtemp=BAS.normmax(ISIM_DE(v,(100,gmax,mult,ALG0),(ntest0[j],ntest0[j]))[:,end])
            elseif meth=="G"
                valtemp=BAS.normmax(ISIM_CPW(v,(100,gmax,mult,ALG0),(ntest0[j],ntest0[j]),(reltol0int,abstol0int))[:,end])
                #valtemp=BAS.normmax(ISIM_CPW(v,(ntest0[j],gmax,mult,ALG0),(reltol1,abstol1),(reltol0int,abstol0int))[:,end])
            elseif meth=="EXP"
                #valtemp=BAS.normmax(ISIM_CIN(v,(ntest0[j],gmax,mult,ALG0),(reltol1,abstol1))[:,end])
                valtemp=BAS.normmax(ISIM_CIN(v,(100,gmax,mult,ALG0),(ntest0[j],ntest0[j]))[:,end])
            end
            valtemp=real(valtemp)+im*abs(imag(valtemp))
            tabletemp[jj]=valtemp
        end
        table[j,:]=[ntest0[j] mean(tabletemp)]
    end
    return(table)
end

meth="RK"
ntest=collect(50:50:300)
ALG0=BAS.BRK3
convRK3q=convergence(ntest,2)
convRK3l=convergence(ntest,2)
convRK3k=convergence(ntest,2)

Plots.savefig("1.pdf")

refA=ISIM(v,(nlim,gmax,mult,ALG0))
refAlast=refA[:,end]

BAS.evplot(refAlast,1)
plot!(real.(hcat(convRK3q[:,2],convRK3l[:,2],convRK3k[:,2])),imag.(hcat(convRK3q[:,2],convRK3l[:,2],convRK3k[:,2])),markershapes=[:circle],label=["RK3 q" "RK3 l" "RK3 k"],xaxis=)

Plots.plot(real.(hcat(convRK3q[:,2],convRK3l[:,2],convRK3k[:,2])),imag.(hcat(convRK3q[:,2],convRK3l[:,2],convRK3k[:,2])),markershapes=[:circle],label=["RK3 q" "RK3 l" "RK3 k"],xaxis=((-0.67,-0.64)),yaxis=((0.875,0.885)))



refB=ISIM(vB,(nlim,gmax,mult,ALG0))
refBlast=refB[:,end]
BAS.evplot(refBlast,1)

meth="RK"
ntest=collect(50:50:300)
ALG0=BAS.BRK3
v=vB
convRK3Bq=convergence(ntest,2)
convRK3Bl=convergence(ntest,2)
convRK3Bk=convergence(ntest,2)
plot!(real.(hcat(convRK3Bq[:,2],convRK3Bl[:,2],convRK3Bk[:,2])),imag.(hcat(convRK3Bq[:,2],convRK3Bl[:,2],convRK3Bk[:,2])),markershapes=[:circle],label=["RK3 q" "RK3 l" "RK3 k"])

#COMPARE
ntest=collect(50:50:300)
ntest=[1e-4; 1e-5; 1e-6; 1e-7; 1e-8; 1e-9]
rep1=1
gmax=12
mult=12
nlim=20000

#point A
omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
v=[omega kappa delta epsilon b tau0]

#RK4 lin
meth="RK"
ALG0=BAS.BRK4
A_RK4_l=convergence(ntest,2)

#RK4 cub
meth="RK"
ALG0=BAS.BRK4
A_RK4_c=convergence(ntest,2)

#AB4
meth="LMS"
ALG0=BAS.BLM4
A_AB4=convergence(ntest,2)

#RK4 dense
meth="DE"
ALG0=MethodOfSteps(RK4())
A_RK4_dense=convergence(ntest,2)

#RK4 dense adapt
meth="DE"
ALG0=MethodOfSteps(RK4())
A_RK4_ad=convergence(ntest,2)

meth="G"
#RK4 Gauss fix
ALG0=MethodOfSteps(RK4())
abstol0int=1e-5
reltol0int=1e-5
A_RK4_G_fix=convergence(ntest,2)

#RK4 Gauss tol
ALG0=MethodOfSteps(RK4())
A_RK4_G_tol=convergence(ntest,2)

meth="EXP"
#RK4 combined fix
ALG0=MethodOfSteps(RK4())
A_RK4_COM_fix=convergence(ntest,2)

#RK4 combined tol
ALG0=MethodOfSteps(RK4())
A_RK4_COM_tol=convergence(ntest,2)

refA=ISIM(v,(nlim,gmax,mult,ALG0))
refAlast=refA[:,end]
BAS.evplot(refAlast,1)
Plots.plot!(real.(hcat(A_RK4_l[:,2],A_RK4_c[:,2],A_AB4[:,2],A_RK4_dense[:,2],A_RK4_ad[:,2])),imag.(hcat(A_RK4_l[:,2],A_RK4_c[:,2],A_AB4[:,2],A_RK4_dense[:,2],A_RK4_ad[:,2])),xlims=(-0.6450495,-0.645047),ylims=(0.884496,0.884501),label=["RK4 l" "RK4 c" "AB4" "RK4 dense" "RK4 ad"],markershapes=[:circle],aspect_ratio=:equal)

#save
A_RK4_l=[50.0 + 0.0im -0.6445595535199072 + 0.8824819771856698im; 100.0 + 0.0im -0.6449165831873603 + 0.8840226271970862im; 150.0 + 0.0im -0.6449885592533006 + 0.8842905616298051im; 200.0 + 0.0im -0.6450141541387961 + 0.8843828626612411im; 250.0 + 0.0im -0.6450260610571735 + 0.8844252845906502im; 300.0 + 0.0im -0.6450325418487527 + 0.8844482376051308im]

A_RK4_c=[50.0 + 0.0im -0.6451001130213107 + 0.8843814773352587im; 100.0 + 0.0im -0.6450504530725782 + 0.8844917756449595im; 150.0 + 0.0im -0.6450478358858032 + 0.8844981886885377im; 200.0 + 0.0im -0.6450474320832078 + 0.8844994012184466im; 250.0 + 0.0im -0.6450473333512814 + 0.884499772003146im; 300.0 + 0.0im -0.6450473022701722 + 0.8844999197693983im]

A_AB4=[50.0 + 0.0im -0.6460197429257889 + 0.880132304175925im; 100.0 + 0.0im -0.6451636271219802 + 0.8842276566733894im; 150.0 + 0.0im -0.6450738962473944 + 0.884446847098084im; 200.0 + 0.0im -0.6450562653179761 + 0.8844833478917484im; 250.0 + 0.0im -0.6450510971728287 + 0.88449325601224im; 300.0 + 0.0im -0.6450491649484135 + 0.8844968005494003im]

A_RK4_dense=[50.0 + 0.0im -0.6451081617106882 + 0.884407757497716im; 100.0 + 0.0im -0.6450514231571132 + 0.8844949367698762im; 150.0 + 0.0im -0.6450481201721086 + 0.8844991147311113im; 200.0 + 0.0im -0.645047551378451 + 0.8844997897677579im; 250.0 + 0.0im -0.6450473942382754 + 0.884499970304035im; 300.0 + 0.0im -0.6450473374324739 + 0.8845000342850908im]

A_RK4_ad=[0.0001 + 0.0im -0.645291222003364 + 0.8840411370410827im; 1.0e-5 + 0.0im -0.645086352204321 + 0.8844437932108455im; 1.0e-6 + 0.0im -0.6450526326964713 + 0.8844934002245792im; 1.0e-7 + 0.0im -0.6450480162034815 + 0.8844992472734826im; 1.0e-8 + 0.0im -0.6450473778278886 + 0.88449999041544im; 1.0e-9 + 0.0im -0.6450472965771835 + 0.8845000796845571im]

#save Gauss+combined

A_RK4_dense=[50.0 + 0.0im -0.6451081617106873 + 0.8844077574977158im; 100.0 + 0.0im -0.6450514231571136 + 0.8844949367698779im; 150.0 + 0.0im -0.6450481201721093 + 0.8844991147311121im; 200.0 + 0.0im -0.6450475513784506 + 0.8844997897677578im; 250.0 + 0.0im -0.6450473942382758 + 0.8844999703040348im; 300.0 + 0.0im -0.6450473374324739 + 0.8845000342850913im]

A_RK4_G_fix=[50.0 + 0.0im -0.64510981661919 + 0.8843997447229607im; 100.0 + 0.0im -0.6450512171033671 + 0.8844952424505492im; 150.0 + 0.0im -0.6450479633547566 + 0.8844991184834243im; 200.0 + 0.0im -0.6450475940333623 + 0.8844998279226539im; 250.0 + 0.0im -0.645047395410334 + 0.8844999358088361im; 300.0 + 0.0im -0.6450473339760631 + 0.8845000267472803im]

A_RK4_G_tol=[0.0001 + 0.0im -0.6452792450712539 + 0.8840650768828278im; 1.0e-5 + 0.0im -0.6450905774572314 + 0.8844447270433251im; 1.0e-6 + 0.0im -0.6450527263416238 + 0.8844933441009999im; 1.0e-7 + 0.0im -0.6450480947269199 + 0.8844992344073045im; 1.0e-8 + 0.0im -0.6450473708131358 + 0.8844999906909308im; 1.0e-9 + 0.0im -0.6450472963684999 + 0.8845000795491178im]

A_RK4_COM_fix=[50.0 + 0.0im -0.645076568903212 + 0.8843901845418866im; 100.0 + 0.0im -0.6450497716791443 + 0.8844933051800713im; 150.0 + 0.0im -0.6450478161018756 + 0.8844987673431192im; 200.0 + 0.0im -0.6450474563687458 + 0.8844996835794448im; 250.0 + 0.0im -0.6450473573369241 + 0.884499924800352im; 300.0 + 0.0im -0.6450473207326882 + 0.8845000130406092im]

A_RK4_COM_tol=[0.0001 + 0.0im -0.6449529846407223 + 0.883591593461136im; 1.0e-5 + 0.0im -0.6450575160481485 + 0.8843816827034289im; 1.0e-6 + 0.0im -0.6450510454238856 + 0.8844842322340551im; 1.0e-7 + 0.0im -0.6450479387554036 + 0.8844979368893322im; 1.0e-8 + 0.0im -0.64504737774267 + 0.8844998180551961im; 1.0e-9 + 0.0im -0.6450472976868367 + 0.8845000568980951im]

refA=ISIM_DE(v,(nlim,gmax,mult,ALG0),(reltol1,abstol1))
refAlast=refA[:,end]
BAS.evplot(refAlast,1)
Plots.plot!(real.(hcat(A_RK4_dense[:,2],A_RK4_G_fix[:,2],A_RK4_G_tol[:,2],A_RK4_COM_fix[:,2],A_RK4_COM_tol[:,2])),imag.(hcat(A_RK4_dense[:,2],A_RK4_G_fix[:,2],A_RK4_G_tol[:,2],A_RK4_COM_fix[:,2],A_RK4_COM_tol[:,2])),label=["RK4 dense" "RK4 G fix" "RK4 G tol" "RK4 com fix" "RK4 com tol"],markershapes=[:circle],aspect_ratio=:equal)
Plots.savefig("1.pdf")
xlims=(-0.6450474,-0.6450472),ylims=(0.88449975,0.8845001)
