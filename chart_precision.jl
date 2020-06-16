using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using Plots
using TimerOutputs
using DifferentialEquations
pyplot()
using PyPlot
pygui(true);
using Revise
using SparseArrays
includet("ISIM_DE_sol.jl")
includet("ISIM_LMS_RK.jl")
includet("ISIM.jl")
includet("ISIM_DE_append.jl")
includet("SDM.jl")
includet("SDM_append.jl")
includet("cont_integration_pw.jl")
includet("ISIM_CIN.jl")
using Main.SIDE
using Main.SILMRK
using Main.SI
using Main.SIDEA
using Main.SD
using Main.SDapp
using Main.ISCPW
using Main.ISCIN


setprecision(BigFloat,128)
setrounding(BigFloat, RoundUp)

prec=Float64
precF=ComplexF64

# prec=BigFloat
# precF=Complex{BigFloat}

omega=1; kappa=0.2; delta=2.5; epsilon=1.0; b=0.5; tau0=2*pi;
v=prec.([omega kappa delta epsilon b tau0])

ALG0=MethodOfSteps(RK4())
gmax=12
mult=12
n=110
n=100
nref=10000

reltol0=0.0
abstol0=1e-7

ALG0=MethodOfSteps(Tsit5())
p1test=collect(0.0:1.0:20.0)
p2test=collect(0.5:0.25:1.0)

rep1=2

function chart_test_ref(p1test1,p2test1)
    p1n=size(p1test1)[1]
    p2n=size(p2test1)[1]
    table=zeros(precF,p1n*p2n)
    for j1=1:p1n
        for j2=1:p2n
            vv=prec.([omega kappa p1test1[j1] epsilon p2test1[j2] tau0])
            valrefn=BAS.normmax(ISIM_DE(vv,(nref,gmax,mult,ALG0),(reltol0,abstol0))[:,end])
            valref=real(valrefn)+im*abs.(imag(valrefn))
            table[(j1-1)*p2n+j2]=valref
        end
    end
    return(table)
end

ref=chart_test_ref(p1test,p2test)

ref=precF.([202.88070223683883 + 3.02424751907893e-13im, 204.74477083406515 + 3.195013698054085e-13im, 206.58812431061347 + 1.56002020085647e-9im, 208.41158415395725 + 1.515575859256655e-13im, 210.21591846767194 + 8.130128514860528e-14im, -8.311841039170684 + 2.825170652975828e-14im, -8.68717667701957 + 1.9567680809018384e-15im, -8.616976164617789 + 8.673702930840447e-14im, -8.223513283054135 + 1.873935034923946e-15im, -7.441630728623712 + 2.9026260561781925e-15im, -2.0355232511293075 + 4.423544863740858e-17im, 1.018176334973702 + 0.736244537956954im, 0.9224622864975563 + 1.103149134877029e-12im, 1.0838509762918753 + 1.3176167242391472im, 1.5029709367016124 + 2.0569927699863393im, -0.7574849522260372 + 1.6984127414499284im, -0.5789532143580676 + 1.1488083194614096im, -0.39595811396837466 + 0.35752862167909893im, -0.7939483981276064 + 0.3629144275953816im, -1.0615466526845705 + 0.7877050891636891im, -1.0298593257848256 + 0.3834126626554426im, -0.66308755589001 + 0.0893021040622366im, -0.11506889195847747 + 0.520930603285111im, -0.12141306430569632 + 1.1627771470080903im, -0.1185631044247557 + 1.631631132052371im, 0.8677747819132811 + 0.906734383146332im, 0.7104911904262871 + 0.5301900207374534im, 0.5322081394269339 + 0.03693290761840187im, 0.7211258824558985 + 0.6187729349556187im, 0.9062220646111344 + 1.0119256422508567im, 0.0744780273967065 + 1.4266969520611084im, 0.07642650289066288 + 1.0429312125182424im, 0.07398077700089317 + 0.5283335953468249im, 0.5260006728803431 + 0.004860826588903953im, 0.8327656755856832 + 0.21919229186362676im])

ref2=[208.41158415395785 + 2.7962614862486877e-13im, 209.3160959477841 + 2.7940930819037035e-13im, 210.21591846767257 + 2.2870594307278225e-14im, -8.223513283054212 + 2.6287999554952535e-14im, -7.893709157794279 + 3.4911309954033243e-16im, -7.441630728623645 + 8.87744738831131e-16im, 1.0838509762918687 + 1.3176167242391488im, 1.3005561047407719 + 1.70956858829995im, 1.5029709367016102 + 2.0569927699863433im, -0.7939483981276066 + 0.3629144275953819im, -0.9387037199868993 + 0.5849609869535156im, -1.0615466526845696 + 0.7877050891636889im, -0.1214130643056939 + 1.1627771470080919im, -0.12085741296980124 + 1.4073768390017678im, -0.11856310442475343 + 1.6316311320523686im, 0.7211258824558977 + 0.618772934955618im, 0.8175924965659769 + 0.8252453743277927im, 0.9062220646111404 + 1.0119256422508642im, 0.5260006728803427 + 0.004860826588903618im, 0.691360578264212 + 0.1106114940985254im, 0.8327656755856809 + 0.2191922918636292im]

ref_wide=[-8.223513283054151 + 6.1608704249316304e-15im, -7.893709157794292 + 1.9572234458142823e-14im, -7.441630728623621 + 2.106084404096542e-14im, 1.0838509762918673 + 1.31761672423915im, 1.3005561047407719 + 1.7095685882999483im, 1.5029709367016169 + 2.0569927699863495im, -0.7939483981276071 + 0.3629144275953815im, -0.9387037199869024 + 0.5849609869535163im, -1.0615466526845658 + 0.7877050891636894im, -0.12141306430569129 + 1.1627771470080923im, -0.12085741296979904 + 1.4073768390017691im, -0.11856310442475682 + 1.631631132052374im, 0.7211258824558996 + 0.6187729349556206im, 0.8175924965659797 + 0.8252453743277921im, 0.9062220646111351 + 1.0119256422508576im, 0.5260006728803425 + 0.0048608265889040275im, 0.6913605782642128 + 0.11061149409852533im, 0.8327656755856798 + 0.21919229186363023im, -0.6811104782600997 + 0.3020236307360482im, -0.7713603739818725 + 0.46597452592834854im, -0.8508199826689224 + 0.6126119474363793im, -0.3857966966255719 + 0.8593314797288635im, -0.41055361326962614 + 1.0330645249262234im, -0.4331666806648917 + 1.1921939538492603im, 0.27736244567372337 + 0.9017557041085913im, 0.300485287753907 + 1.0688273434090412im, 0.3229317363482631 + 1.2222169517905588im, 0.6514325484810453 + 0.4343746888279133im, 0.7201053050910129 + 0.5882723114752793im, 0.7839239050690565 + 0.7264517851443296im, 0.5226225393099717 + 0.010624808622531018im, 0.6432585146128283 + 0.10768378874731592im, 0.7462095918528924 + 0.2114146781298671im, -0.4472882992534196 + 0.07391648519387589im, -0.5647482705615083 + 0.027648684354738218im, -0.6619090484734284 + 0.11717471275656266im, -0.6302075668418881 + 0.2467510949761118im, -0.6968642572694477 + 0.38050376848052425im, -0.7574417037567236 + 0.49915806555758696im, -0.4729864004506068 + 0.668947238920718im, -0.5012619186761441 + 0.806004566250118im, -0.5280125714027393 + 0.9312690702248302im, -0.03699701302570264 + 0.865210666437161im, -0.036851605464390036 + 0.9999178784222988im, -0.03644201030491885 + 1.1242145076720462im, 0.4064706369337134 + 0.7154752635868085im, 0.431955500298114 + 0.8467354342595359im, 0.4568620648885753 + 0.9671190635323216im, 0.6151461953958655 + 0.3410438712489087im, 0.6683721757951733 + 0.4673100561147543im, 0.7188725245494458 + 0.5802816841894987im, 0.524340413574029 + 0.011938424421998025im, 0.6142745595824829 + 0.1007848627385824im, 0.6934700805148665 + 0.19688357756107205im, 0.2874126195533257 + 0.13549634781191316im, 0.4129678671616912 + 0.072361792789366im, 0.5126723972776506 + 0.01322336763592904im, -0.4688616301935554 + 0.07482635156234517im, -0.5586859020313967 + 0.0283956808675484im, -0.6357599737436984 + 0.11523220140784154im, -0.6018043477059453 + 0.20710431017902894im, -0.6540626840058374 + 0.3215913626560386im, -0.7030946563855786 + 0.42290933605406134im]

function chart_test(p1test1,p2test1,rep,valreftable)
    p1n=size(p1test1)[1]
    p2n=size(p2test1)[1]
    table=zeros(prec,p1n*p2n,3)

    tabtempcomp=zeros(prec,rep)
    tabtemperr=zeros(precF,rep)
    tabtempeval=zeros(prec,rep)
    reset_timer!(BAS.to)
    reset_timer!(SI.toISIM)
    reset_timer!(SILMRK.toSILMRK)
    reset_timer!(SIDE.toSIDE)
    reset_timer!(SD.toSDM)
    reset_timer!(ISCIN.toISCIN)
    for j1=1:p1n
        for j2=1:p2n
            vv=prec.([omega kappa p1test1[j1] epsilon p2test1[j2] tau0])
                for j4=1:rep
                    if meth=="RK"
                        valtemp=@timed ISIM(vv,(n,gmax,mult,ALG0))
                        if typeof(ALG0) <: MethodOfSteps
                            tempeval=TimerOutputs.ncalls(SI.toISIM["fmult_Julia"])
                            reset_timer!(SI.toISIM)
                        elseif typeof(ALG0) <: Tuple
                            tempeval=TimerOutputs.ncalls(BAS.to["fmult"])
                            reset_timer!(BAS.to)
                        end
                        elseif meth=="LMS"
                            valtemp=@timed ISIM_LMS_RK(vv,(n,gmax,mult,ALG0))
                            tempeval=TimerOutputs.ncalls(SILMRK.toSILMRK["fmult_LM"])
                            reset_timer!(SILMRK.toSILMRK)
                        elseif meth=="DE"
                            valtemp=@timed ISIM_DE(vv,(n,gmax,mult,ALG0),(reltol0,abstol0))
                            tempeval=TimerOutputs.ncalls(SIDE.toSIDE["fmult_DE"])
                            reset_timer!(SIDE.toSIDE)
                        elseif meth=="DEA"
                            valtemp=@timed ISIM_DE(vv,(n,gmaxtest1[j2],multtest1[j1],ALG0),(ntest1[j3],ntest1[j3]))
                            tempeval=TimerOutputs.ncalls(SIDE.toSIDE["fmult_DE"])
                            reset_timer!(SIDE.toSIDE)
                        elseif meth=="SD"
                            valtemp=@timed SDM(vv,n)
                            tempeval=TimerOutputs.ncalls(SD.toSDM["f_SDM"])
                            reset_timer!(SD.toSDM)
                        elseif meth=="SDsp"
                            valtemp=@timed SDM_sp(vv,n)
                            tempeval=TimerOutputs.ncalls(SD.toSDM["f_SDM"])
                            reset_timer!(SD.toSDM)
                        elseif meth=="SDfull"
                            valtemp=@timed SDM_full(vv,(n,gmax,mult,ALG0))
                            tempeval=TimerOutputs.ncalls(SD.toSDM["f_SDM"])
                            reset_timer!(SD.toSDM)
                        elseif meth=="ISI_cont"
                            valtemp=@timed ISIM_CI(vv,(abstol0,gmax,mult,ALG0))
                            tempeval=TimerOutputs.ncalls(ISCI.toISCI["fmult_CI"])
                            reset_timer!(ISCI.toISCI)
                        elseif meth=="CIN"
                            valtemp=@timed ISIM_CIN(vv,(n,gmax,mult,ALG0),(reltol0,abstol0))
                            tempeval=TimerOutputs.ncalls(ISCIN.toISCIN["fmult_CIN"])
                            reset_timer!(ISCIN.toISCIN)
                end
                    tabtemperr[j4]=BAS.normmax((valtemp[1])[:,end])
                    tabtempcomp[j4]=valtemp[2]
                    tabtempeval[j4]=tempeval
                end
                vecttempabs=real.(tabtemperr)+im*abs.(imag.(tabtemperr))
                tableerr=real(abs(mean(vecttempabs-(valreftable[(j1-1)*p2n+j2])*ones(precF,rep,1))))
                #tableerr=real(mean(abs.(vecttempabs-(valreftable[(j1-1)*p2n+j2])*ones(precF,rep,1))/abs(valreftable[(j1-1)*p2n+j2])))
                tablecomp=mean(tabtempcomp)
                tableeval=mean(tabtempeval)
                table[(j1-1)*p2n+j2,:]=[tablecomp tableerr tableeval]
        end
    end
    return(table)
end

meth="DE"
test1_ad=chart_test(p1test,p2test,rep1,ref_wide)
test1_fix=chart_test(p1test,p2test,rep1,ref_wide)
test1_ad[1,1]=0.29256853153165651
Plots.plot(hcat(test1_ad[:,1],test1_fix[:,1]),hcat(test1_ad[:,2],test1_fix[:,2]),xscale=:log10, yscale=:log10,label=["Adaptive" "Fix"],markershapes= [:circle],xlabel = "CPU time",ylabel = "absolute error")

Plots.plot(hcat(test1_ad[:,3],test1_fix[:,3]),hcat(test1_ad[:,2],test1_fix[:,2]),xscale=:log10, yscale=:log10,label=["Adaptive" "Fix"],markershapes= [:circle],xlabel = "fcalls",ylabel = "relative error")


sum(test1_ad[:,1])
sum(test1_fix[:,1])

sum(test1_ad[:,3])
sum(test1_fix[:,3])

meth="CIN"
test1_ad_abs=chart_test(p1test,p2test,rep1,ref_wide)
test1_fix_abs=chart_test(p1test,p2test,rep1,ref_wide)


Plots.plot(hcat(test1_ad_abs[:,1],test1_fix_abs[:,1]),hcat(test1_ad_abs[:,2],test1_fix_abs[:,2]),xscale=:log10, yscale=:log10,label=["Adaptive" "Fix"],markershapes= [:circle],xlabel = "CPU time",ylabel = "absolute error")
Plots.plot(hcat(test1_ad_abs[:,3],test1_fix_abs[:,3]),hcat(test1_ad_abs[:,2],test1_fix_abs[:,2]),xscale=:log10, yscale=:log10,label=["Adaptive" "Fix"],markershapes= [:circle],xlabel = "function calls",ylabel = "absolute error")
sum(test1_ad_abs[:,1])
sum(test1_fix_abs[:,1])

sum(test1_ad_abs[:,3])
sum(test1_fix_abs[:,3])

Plots.savefig("1.pdf")

###################x
#Note
#T5 jobb hibakezelés
#abstol=1e-7
#konstans n=110
