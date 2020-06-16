using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using Plots
pyplot()
using PyPlot
pygui(true);

SDtest=[20.0 20.0 100.0 0.04039099408 0.0009545288260737481 100.0; 20.0 20.0 130.0 0.085684214 0.0005648921911916016 130.0; 20.0 20.0 160.0 0.14415177796 0.00037294375148073174 160.0; 20.0 20.0 190.0 0.19710940594 0.0002644804072319586 190.0; 20.0 20.0 220.0 0.28288800384 0.00019727225237595907 220.0; 20.0 20.0 250.0 0.41294403996 0.00015277006872441035 250.0; 20.0 20.0 280.0 0.64063276398 0.00012178867829380278 280.0]

SDsptest=[20.0 20.0 100.0 0.030088751457142862 0.0009545288260744258 100.0; 20.0 20.0 130.0 0.04386523567142857 0.0005648921911907693 130.0; 20.0 20.0 160.0 0.07785652 0.000372943751477933 160.0; 20.0 20.0 190.0 0.10671704007142856 0.00026448040723235053 190.0; 20.0 20.0 220.0 0.1441917871857143 0.00019727225238061892 220.0; 20.0 20.0 250.0 0.18978978264285717 0.00015277006872479394 250.0; 20.0 20.0 280.0 0.2790684242285715 0.00012178867829161501 280.0]

SDfulltest=[20.0 20.0 100.0 0.10177594202000001 0.0009545288260757554 100.0; 20.0 20.0 130.0 0.252632326 0.0005648921911906949 130.0; 20.0 20.0 160.0 0.44659127792 0.0003729437514799161 160.0; 20.0 20.0 190.0 0.9057419859 0.000264480407234677 190.0; 20.0 20.0 220.0 1.47638074598 0.00019727225237897085 220.0; 20.0 20.0 250.0 2.11054196988 0.00015277006872672019 250.0; 20.0 20.0 280.0 3.0442554518400002 0.00012178867829221086 280.0]

Plots.plot(SDtest[:,3],hcat(SDtest[:,4],SDsptest[:,4],SDfulltest[:,4]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["SDM" "SDM sparse matrices" "SDM full"],markershapes= [:circle],xlabel = "number of steps",ylabel = "T_comp")


SDtest=[20.0 20.0 100.0 0.07924058803999999 0.0009545288260757554 100.0; 20.0 20.0 130.0 0.13934559187999998 0.0005648921911906949 130.0; 20.0 20.0 160.0 0.3242773879 0.0003729437514799161 160.0; 20.0 20.0 190.0 1.0092619919599999 0.000264480407234677 190.0; 20.0 20.0 220.0 1.60632839998 0.00019727225237897085 220.0; 20.0 20.0 250.0 2.2122149360199996 0.00015277006872672019 250.0; 20.0 20.0 280.0 3.3776225240999995 0.00012178867829221086 280.0]

SDISI_1=[4.0 5.0 100.0 0.06203822861333335 0.0006988645520868372 100.0; 4.0 5.0 130.0 0.13219355999333332 0.00035607460554530244 130.0; 4.0 5.0 160.0 0.2586420046866667 0.00023104947531045917 160.0; 4.0 5.0 190.0 0.7508639646733333 0.00015500397821371188 190.0; 4.0 5.0 220.0 1.3826498433533332 0.0002228389146481766 220.0; 4.0 5.0 250.0 1.9913120780199998 0.0002135201851868623 250.0; 4.0 5.0 280.0 2.9603898666866666 0.000229478003331117 280.0]

SDISI_2=[8.0 10.0 100.0 0.06714141198666666 0.000954528825627201 100.0; 8.0 10.0 130.0 0.13417358199333335 0.00056489219068048 130.0; 8.0 10.0 160.0 0.3243712113333333 0.00037294375090853217 160.0; 8.0 10.0 190.0 0.9317269639733333 0.000264480406582696 190.0; 8.0 10.0 220.0 1.5230886667199999 0.00019727225190899455 220.0; 8.0 10.0 250.0 2.067666549953333 0.00015277006810723127 250.0; 8.0 10.0 280.0 3.1028980927933327 0.00012178867779523654 280.0]

Plots.plot(SDtest[:,3],hcat(SDtest[:,4],SDISI_1[:,4],SDISI_2[:,4]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["SDM" "SDM ISI 1" "SDM ISI 2"],markershapes= [:circle],xlabel = "number of steps",ylabel = "T_comp")

#RK vs SD comparison

ISI_RK3_quad=[12.0 12.0 100.0 0.23014466030000005 0.00031210529994779373 3564.0; 12.0 12.0 200.0 0.48674575000000003 3.8378513663681676e-5 7164.0; 12.0 12.0 300.0 0.6855867401 1.1308499479840472e-5 10764.0; 12.0 12.0 400.0 0.9054576199 4.757345784885874e-6 14364.0; 12.0 12.0 500.0 1.1091433299000002 2.4315143379471994e-6 17964.0; 12.0 12.0 600.0 1.3117898399999999 1.4053946934296056e-6 21564.0; 12.0 12.0 700.0 1.5718381799999999 8.841781230909413e-7 25164.0; 12.0 12.0 800.0 1.7319335000999998 5.918440816132366e-7 28764.0; 12.0 12.0 900.0 1.9717320400999996 4.1535793498223545e-7 32364.0; 12.0 12.0 1000.0 2.1465615401 3.0257389980958355e-7 35964.0; 12.0 12.0 1100.0 2.3549832998 2.271585030333033e-7 39564.0; 12.0 12.0 1200.0 2.54855922 1.7483252362936412e-7 43164.0; 12.0 12.0 1300.0 2.8132279302 1.3739455669622175e-7 46764.0; 12.0 12.0 1400.0 3.1051707198999994 1.0990479272545953e-7 50364.0; 12.0 12.0 1500.0 3.3204282501 8.926668736465982e-8 53964.0]

Plots.plot(vcat(SDfulltest[:,3],SDsptest[:,3],ISI_RK3_quad[:,3]),vcat(SDfulltest[:,4],SDsptest[:,4],ISI_RK3_quad[:,4]),markershapes= [:circle],xscale=:log10, yscale=:log10)

open("temp.txt", "w") do io
           writedlm(io, SDfulltest)
       end

Plots.savefig("1.pdf")
