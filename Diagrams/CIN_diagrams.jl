using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using Plots

#fixed step

BS3_m=[12.0 12.0 100.0 5.4202212 0.00032193522912318264 3648.0; 12.0 12.0 200.0 1.2652883 4.0444422070497775e-5 7248.0; 12.0 12.0 300.0 1.9040851 1.1883386496116983e-5 10848.0; 12.0 12.0 400.0 2.4571636 5.066546938411943e-6 14448.0; 12.0 12.0 500.0 3.427779 2.5898418896298795e-6 18012.0; 12.0 12.0 600.0 3.9240238 1.500951432448713e-6 21612.0; 12.0 12.0 700.0 4.8960819 9.397261526301976e-7 25248.0; 12.0 12.0 800.0 5.40934 6.3390548985412e-7 28812.0; 12.0 12.0 900.0 6.1894893 4.4471452222007574e-7 32448.0; 12.0 12.0 1000.0 6.8594777 3.2147508997492525e-7 36048.0; 12.0 12.0 1100.0 7.358373599 2.431380165367002e-7 39612.0; 12.0 12.0 1200.0 7.9715474 1.8653965203696492e-7 43248.0; 12.0 12.0 1300.0 8.778891601 1.473755106990615e-7 46848.0; 12.0 12.0 1400.0 9.460292901 1.1784633400977696e-7 50412.0; 12.0 12.0 1500.0 10.371734099 9.540714137450083e-8 54048.0]

RK4_m=[12.0 12.0 100.0 0.9340539 7.2277943917968625e-6 4860.0; 12.0 12.0 200.0 1.730779401 4.47524484730275e-7 9660.0; 12.0 12.0 300.0 2.5363204 8.777491493501288e-8 14460.0; 12.0 12.0 400.0 3.410826499 2.7730702964979433e-8 19260.0; 12.0 12.0 500.0 4.2131247 1.1730752454589376e-8 24012.0; 12.0 12.0 600.0 5.2310312 5.5796024739440126e-9 28812.0; 12.0 12.0 700.0 6.219534199 2.970673835096275e-9 33660.0; 12.0 12.0 800.0 7.1061633 1.7751120563458198e-9 38412.0; 12.0 12.0 900.0 7.909546399 1.0990850870584358e-9 43260.0; 12.0 12.0 1000.0 8.8606973 7.154252378157549e-10 48060.0; 12.0 12.0 1100.0 9.502960101 4.949964144825488e-10 52812.0; 12.0 12.0 1200.0 10.637608399 3.47940568224135e-10 57660.0; 12.0 12.0 1300.0 11.207043801 2.659892331067075e-10 62460.0; 12.0 12.0 1400.0 12.131092301 1.873695193243816e-10 67212.0; 12.0 12.0 1500.0 12.977660299 1.3987905984266176e-10 72060.0]

T5_m=[12.0 12.0 100.0 12.515948101 4.8909812517218314e-9 7284.0; 12.0 12.0 200.0 2.8347748 1.3897814258822834e-10 14484.0; 12.0 12.0 300.0 4.585743 1.6731976927965135e-11 21684.0; 12.0 12.0 400.0 6.1060388 4.1340731592700695e-12 28884.0; 12.0 12.0 500.0 7.5595034 1.5882663133308481e-12 36012.0; 12.0 12.0 600.0 8.8426515 1.0725467706705903e-12 43212.0; 12.0 12.0 700.0 10.378860101 1.015778721241924e-12 50484.0; 12.0 12.0 800.0 11.8651029 8.09056349429948e-13 57612.0; 12.0 12.0 900.0 12.9144848 1.056536530130928e-12 64884.0; 12.0 12.0 1000.0 14.776865899 1.0300775823974865e-12 72084.0; 12.0 12.0 1100.0 15.880218601 7.509932363253416e-13 79212.0; 12.0 12.0 1200.0 17.5985539 1.0355907571121297e-12 86484.0; 12.0 12.0 1300.0 19.61715 1.036608882560477e-12 93684.0; 12.0 12.0 1400.0 20.279926999 7.54116208806188e-13 100812.0; 12.0 12.0 1500.0 21.933035 9.96789003922704e-13 108084.0]


Plots.plot(BS3_m[:,3],hcat(BS3_m[:,5],RK4_m[:,5],T5_m[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3" "RK4" "T5"],markershapes= [:circle],xlabel = "segment number",ylabel = "eigenvalue error")

Plots.plot(hcat(BS3_m[:,6],RK4_m[:,6],T5_m[:,6]),hcat(BS3_m[:,5],RK4_m[:,5],T5_m[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3" "RK4" "T5"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#conventional tolerance (reltol=0)
BS3_tol=[12.0 12.0 0.001 0.28951762 0.03306764949962772 1213.8; 12.0 12.0 0.0001 0.5328538798 0.005324081829813907 2613.0; 12.0 12.0 1.0e-5 0.8801450599999999 0.0007982988800392238 4672.2; 12.0 12.0 1.0e-6 1.55659968 0.00011333906960821929 8472.0; 12.0 12.0 1.0e-7 3.04644494 1.6491257787971777e-5 16615.8; 12.0 12.0 1.0e-8 5.807880120399999 2.36493367372889e-6 32526.6; 12.0 12.0 1.0e-9 11.27919142 3.4397507983519823e-7 63676.2; 12.0 12.0 1.0e-10 21.9591871596 4.909193573110624e-8 126045.0]

RK4_tol=[12.0 12.0 0.001 0.6209394 0.0050121672017777626 3043.2; 12.0 12.0 0.0001 0.94959136 0.0007919633897420541 4868.4; 12.0 12.0 1.0e-5 1.4091622797999999 0.00010864063171245469 7482.0; 12.0 12.0 1.0e-6 2.4086187198 1.4966124693496467e-5 12690.0; 12.0 12.0 1.0e-7 3.9830188 2.0538866516260254e-6 21865.2; 12.0 12.0 1.0e-8 6.6694332398 2.757012811869218e-7 35253.6; 12.0 12.0 1.0e-9 10.52901362 3.6489037301457284e-8 56664.0; 12.0 12.0 1.0e-10 17.508842939799997 4.931064726020952e-9 93326.4]

T5_tol=[12.0 12.0 0.001 0.2809635802 0.006932634357487853 1042.8; 12.0 12.0 0.0001 0.37312362 0.000778220322144805 1610.4; 12.0 12.0 1.0e-5 0.6123046598 6.240231848948544e-5 2745.6; 12.0 12.0 1.0e-6 0.8928662402 5.53088073080319e-6 4393.2; 12.0 12.0 1.0e-7 1.3388432397999999 4.434844152332751e-7 6859.2; 12.0 12.0 1.0e-8 1.9976061999999999 4.0852617461114045e-8 10318.8; 12.0 12.0 1.0e-9 3.0767641595999997 4.0784336072355435e-9 15247.2; 12.0 12.0 1.0e-10 4.53864224 4.694373700260091e-10 22840.8]

Plots.plot(BS3_tol[:,3],hcat(BS3_tol[:,5],RK4_tol[:,5],T5_tol[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 tol" "RK4 tol" "T5 tol"],markershapes= [:circle],xlabel = "absolute tolerance",ylabel = "eigenvalue error")

Plots.plot(hcat(BS3_tol[:,6],RK4_tol[:,6],T5_tol[:,6]),hcat(BS3_tol[:,5],RK4_tol[:,5],T5_tol[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 tol" "RK4 tol" "T5 tol"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#array tolerance (reltol=0)

BS3_arr_tol=[12.0 12.0 0.001 1.8277881598 0.0731183259384573 1301.4; 12.0 12.0 0.0001 0.45149488 0.0116540271353403 2046.6; 12.0 12.0 1.0e-5 0.81888244 0.0019407721508028465 3973.2; 12.0 12.0 1.0e-6 1.3530057404 0.0003268532679074574 6765.0; 12.0 12.0 1.0e-7 2.3702864200000002 5.5364897888614026e-5 12260.2; 12.0 12.0 1.0e-8 4.3512787602000005 9.631741456369355e-6 21928.0; 12.0 12.0 1.0e-9 7.815700619999999 1.6954619685836212e-6 40513.0; 12.0 12.0 1.0e-10 13.801609379800002 3.023844258412851e-7 72286.0]

RK4_arr_tol=[12.0 12.0 0.001 2.5699097197999996 0.007008883226146747 2848.8; 12.0 12.0 0.0001 0.9539130199999999 0.0012323994896607116 4448.4; 12.0 12.0 1.0e-5 1.4166059404 0.0001943882881934563 6904.8; 12.0 12.0 1.0e-6 2.25520018 2.9955474822033787e-5 11175.6; 12.0 12.0 1.0e-7 3.3538978800000003 4.261272302653382e-6 16901.8; 12.0 12.0 1.0e-8 5.5642829396 6.605466897314098e-7 27499.0; 12.0 12.0 1.0e-9 8.7265345802 1.0688583828737712e-7 44348.2; 12.0 12.0 1.0e-10 14.146337259800001 1.6522192685468916e-8 71915.8]

T5_arr_tol=[12.0 12.0 0.001 0.3333618796 0.010784893152152104 1292.4; 12.0 12.0 0.0001 0.44186272019999995 0.0013741630232497015 1803.6; 12.0 12.0 1.0e-5 0.6800396806 0.0001661523172588209 2978.4; 12.0 12.0 1.0e-6 0.9029716198000001 1.5038114737138623e-5 4172.4; 12.0 12.0 1.0e-7 1.1669062397999999 1.1233963374816997e-6 5536.6; 12.0 12.0 1.0e-8 1.78867008 1.2024753044378075e-7 8342.2; 12.0 12.0 1.0e-9 2.4814776198 1.2684229972491945e-8 12293.8; 12.0 12.0 1.0e-10 3.7832345004000003 1.4610265467078433e-9 18453.4]

Plots.plot(BS3_arr_tol[:,3],hcat(BS3_arr_tol[:,5],RK4_arr_tol[:,5],T5_arr_tol[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 tol" "RK4 tol" "T5 tol"],markershapes= [:circle],xlabel = "absolute tolerance",ylabel = "eigenvalue error")

Plots.plot(hcat(BS3_arr_tol[:,6],RK4_arr_tol[:,6],T5_arr_tol[:,6]),hcat(BS3_arr_tol[:,5],RK4_arr_tol[:,5],T5_arr_tol[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 tol" "RK4 tol" "T5 tol"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

Plots.plot(hcat(BS3_tol[:,6],RK4_tol[:,6],T5_tol[:,6],BS3_arr_tol[:,6],RK4_arr_tol[:,6],T5_arr_tol[:,6]),hcat(BS3_tol[:,5],RK4_tol[:,5],T5_tol[:,5],BS3_arr_tol[:,5],RK4_arr_tol[:,5],T5_arr_tol[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 tol" "RK4 tol" "T5 tol" "BS3 arr tol" "RK4 arr tol" "T5 arra tol"],markershapes= [:circle],xlabel = "absolute tolerance",ylabel = "eigenvalue error")

#TIME compare (reltol=0)

RK4_arr=[12.0 12.0 0.001 0.6077114997999999 0.008466023285516591 2829.0; 12.0 12.0 0.0001 0.89561984 0.0012615482919412545 4500.6; 12.0 12.0 1.0e-5 1.3558341000000003 0.00019041863584380462 6877.8; 12.0 12.0 1.0e-6 2.2463938701 2.9992312313401096e-5 11074.8; 12.0 12.0 1.0e-7 3.2164193699999997 4.346965798567743e-6 16784.8; 12.0 12.0 1.0e-8 5.1768124298 6.74137651810476e-7 27201.4; 12.0 12.0 1.0e-9 8.32715561 1.0599516993049487e-7 44340.4; 12.0 12.0 1.0e-10 13.318502900000002 1.672304205691071e-8 71527.6]

RK4_simple=[12.0 12.0 0.001 0.6902720197 0.004982352558851089 3083.4; 12.0 12.0 0.0001 0.9571275502999999 0.000765013618888943 4932.0; 12.0 12.0 1.0e-5 1.4270935599999999 0.00010893851719490713 7440.6; 12.0 12.0 1.0e-6 2.4241668902 1.4861193221475773e-5 12885.6; 12.0 12.0 1.0e-7 3.9329043806 2.02936458556892e-6 21666.0; 12.0 12.0 1.0e-8 6.38555144 2.7363374518493543e-7 35105.4; 12.0 12.0 1.0e-9 10.1119114999 3.680444648238052e-8 56299.8; 12.0 12.0 1.0e-10 16.48581269 4.941768168497353e-9 93759.6]

Plots.plot(hcat(RK4_simple[:,4],RK4_arr[:,4]),hcat(RK4_simple[:,5],RK4_arr[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK4 simple" "RK4 array"],markershapes= [:circle],xlabel = "CPU time",ylabel = "eigenvalue error")

#reltol=abstol
println(A_RK4_T)
BS3_tolE=[100.0 1.0e-5 0.001 4.4849535000000005 0.009936011644361632 1663.5; 100.0 1.0e-5 0.0001 0.6061252505 0.001463079911687629 2892.0; 100.0 1.0e-5 1.0e-5 1.04411735 0.00020232639421881703 5236.5; 100.0 1.0e-5 1.0e-6 1.8690299505 2.8808763914268573e-5 9607.5; 100.0 1.0e-5 1.0e-7 3.5636830505 4.168641751292576e-6 17536.5; 100.0 1.0e-5 1.0e-8 6.6818386995 5.966206086112762e-7 33207.0; 100.0 1.0e-5 1.0e-9 12.8468287 8.443813749720747e-8 64272.0; 100.0 1.0e-5 1.0e-10 24.515339799499998 1.1942692776907183e-8 125098.5]

RK4_tolE=[100.0 1.0e-5 0.001 2.2609939 0.0007801692093169503 3621.0; 100.0 1.0e-5 0.0001 1.2691573505 0.0001296609827582439 5643.0; 100.0 1.0e-5 1.0e-5 1.7863811004999999 1.763814177265187e-5 8802.0; 100.0 1.0e-5 1.0e-6 3.0876745999999997 2.3281260662131153e-6 14196.0; 100.0 1.0e-5 1.0e-7 5.1822985500000005 3.0451950934990216e-7 22989.0; 100.0 1.0e-5 1.0e-8 7.6238895499999995 3.978185207225552e-8 37485.0; 100.0 1.0e-5 1.0e-9 12.22855955 5.224183243803339e-9 60975.0; 100.0 1.0e-5 1.0e-10 20.0156892005 7.012891281827744e-10 98082.0]

T5_tolE=[100.0 1.0e-5 0.001 3.01637175 0.0007277444408565697 1551.0; 100.0 1.0e-5 0.0001 0.5236946505 9.985131573409602e-5 2205.0; 100.0 1.0e-5 1.0e-5 0.7232553500000001 7.910562523421285e-6 3090.0; 100.0 1.0e-5 1.0e-6 0.9641905 7.013071509609804e-7 4428.0; 100.0 1.0e-5 1.0e-7 1.390708301 6.147529526315836e-8 6570.0; 100.0 1.0e-5 1.0e-8 2.02968875 5.703674490367212e-9 9861.0; 100.0 1.0e-5 1.0e-9 3.1379602495 5.853450143690869e-10 15042.0; 100.0 1.0e-5 1.0e-10 5.311606250500001 6.834935732329056e-11 22320.0]

Plots.plot(BS3_tolE[:,3],hcat(BS3_tolE[:,5],RK4_tolE[:,5],T5_tolE[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 tol" "RK4 tol" "T5 tol"],markershapes= [:circle],xlabel = "absolute tolerance",ylabel = "eigenvalue error")

Plots.plot(hcat(BS3_tolE[:,6],RK4_tolE[:,6],T5_tolE[:,6]),hcat(BS3_tolE[:,5],RK4_tolE[:,5],T5_tolE[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 tol" "RK4 tol" "T5 tol"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

Plots.savefig("1.pdf")