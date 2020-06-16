using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using Plots

valrefZ1=0.6144530338830984 + 0.7019781358042942*im

Z1_RK4_dense=[12.0 12.0 100.0 5.77020954 0.002221334852928402 4764.0; 12.0 12.0 200.0 9.948554940000001 0.00043260402194843716 9564.0; 12.0 12.0 300.0 14.50691018 0.00015727322138187888 14364.0; 12.0 12.0 400.0 20.147090519400003 0.0008009405039744316 19164.0; 12.0 12.0 500.0 26.044938960000003 0.00038112910636208224 23964.0; 12.0 12.0 600.0 31.68277306 0.00011297781359395052 28764.0; 12.0 12.0 700.0 35.2966787598 1.12196636087079e-5 33612.0; 12.0 12.0 800.0 36.716746900000004 6.759048376864004e-5 38364.0; 12.0 12.0 900.0 44.62021292 0.00023847663716852226 43212.0; 12.0 12.0 1000.0 47.514161179400006 0.00026842641799718325 47964.0; 12.0 12.0 1100.0 50.6778193 3.695043478496668e-5 52812.0; 12.0 12.0 1200.0 53.57501441979999 5.818904436475282e-5 57612.0; 12.0 12.0 1300.0 59.64970154 4.392099478286532e-5 62412.0; 12.0 12.0 1400.0 64.5032650202 2.2657768531885955e-5 67164.0; 12.0 12.0 1500.0 69.8111526 0.00010938814498588108 72012.0]

Z1_RK4_dense_ad=[12.0 12.0 0.01 2.169067340000001 0.017730343186051416 1382.4; 12.0 12.0 0.001 3.1287354400000003 0.0033837759867855447 2868.0; 12.0 12.0 0.0001 4.9610001804 0.0006569878436247114 4579.2; 12.0 12.0 1.0e-5 7.55385928 0.00010695122898313669 6686.4; 12.0 12.0 1.0e-6 10.70812114 1.4773525527794878e-5 9321.6; 12.0 12.0 1.0e-7 14.3402794 1.936599197318543e-6 13567.2; 12.0 12.0 1.0e-8 20.2002883002 2.500757956057049e-7 20355.6; 12.0 12.0 1.0e-9 32.63820054 3.1485625254000436e-8 30790.8; 12.0 12.0 1.0e-10 47.684032880400004 3.971291931979007e-9 47545.2]

Z1_RK4_G_5=[100.0 1.0e-5 0.01 1.810759466333333 257.6646585975242 1764.0; 100.0 1.0e-5 0.001 2.7271731999999997 0.0042823920000986705 2826.0; 100.0 1.0e-5 0.0001 4.2258629 0.0006962517755027586 4650.0; 100.0 1.0e-5 1.0e-5 6.0081492999999995 0.010198689913873751 6770.0; 100.0 1.0e-5 1.0e-6 8.467860366333333 0.0029359473783967206 9690.0; 100.0 1.0e-5 1.0e-7 13.440883632999999 3.508927779528454e-6 14366.0; 100.0 1.0e-5 1.0e-8 20.179510299333334 3.259247564921242e-7 21030.0; 100.0 1.0e-5 1.0e-9 29.471676066333334 2.216202894775251e-6 31240.0; 100.0 1.0e-5 1.0e-10 43.7613058 0.8262747773504054 48376.0]

Z1_RK4_G_10=[100.0 1.0e-10 0.01 2.02628810066667 205.88446366599362 1768.0; 100.0 1.0e-10 0.001 4.839879266 0.0033612770772838137 2910.0; 100.0 1.0e-10 0.0001 6.167620266666667 0.0006559657844515171 4542.0; 100.0 1.0e-10 1.0e-5 7.599073800333333 0.0001078221237974052 6680.0; 100.0 1.0e-10 1.0e-6 9.6753914 1.4881974481850575e-5 9258.0; 100.0 1.0e-10 1.0e-7 14.423140466666666 1.955077083299165e-6 13792.0; 100.0 1.0e-10 1.0e-8 19.74514376666667 2.6392032740052447e-7 20816.0; 100.0 1.0e-10 1.0e-9 28.2149687 4.164570106620426e-8 31524.0; 100.0 1.0e-10 1.0e-10 44.75582820033333 3.84712178390306e-9 48160.0]

Z1_RK4_CX=[100.0 1.0e-5 0.01 1.9242923003333334 0.052823970456582804 1519.0; 100.0 1.0e-5 0.001 3.889181800333333 0.009050086186777998 3043.0; 100.0 1.0e-5 0.0001 5.4829784 0.0014331547638058591 5091.0; 100.0 1.0e-5 1.0e-5 9.122633767 0.0002329344003433715 8053.0; 100.0 1.0e-5 1.0e-6 13.295807566666667 3.261592273764551e-5 12237.0; 100.0 1.0e-5 1.0e-7 19.90717776633333 4.135717226469159e-6 17775.0; 100.0 1.0e-5 1.0e-8 30.158499600333332 5.071954893678238e-7 26465.0; 100.0 1.0e-5 1.0e-9 49.37001143333333 6.302796997970505e-8 39539.0; 100.0 1.0e-5 1.0e-10 78.31409326666667 7.895334062879467e-9 61657.0]

Z1_RK4_C=[100.0 1.0e-10 0.01 4.9835922806 0.04714952907484077 1519.4; 100.0 1.0e-10 0.001 4.1015255002 0.008964476596718205 3003.8; 100.0 1.0e-10 0.0001 6.050625040600001 0.0014441046669851467 5071.4; 100.0 1.0e-10 1.0e-5 9.0845360404 0.0002340587435233721 8090.6; 100.0 1.0e-10 1.0e-6 13.637384039599999 3.2427513418503505e-5 12235.4; 100.0 1.0e-10 1.0e-7 19.766825439999998 4.09628859648418e-6 17693.0; 100.0 1.0e-10 1.0e-8 31.8798741598 5.068126902300991e-7 26477.0; 100.0 1.0e-10 1.0e-9 48.4450386804 6.304431900633036e-8 39915.8; 100.0 1.0e-10 1.0e-10 75.85188222 7.857324113298018e-9 61466.6]

Z1_RK4_C_norm=[100.0 1.0e-5 0.01 4.0085133495 0.006420172960043231 2304.0; 100.0 1.0e-5 0.001 4.72478795 0.001314043456039956 3957.0; 100.0 1.0e-5 0.0001 6.99143565 0.00023741081311962588 6450.0; 100.0 1.0e-5 1.0e-5 9.783340450499999 3.4761350184988844e-5 9129.0; 100.0 1.0e-5 1.0e-6 14.491105300000001 4.690597933737881e-6 13365.0; 100.0 1.0e-5 1.0e-7 21.3316559 6.112996172289649e-7 19971.0]

Plots.plot(Z1_RK4_dense[:,6],Z1_RK4_dense[:,5],xscale=:log10,yscale=:log10,label="RK4 naive",markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")
Plots.plot!(hcat(Z1_RK4_dense_ad[:,6],Z1_RK4_G_5[:,6],Z1_RK4_G_10[:,6]),hcat(Z1_RK4_dense_ad[:,5],Z1_RK4_G_5[:,5],Z1_RK4_G_10[:,5]),xscale=:log10,yscale=:log10,label=["RK4 ad" "RK4 G int5" "RK4 G int10"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")
Plots.plot!(Z1_RK4_C_norm[:,6],Z1_RK4_C_norm[:,5],xscale=:log10,yscale=:log10,label="RK4 C norm",markershapes= [:circle])

Plots.plot(Z1_RK4_dense[:,4],Z1_RK4_dense[:,5],xscale=:log10,yscale=:log10,label="RK4 naive",markershapes= [:circle],xlabel = "CPU time",ylabel = "eigenvalue error")
Plots.plot!(hcat(Z1_RK4_dense_ad[:,4],Z1_RK4_G_5[:,4],Z1_RK4_G_10[:,4]),hcat(Z1_RK4_dense_ad[:,5],Z1_RK4_G_5[:,5],Z1_RK4_G_10[:,5]),xscale=:log10,yscale=:log10,label=["RK4 ad" "RK4 G int5" "RK4 G int10"],markershapes= [:circle])
Plots.plot!(Z1_RK4_C_norm[:,4],Z1_RK4_C_norm[:,5],xscale=:log10,yscale=:log10,label="RK4 C norm",markershapes= [:circle])

Plots.savefig("1.pdf")
println(Z1_RK4_C)