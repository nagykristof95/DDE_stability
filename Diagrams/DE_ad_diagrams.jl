using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using Plots

#point A

A_BS3_DE_ad=[1.0e-8 0.1 1.3018628995 0.4245547347305121 5202.0; 1.0e-8 0.01 0.1172861495 0.1400004917926302 7740.0; 1.0e-8 0.001 0.14431020049999999 0.02850000085693652 11772.0; 1.0e-8 0.0001 0.15484815 0.004628180474364192 20484.0; 1.0e-8 1.0e-5 0.2252525 0.0006486867554772441 36828.0; 1.0e-8 1.0e-6 0.3612575995 9.068171027724766e-5 69102.0; 1.0e-8 1.0e-7 0.5975825 1.2575474424645247e-5 132390.0; 1.0e-8 1.0e-8 1.0429526 1.7899972673610884e-6 248670.0; 1.0e-8 1.0e-9 1.7981318000000002 3.0047090099659734e-7 443448.0; 1.0e-8 1.0e-10 3.114386799 7.778338200814407e-8 681714.0; 1.0e-8 1.0e-11 4.3498734005 2.5823185503080755e-8 984168.0; 1.0e-8 1.0e-12 5.5488158495 8.951722560658586e-9 1.278e6]

A_RK4_DE_ad=[1.0e-8 0.1 1.2985716999999999 0.08891770279187179 10872.0; 1.0e-8 0.01 0.1393931 0.018812923204422974 15768.0; 1.0e-8 0.001 0.18865875 0.0032116073629830295 23364.0; 1.0e-8 0.0001 0.214366651 0.00046487430243367066 36000.0; 1.0e-8 1.0e-5 0.3273244995 6.082867524071127e-5 59616.0; 1.0e-8 1.0e-6 0.43816374999999996 8.103743761166892e-6 98460.0; 1.0e-8 1.0e-7 0.6541939995 1.0571001118922665e-6 157644.0; 1.0e-8 1.0e-8 1.04215945 1.4084025941326644e-7 254628.0; 1.0e-8 1.0e-9 1.5020483 2.3261832783155062e-8 384840.0; 1.0e-8 1.0e-10 2.07973435 6.151593196135515e-9 538020.0; 1.0e-8 1.0e-11 2.856787899 1.6687894572289686e-9 709524.0; 1.0e-8 1.0e-12 3.77096625 5.135099513108336e-10 840132.0]

A_T5_DE_ad=[1.0e-8 0.1 1.4585632000000002 0.04690347788022724 7380.0; 1.0e-8 0.01 0.10855580000000001 0.009178721841495692 9252.0; 1.0e-8 0.001 0.13516745 0.002033728221211427 12672.0; 1.0e-8 0.0001 0.1475558495 0.00021383773445882897 16704.0; 1.0e-8 1.0e-5 0.16851850000000002 2.3701067688794252e-5 23508.0; 1.0e-8 1.0e-6 0.21464205 2.214156473093717e-6 33156.0; 1.0e-8 1.0e-7 0.2804917 2.153795829700418e-7 52452.0; 1.0e-8 1.0e-8 0.35704464999999996 2.384024852603496e-8 73692.0; 1.0e-8 1.0e-9 0.45924255 4.48221835549218e-9 100224.0; 1.0e-8 1.0e-10 0.55069205 1.4474401660623592e-9 121824.0; 1.0e-8 1.0e-11 0.6311352 4.322433219675629e-10 146196.0; 1.0e-8 1.0e-12 0.6989491505000001 1.514582460595386e-10 161352.0]

Plots.plot(A_T5_DE_ad[:,2],hcat(A_T5_DE_ad[:,4],A_RK4_DE_ad[:,4],A_BS3_DE_ad[:,4]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["T5 ad" "RK4 ad" "BS3 ad"],markershapes= [:circle],xlabel = "t_abs",ylabel = "eigenvalue error")

Plots.plot(hcat(A_T5_DE_ad[:,5],A_RK4_DE_ad[:,5],A_BS3_DE_ad[:,5]),hcat(A_T5_DE_ad[:,4],A_RK4_DE_ad[:,4],A_BS3_DE_ad[:,4]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["T5 ad" "RK4 ad" "BS3 ad"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#sampling effect
A_RK4_DE_ad_m8=[1.0e-8 1.0e-8 3.0 1.5230791600000002 1.793616523543577e-7 238852.80000000002; 1.0e-8 1.0e-8 4.0 1.18760978 1.5996314183371857e-7 243705.59999999998; 1.0e-8 1.0e-8 5.0 1.3199360402000002 1.4987564663919044e-7 249580.80000000002; 1.0e-8 1.0e-8 6.0 1.05533568 1.4427263268959933e-7 250372.80000000002; 1.0e-8 1.0e-8 7.0 1.0999734994 1.4327423322770637e-7 251524.80000000002; 1.0e-8 1.0e-8 8.0 1.2333499202 1.4228260953562246e-7 251236.80000000002; 1.0e-8 1.0e-8 9.0 1.29462138 1.4291993771266416e-7 250876.80000000002; 1.0e-8 1.0e-8 10.0 1.1274523198000002 1.4127054093337787e-7 250387.19999999998; 1.0e-8 1.0e-8 25.0 1.1581617202 1.4083719338813667e-7 251827.19999999998; 1.0e-8 1.0e-8 50.0 1.42178236 1.4325349490579134e-7 251193.59999999998; 1.0e-8 1.0e-8 100.0 1.4074006994 1.416757019948529e-7 251668.80000000002; 1.0e-8 1.0e-8 200.0 1.3277923798 1.4161129636891202e-7 251481.59999999998; 1.0e-8 1.0e-8 300.0 1.4912433604 1.421408919513572e-7 251481.59999999998; 1.0e-8 1.0e-8 400.0 1.55273964 1.4080107582722143e-7 251841.59999999998; 1.0e-8 1.0e-8 500.0 1.5271806004 1.4142290789005456e-7 249854.40000000002; 1.0e-8 1.0e-8 1000.0 1.7234010402 1.4192674027242698e-7 252129.59999999998; 1.0e-8 1.0e-8 1500.0 1.9870818402000001 1.4041974671312645e-7 252403.19999999998; 1.0e-8 1.0e-8 2000.0 2.2476683601999996 1.406500547168861e-7 251956.80000000002; 1.0e-8 1.0e-8 2500.0 2.5015841998000004 1.4172176185783652e-7 252691.19999999998; 1.0e-8 1.0e-8 3000.0 2.7018435208 1.410037061708157e-7 251740.80000000002]

B_RK4_DE_ad_m8=[1.0e-8 1.0e-8 3.0 1.19468416 1.751559702839068e-7 265276.80000000005; 1.0e-8 1.0e-8 4.0 1.3691911402 1.602024170132031e-7 271238.4; 1.0e-8 1.0e-8 5.0 1.5789598196 1.525711175826718e-7 274075.19999999995; 1.0e-8 1.0e-8 6.0 1.3485189001999998 1.5236457386720283e-7 274968.0; 1.0e-8 1.0e-8 7.0 1.8240858998 1.5168335588006153e-7 276868.80000000005; 1.0e-8 1.0e-8 8.0 1.3298400006 1.5112109105271452e-7 275918.4; 1.0e-8 1.0e-8 9.0 1.18751698 1.5160312209781738e-7 275904.0; 1.0e-8 1.0e-8 10.0 1.3574199202 1.5284112087293796e-7 275644.80000000005; 1.0e-8 1.0e-8 25.0 1.38462506 1.521578332894387e-7 275472.0; 1.0e-8 1.0e-8 50.0 1.10789616 1.5253395051232535e-7 276364.80000000005; 1.0e-8 1.0e-8 100.0 1.4540297596 1.525513423573028e-7 276134.4; 1.0e-8 1.0e-8 200.0 1.3064370598 1.522165329643591e-7 276379.19999999995; 1.0e-8 1.0e-8 300.0 1.27232214 1.5203810673978287e-7 276998.4; 1.0e-8 1.0e-8 400.0 1.3712013602 1.5187816158648295e-7 277185.6; 1.0e-8 1.0e-8 500.0 1.4475517199999999 1.5220352251307008e-7 276523.19999999995; 1.0e-8 1.0e-8 1000.0 1.8523992797999997 1.4972651262229372e-7 278467.19999999995; 1.0e-8 1.0e-8 1500.0 1.8145625002000003 1.4996440676182588e-7 278280.0; 1.0e-8 1.0e-8 2000.0 2.2151968798 1.5211176951115335e-7 276811.19999999995; 1.0e-8 1.0e-8 2500.0 2.5347850997999997 1.5229632213848378e-7 277056.0; 1.0e-8 1.0e-8 3000.0 2.5604721802 1.5153635092695076e-7 276796.80000000005]

A_RK4_DE_ad_m7=[1.0e-7 1.0e-7 3.0 0.6577359802 1.3766808658545381e-6 146808.0; 1.0e-7 1.0e-7 4.0 0.64176006 1.2439125042799195e-6 149961.59999999998; 1.0e-7 1.0e-7 5.0 0.6395842201999999 1.1654528508710022e-6 153158.40000000002; 1.0e-7 1.0e-7 6.0 0.6312105802000001 1.1286778615229172e-6 153201.59999999998; 1.0e-7 1.0e-7 7.0 0.6493082402000001 1.0869005552772467e-6 153374.40000000002; 1.0e-7 1.0e-7 8.0 0.6601956604 1.107776282491782e-6 153806.40000000002; 1.0e-7 1.0e-7 9.0 0.6909021604 1.115628753970804e-6 153187.2; 1.0e-7 1.0e-7 10.0 0.65686096 1.112720269474884e-6 153576.0; 1.0e-7 1.0e-7 25.0 0.7492210396 1.1151420272518578e-6 154324.8; 1.0e-7 1.0e-7 50.0 0.6684278204 1.1083522166780087e-6 154656.0; 1.0e-7 1.0e-7 100.0 0.6747090196000001 1.108006703839939e-6 154137.59999999998; 1.0e-7 1.0e-7 200.0 0.7606153996 1.1027263306236048e-6 155577.59999999998; 1.0e-7 1.0e-7 300.0 0.8000994996000002 1.107795433067466e-6 154152.0; 1.0e-7 1.0e-7 400.0 0.8600881198 1.0946299569630452e-6 155016.0; 1.0e-7 1.0e-7 500.0 0.8818279800000001 1.1171296792286267e-6 153864.0; 1.0e-7 1.0e-7 1000.0 1.09660146 1.0981056623915386e-6 154800.0; 1.0e-7 1.0e-7 1500.0 1.2895725598 1.1133291359844354e-6 155088.0; 1.0e-7 1.0e-7 2000.0 1.5808850797999998 1.1083096327284215e-6 154339.2; 1.0e-7 1.0e-7 2500.0 1.7350152602000002 1.0996656232828543e-6 155534.40000000002; 1.0e-7 1.0e-7 3000.0 2.0038747798 1.100191703501389e-6 155044.8]

B_RK4_DE_ad_m7=[1.0e-7 1.0e-7 3.0 0.7654632202 1.3266709193054647e-6 161640.0; 1.0e-7 1.0e-7 4.0 0.7670076398 1.205130807666481e-6 166118.40000000002; 1.0e-7 1.0e-7 5.0 0.7722214002 1.202850375217388e-6 169142.40000000002; 1.0e-7 1.0e-7 6.0 0.7129445998 1.1770044318015687e-6 168998.40000000002; 1.0e-7 1.0e-7 7.0 0.7919431602 1.1868417157833179e-6 169560.0; 1.0e-7 1.0e-7 8.0 0.7303591008 1.1930356234876818e-6 169977.59999999998; 1.0e-7 1.0e-7 9.0 1.0580429604 1.191961987892985e-6 167947.2; 1.0e-7 1.0e-7 10.0 0.8148082196 1.1879024862811398e-6 169344.0; 1.0e-7 1.0e-7 25.0 0.8549851599999998 1.1865468048627462e-6 169401.59999999998; 1.0e-7 1.0e-7 50.0 0.76139398 1.1920747411976947e-6 169444.8; 1.0e-7 1.0e-7 100.0 0.8022676996 1.1840433534557324e-6 170870.40000000002; 1.0e-7 1.0e-7 200.0 0.8576324200000001 1.188604613212792e-6 170697.59999999998; 1.0e-7 1.0e-7 300.0 0.8606759404 1.1854330636753956e-6 170078.40000000002; 1.0e-7 1.0e-7 400.0 0.9167928798 1.1840632793171453e-6 169905.59999999998; 1.0e-7 1.0e-7 500.0 0.93292348 1.1886703918572723e-6 169531.2; 1.0e-7 1.0e-7 1000.0 1.1552085002 1.1874958385728272e-6 169286.40000000002; 1.0e-7 1.0e-7 1500.0 1.3708232398 1.192862354906349e-6 169228.8; 1.0e-7 1.0e-7 2000.0 1.5898113600000001 1.1908955469365201e-6 169185.59999999998; 1.0e-7 1.0e-7 2500.0 1.9313150602 1.18955458720659e-6 169531.2; 1.0e-7 1.0e-7 3000.0 2.0357141594 1.1871521265424033e-6 170856.0]

A_RK4_DE_ad_m6=[1.0e-6 1.0e-6 3.0 0.3914056402 1.1012708196240696e-5 87969.6; 1.0e-6 1.0e-6 4.0 0.3824554402 9.5575168601401e-6 90648.0; 1.0e-6 1.0e-6 5.0 0.3913520202 9.108880975281281e-6 92548.79999999999; 1.0e-6 1.0e-6 6.0 0.38696105999999997 8.885755140619754e-6 91886.4; 1.0e-6 1.0e-6 7.0 0.3798095998 8.55677097673478e-6 92548.79999999999; 1.0e-6 1.0e-6 8.0 0.38491734 8.686549353467055e-6 92808.0; 1.0e-6 1.0e-6 9.0 0.40270316 8.701737696607564e-6 92865.6; 1.0e-6 1.0e-6 10.0 0.3908863 8.790204498467926e-6 92073.6; 1.0e-6 1.0e-6 25.0 0.40393334000000003 8.814113240404552e-6 92433.6; 1.0e-6 1.0e-6 50.0 0.48157956 8.744584550011041e-6 92347.20000000001; 1.0e-6 1.0e-6 100.0 0.5427241398 8.606433900462525e-6 93024.0; 1.0e-6 1.0e-6 200.0 0.5293874198 8.645197519582844e-6 92966.4; 1.0e-6 1.0e-6 300.0 0.5908004593999999 8.751079978815042e-6 92404.79999999999; 1.0e-6 1.0e-6 400.0 0.9412951996000001 8.603116819992848e-6 92779.20000000001; 1.0e-6 1.0e-6 500.0 0.7695598196 8.611095974325918e-6 92289.6; 1.0e-6 1.0e-6 1000.0 0.8364510598000001 8.610475587353096e-6 93139.20000000001; 1.0e-6 1.0e-6 1500.0 1.0502910602 8.534442550028404e-6 92894.4; 1.0e-6 1.0e-6 2000.0 1.3293679597999999 8.583352682757355e-6 92995.20000000001; 1.0e-6 1.0e-6 2500.0 1.4624705796000002 8.77274791434814e-6 92923.20000000001; 1.0e-6 1.0e-6 3000.0 1.6781204203999998 8.620832859272814e-6 92779.20000000001]

B_RK4_DE_ad_m6=[1.0e-6 1.0e-6 3.0 0.41917588000000006 1.040005110029359e-5 100728.0; 1.0e-6 1.0e-6 4.0 0.4304889998 9.582193735823036e-6 103924.79999999999; 1.0e-6 1.0e-6 5.0 0.4414099998 9.332463504081461e-6 103622.40000000001; 1.0e-6 1.0e-6 6.0 0.5379343798 9.309930887679349e-6 105062.40000000001; 1.0e-6 1.0e-6 7.0 0.50015012 9.249207885278508e-6 104601.59999999999; 1.0e-6 1.0e-6 8.0 0.48829571979999997 9.260633160672105e-6 105278.40000000001; 1.0e-6 1.0e-6 9.0 0.5141686402000001 9.281429944937001e-6 104860.79999999999; 1.0e-6 1.0e-6 10.0 0.46372807979999997 9.24550489170476e-6 104558.40000000001; 1.0e-6 1.0e-6 25.0 0.4521661 9.270204526828285e-6 105609.59999999999; 1.0e-6 1.0e-6 50.0 0.46973592020000005 9.246944890519932e-6 105465.59999999999; 1.0e-6 1.0e-6 100.0 0.5461482202000001 9.204557610176815e-6 104558.40000000001; 1.0e-6 1.0e-6 200.0 0.58359942 9.256190003639635e-6 105969.59999999999; 1.0e-6 1.0e-6 300.0 0.67215662 9.250770073074358e-6 105566.40000000001; 1.0e-6 1.0e-6 400.0 0.6351241804000001 9.257242512898954e-6 104457.59999999999; 1.0e-6 1.0e-6 500.0 0.71255704 9.2554063697951e-6 104947.20000000001; 1.0e-6 1.0e-6 1000.0 0.9621857803999999 9.270127115463977e-6 104400.0; 1.0e-6 1.0e-6 1500.0 1.1093156998000002 9.264572116098938e-6 104644.79999999999; 1.0e-6 1.0e-6 2000.0 1.4033159198 9.248335722935331e-6 105609.59999999999; 1.0e-6 1.0e-6 2500.0 1.5033980398 9.284414181062744e-6 105724.79999999999; 1.0e-6 1.0e-6 3000.0 1.7970735202 9.252409393528691e-6 104788.79999999999]

Plots.plot(A_RK4_DE_ad_m8[:,3],hcat(A_RK4_DE_ad_m8[:,5],B_RK4_DE_ad_m8[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["A RK4 1e-8" "B RK4 1e-8"],markershapes= [:circle],xlabel = "sampling segment number",ylabel = "eigenvalue error")

Plots.plot(A_RK4_DE_ad_m6[:,3],hcat(A_RK4_DE_ad_m6[:,5],B_RK4_DE_ad_m6[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["A RK4 1e-6" "B RK4 1e-6"],markershapes= [:circle],xlabel = "sampling segment number",ylabel = "eigenvalue error")

Plots.plot(A_RK4_DE_ad_m6[:,3],hcat(A_RK4_DE_ad_m6[:,5],B_RK4_DE_ad_m6[:,5],A_RK4_DE_ad_m8[:,5],B_RK4_DE_ad_m8[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["A RK4 1e-6" "B RK4 1e-6" "A RK4 1e-8" "B RK4 1e-8"],markershapes= [:circle],xlabel = "sampling segment number",ylabel = "eigenvalue error")

Plots.plot(A_RK4_DE_ad_m6[:,3],hcat(A_RK4_DE_ad_m6[:,5],A_RK4_DE_ad_m7[:,5],A_RK4_DE_ad_m8[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["A RK4 1e-6" "A RK4 1e-7" "A RK4 1e-8"],markershapes= [:circle],xlabel = "sampling segment number",ylabel = "eigenvalue error")

#reltol & abstol together
A_BS3_DE_ad_tog=[0.1 1.0e-7 100.0 0.1070325194 0.4841976642705301 5112.0; 0.01 1.0e-7 100.0 0.11140026 0.14743377510645497 7286.400000000001; 0.001 1.0e-7 100.0 0.1306552602 0.033334303628016544 11340.0; 0.0001 1.0e-7 100.0 0.1590048202 0.005107889339486951 19389.6; 1.0e-5 1.0e-7 100.0 0.2252371802 0.0007129582584277991 35352.0; 1.0e-6 1.0e-7 100.0 0.33907864 9.645129064646591e-5 66643.20000000001; 1.0e-7 1.0e-7 100.0 0.58878056 1.3136414142837385e-5 128347.20000000001; 1.0e-8 1.0e-7 100.0 1.1116801398 1.7924652291541132e-6 250473.59999999998; 1.0e-9 1.0e-7 100.0 2.1089492598 2.459017352072497e-7 487785.60000000003; 1.0e-10 1.0e-7 100.0 4.3797708602 3.382048889813197e-8 955159.2000000001; 1.0e-11 1.0e-7 100.0 8.1761995196 4.678134097966194e-9 1.8841247999999998e6; 1.0e-12 1.0e-7 100.0 16.4648836 6.4702251469046e-10 3.7029528000000003e6]

A_RK4_DE_ad_tog=[0.1 1.0e-7 100.0 0.13286368019999997 0.11194311856012124 10411.2; 0.01 1.0e-7 100.0 0.14299854 0.02216573531768844 14745.599999999999; 0.001 1.0e-7 100.0 0.17227363960000003 0.0035905867117381313 22507.199999999997; 0.0001 1.0e-7 100.0 0.21630598019999997 0.0005125126472560284 34545.600000000006; 1.0e-5 1.0e-7 100.0 0.3301012996 6.764950639098015e-5 58924.799999999996; 1.0e-6 1.0e-7 100.0 0.45470125980000003 8.849106415397762e-6 91684.79999999999; 1.0e-7 1.0e-7 100.0 0.7184975 1.0952819758178724e-6 154713.59999999998; 1.0e-8 1.0e-7 100.0 1.0615118998 1.423803186587993e-7 251078.40000000002; 1.0e-9 1.0e-7 100.0 1.7514673002000003 1.814231106363586e-8 409190.39999999997; 1.0e-10 1.0e-7 100.0 2.7128933202 2.314637604266228e-9 671904.0; 1.0e-11 1.0e-7 100.0 4.75764722 3.008489589882364e-10 1.1181744e6; 1.0e-12 1.0e-7 100.0 7.613931619800001 3.929679055472035e-11 1.845432e6]

A_T5_DE_ad_tog=[0.1 1.0e-7 100.0 0.11728488020000001 0.07867771106184122 7286.400000000001; 0.01 1.0e-7 100.0 0.12043235940000001 0.011836258943494197 9100.8; 0.001 1.0e-7 100.0 0.1301070198 0.0018982505935387286 11908.8; 0.0001 1.0e-7 100.0 0.14835440020000001 0.0002999430762640492 16416.0; 1.0e-5 1.0e-7 100.0 0.16431718039999998 2.8542295329834144e-5 22204.800000000003; 1.0e-6 1.0e-7 100.0 0.2077842802 2.646472718751679e-6 31838.399999999998; 1.0e-7 1.0e-7 100.0 0.2790689598 2.420042227430653e-7 50414.399999999994; 1.0e-8 1.0e-7 100.0 0.3744239 2.383129211152763e-8 74246.4; 1.0e-9 1.0e-7 100.0 0.5430784396 2.575202921506168e-9 113616.0; 1.0e-10 1.0e-7 100.0 0.7222152599999999 2.9371658850556877e-10 168422.40000000002; 1.0e-11 1.0e-7 100.0 1.1001593202000002 3.4863617219561184e-11 253627.19999999998; 1.0e-12 1.0e-7 100.0 1.4918969002 4.414630335200242e-12 373636.80000000005]

Plots.plot(A_BS3_DE_ad_tog[:,1],hcat(A_BS3_DE_ad_tog[:,5],A_RK4_DE_ad_tog[:,5],A_T5_DE_ad_tog[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3" "RK4" "T5"],markershapes= [:circle],xlabel = "reltol=abstol",ylabel = "eigenvalue error")

Plots.plot(A_BS3_DE_ad_tog[:,1],hcat(A_BS3_DE_ad_tog[:,5],A_RK4_DE_ad_tog[:,5],A_T5_DE_ad_tog[:,5]),title="Largest eigenvalue error A",label=["BS3" "RK4" "T5"],markershapes= [:circle],xlabel = "reltol=abstol",ylabel = "eigenvalue error")

Plots.plot(A_BS3_DE_ad_tog[:,1],hcat(10^(-3)*A_BS3_DE_ad_tog[:,5],10^(-2)*A_RK4_DE_ad_tog[:,5],10^(-1)*A_T5_DE_ad_tog[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3" "RK4" "T5"],markershapes= [:circle],xlabel = "reltol=abstol",ylabel = "eigenvalue error")


println(A_RK4_DE_ad)


open("temp.txt", "w") do io
           writedlm(io, A_T5_J_dense_BF128)
       end

Plots.savefig("1.pdf")
