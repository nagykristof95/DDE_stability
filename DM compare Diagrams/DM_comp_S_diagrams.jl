using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using Plots

valref=-0.6450472841316337 + 0.8845000929140697*im

RK4_fix=[12.0 12.0 100.0 0.5701608505 5.611910075416179e-6 4764.0; 12.0 12.0 200.0 0.77254 4.041268951334547e-7 9564.0; 12.0 12.0 300.0 1.0524615495 7.92359553561787e-8 14364.0; 12.0 12.0 400.0 1.4469343499999998 2.49782207664564e-8 19212.0; 12.0 12.0 500.0 1.58016135 1.0208942140164447e-8 23964.0; 12.0 12.0 600.0 1.91384115 4.916497202538448e-9 28812.0; 12.0 12.0 700.0 2.46620025 2.651486051443305e-9 33612.0; 12.0 12.0 800.0 2.4257658005 1.5533955907523152e-9 38412.0; 12.0 12.0 900.0 2.5821023 9.69689508002345e-10 43164.0; 12.0 12.0 1000.0 2.9413020999999997 6.363962741836754e-10 47964.0; 12.0 12.0 1100.0 3.5042309495 4.349363406672053e-10 52764.0; 12.0 12.0 1200.0 3.6218474499999997 3.072360143262929e-10 57564.0; 12.0 12.0 1300.0 4.4128320495 2.2336314971558443e-10 62364.0; 12.0 12.0 1400.0 4.26291805 1.663307204130392e-10 67164.0; 12.0 12.0 1500.0 5.375166549999999 1.2652043856067543e-10 71964.0]

RK4_ad=[12.0 12.0 3.1622776601683795e-5 0.2003487504999997 0.00018722759496548033 3498.0; 12.0 12.0 1.0000000000000003e-5 0.3429037995 6.781306680751595e-5 4413.0; 12.0 12.0 3.162277660168379e-6 0.4116517495 2.450242158797252e-5 5814.0; 12.0 12.0 1.0000000000000004e-6 0.467455101 8.899470109432265e-6 7407.0; 12.0 12.0 3.162277660168379e-7 0.5227270495 3.0851828926254177e-6 9285.0; 12.0 12.0 1.0000000000000004e-7 0.6479187 1.1186319480769324e-6 11907.0; 12.0 12.0 3.162277660168379e-8 0.79545815 4.0056456488979225e-7 15348.0; 12.0 12.0 1.0000000000000004e-8 1.0940742005 1.3990131862391665e-7 19818.0; 12.0 12.0 3.1622776601683795e-9 1.3287096004999999 4.9923513089516425e-8 25515.0; 12.0 12.0 1.0000000000000009e-9 1.8612521000000002 1.7820628147263774e-8 32691.0; 12.0 12.0 3.1622776601683795e-10 2.23300605 6.444857454919429e-9 41889.0; 12.0 12.0 1.0000000000000006e-10 2.753883549 2.3561730444315482e-9 53673.0; 12.0 12.0 3.1622776601683794e-11 3.53709545 8.400829697948713e-10 68997.0; 12.0 12.0 1.0000000000000006e-11 4.6258949000000005 2.9689594633969563e-10 89178.0; 12.0 12.0 3.1622776601683794e-12 5.94939045 1.0625009416170822e-10 115560.0]

T5_ad=[12.0 12.0 3.1622776601683795e-5 1.722410060399999 9.725532441370628e-5 1544.4; 12.0 12.0 1.0000000000000003e-5 0.3949656202 2.854147295234508e-5 1788.0; 12.0 12.0 3.162277660168379e-6 0.4131226 9.076327288558391e-6 2167.2; 12.0 12.0 1.0000000000000004e-6 0.45730770000000004 2.635577646276202e-6 2538.0; 12.0 12.0 3.162277660168379e-7 0.5075906993999999 7.98167497415004e-7 3080.4; 12.0 12.0 1.0000000000000004e-7 0.5869997600000001 2.414948393988591e-7 3741.6; 12.0 12.0 3.162277660168379e-8 0.6608166404 7.543628971680325e-8 4570.8; 12.0 12.0 1.0000000000000004e-8 0.8092471 2.3814353606039958e-8 5965.2; 12.0 12.0 3.1622776601683795e-9 0.9506194398 7.741855203133844e-9 7164.0; 12.0 12.0 1.0000000000000009e-9 1.10205238 2.5727901247567735e-9 8695.2; 12.0 12.0 3.1622776601683795e-10 1.3623995397999997 8.60490540825145e-10 10768.8; 12.0 12.0 1.0000000000000006e-10 1.6344414 2.9254798637572813e-10 13238.4; 12.0 12.0 3.1622776601683794e-11 2.0124462802 1.0094331672148904e-10 16389.6; 12.0 12.0 1.0000000000000006e-11 2.19221926 3.496228590658986e-11 19764.0; 12.0 12.0 3.1622776601683794e-12 2.7173286400000003 1.2157095388851085e-11 24256.8]

RK4_G_ad=[100.0 1.0e-5 3.1622776601683795e-5 0.622081 0.00017391247700356991 3510.0; 100.0 1.0e-5 1.0000000000000003e-5 0.718598251 6.527274980912781e-5 4398.0; 100.0 1.0e-5 3.162277660168379e-6 0.7462067 2.469276788267424e-5 6000.0; 100.0 1.0e-5 1.0000000000000004e-6 0.749655151 8.646228315755335e-6 7443.0; 100.0 1.0e-5 3.162277660168379e-7 0.8607867005000001 2.8369713081074205e-6 9375.0; 100.0 1.0e-5 1.0000000000000004e-7 0.945126299 1.180132382513211e-6 11943.0; 100.0 1.0e-5 3.162277660168379e-8 1.0904915 4.0253229660977405e-7 15552.0; 100.0 1.0e-5 1.0000000000000004e-8 1.61145555 1.279394680646319e-7 20088.0; 100.0 1.0e-5 3.1622776601683795e-9 1.9056050505000002 4.4721209623527605e-8 25476.0; 100.0 1.0e-5 1.0000000000000009e-9 2.0385994 1.679471446377304e-8 32655.0; 100.0 1.0e-5 3.1622776601683795e-10 2.4039831 6.640554009472552e-9 42354.0; 100.0 1.0e-5 1.0000000000000006e-10 2.8413189 2.345825886917558e-9 53862.0; 100.0 1.0e-5 3.1622776601683794e-11 4.0038067 7.969139917132185e-10 70050.0; 100.0 1.0e-5 1.0000000000000006e-11 4.7116430995 2.8687188447948505e-10 88431.0; 100.0 1.0e-5 3.1622776601683794e-12 6.1218765 1.1380971794278199e-10 115023.0]

T5_G_ad=[100.0 1.0e-5 3.1622776601683795e-5 1.105307219999999 9.640348462357049e-5 1545.6; 100.0 1.0e-5 1.0000000000000003e-5 1.3342457196000002 2.877472552639936e-5 1792.8; 100.0 1.0e-5 3.162277660168379e-6 1.3849489199999998 9.17561169918838e-6 2163.6; 100.0 1.0e-5 1.0000000000000004e-6 1.2735504604 2.5823667323669474e-6 2548.8; 100.0 1.0e-5 3.162277660168379e-7 1.2985971403999998 1.13074758506319e-6 3085.2; 100.0 1.0e-5 1.0000000000000004e-7 1.3644414596 4.431201335833715e-7 3757.2; 100.0 1.0e-5 3.162277660168379e-8 1.3980902198 1.3443161588253424e-7 4593.6; 100.0 1.0e-5 1.0000000000000004e-8 1.5608347197999999 3.843673186020572e-8 5920.8; 100.0 1.0e-5 3.1622776601683795e-9 1.7018543003999997 9.9329777076434037e-9 7165.2; 100.0 1.0e-5 1.0000000000000009e-9 1.7257789398 4.545993759080599e-9 8732.4; 100.0 1.0e-5 3.1622776601683795e-10 1.9362011398000003 4.9787805301052556e-9 10770.0; 100.0 1.0e-5 1.0000000000000006e-10 2.1766102396 1.7680257175694256e-9 13304.4; 100.0 1.0e-5 3.1622776601683794e-11 2.6201940192000004 2.91932332150498e-10 16358.4; 100.0 1.0e-5 1.0000000000000006e-11 2.9262457798 7.299339574644356e-11 19854.0; 100.0 1.0e-5 3.1622776601683794e-12 3.1927341996 1.661472843422859e-11 24434.4]

RK4_C_ad=[100.0 1.0e-5 3.1622776601683795e-5 0.8314084 0.0003254335495273716 5430.0; 100.0 1.0e-5 1.0000000000000003e-5 1.4352561504999999 0.00011917901376041757 6690.0; 100.0 1.0e-5 3.162277660168379e-6 1.7795218505000001 4.448845636420501e-5 8820.0; 100.0 1.0e-5 1.0000000000000004e-6 2.155291301 1.6097276646223648e-5 11106.0; 100.0 1.0e-5 3.162277660168379e-7 2.7766238005 6.020904755068209e-6 14475.0; 100.0 1.0e-5 1.0000000000000004e-7 3.7059437500000003 2.16881718536129e-6 18525.0; 100.0 1.0e-5 3.162277660168379e-8 4.9869304505 7.926550927444194e-7 23595.0; 100.0 1.0e-5 1.0000000000000004e-8 5.8941536 2.8879723447632937e-7 29976.0; 100.0 1.0e-5 3.1622776601683795e-9 7.5346204 1.052617716676258e-7 38001.0; 100.0 1.0e-5 1.0000000000000009e-9 9.34531475 3.756318875257958e-8 49050.0; 100.0 1.0e-5 3.1622776601683795e-10 12.18196125 1.3997429002826038e-8 62808.0; 100.0 1.0e-5 1.0000000000000006e-10 16.590567751000002 5.036402030047819e-9 78354.0; 100.0 1.0e-5 3.1622776601683794e-11 20.783578 1.8776427398439445e-9 100620.0; 100.0 1.0e-5 1.0000000000000006e-11 24.41934405 6.81835561768875e-10 128367.0; 100.0 1.0e-5 3.1622776601683794e-12 30.6685251995 2.5019039837740573e-10 163476.0]

RK4_C_ad_norm=[100.0 1.0e-5 3.1622776601683795e-5 11.6327558 4.7400455992821485e-5 7052.0; 100.0 1.0e-5 1.0000000000000003e-5 1.925683533 1.7279997328956706e-5 8918.0; 100.0 1.0e-5 3.162277660168379e-6 2.4044441659999998 6.282563817593964e-6 11294.0; 100.0 1.0e-5 1.0000000000000004e-6 3.0018388666666667 2.317598388871734e-6 14418.0; 100.0 1.0e-5 3.162277660168379e-7 3.5984431666666663 8.370354862730213e-7 18106.0; 100.0 1.0e-5 1.0000000000000004e-7 4.652442100666666 3.028172633587245e-7 23076.0; 100.0 1.0e-5 3.162277660168379e-8 5.809987100333333 1.0867905287453863e-7 29414.0; 100.0 1.0e-5 1.0000000000000004e-8 7.441262266333332 3.966824312197939e-8 37504.0; 100.0 1.0e-5 3.1622776601683795e-9 9.895482332999999 1.4555692136797187e-8 47384.0; 100.0 1.0e-5 1.0000000000000009e-9 11.964774366666665 5.283625698207438e-9 60346.0; 100.0 1.0e-5 3.1622776601683795e-10 14.944468399666667 1.92013920017119e-9 77198.0; 100.0 1.0e-5 1.0000000000000006e-10 20.276550567 7.001250047461759e-10 99476.0; 100.0 1.0e-5 3.1622776601683794e-11 24.126524699666664 2.56264213678501e-10 126476.0; 100.0 1.0e-5 1.0000000000000006e-11 31.465646000333333 9.376757597343939e-11 161952.0; 100.0 1.0e-5 3.1622776601683794e-12 41.238248667 3.504462493441898e-11 208930.0]

T5_C_ad=[100.0 1.0e-5 3.1622776601683795e-5 1.279056540400002 0.00024722191855550374 1936.8; 100.0 1.0e-5 1.0000000000000003e-5 1.2453084202 7.689244789456937e-5 2468.4; 100.0 1.0e-5 3.162277660168379e-6 1.63913276 2.1084870288485964e-5 3120.0; 100.0 1.0e-5 1.0000000000000004e-6 1.8001416398 6.5736335053958294e-6 3842.4; 100.0 1.0e-5 3.162277660168379e-7 2.2793994603999996 1.8799504026660795e-6 5042.4; 100.0 1.0e-5 1.0000000000000004e-7 2.60137246 5.004380740246737e-7 5818.8; 100.0 1.0e-5 3.162277660168379e-8 3.26607176 1.5620220312534998e-7 7246.8; 100.0 1.0e-5 1.0000000000000004e-8 3.9956549202000007 4.773366762844681e-8 8977.2; 100.0 1.0e-5 3.1622776601683795e-9 4.6872603402 1.4175928751522023e-8 10958.4; 100.0 1.0e-5 1.0000000000000009e-9 5.614051959999999 4.574634488479689e-9 13116.0; 100.0 1.0e-5 3.1622776601683795e-10 6.9868081596 1.493352092015288e-9 16339.2; 100.0 1.0e-5 1.0000000000000006e-10 8.5644696398 5.024664503090551e-10 20055.6; 100.0 1.0e-5 3.1622776601683794e-11 10.4799870404 1.7350067149155702e-10 24514.8; 100.0 1.0e-5 1.0000000000000006e-11 12.769900019600001 5.717970689486481e-11 30189.6; 100.0 1.0e-5 3.1622776601683794e-12 15.827999560200002 1.9932706035183753e-11 36968.4]

T5_C_ad_norm=[100.0 1.0e-5 3.1622776601683795e-5 8.5067759495 2.601976198578612e-5 2616.0; 100.0 1.0e-5 1.0000000000000003e-5 0.6838879495000001 7.904326019416686e-6 3117.0; 100.0 1.0e-5 3.162277660168379e-6 0.8271648 2.4284850495941098e-6 3756.0; 100.0 1.0e-5 1.0000000000000004e-6 0.9659769 7.058107429990469e-7 4596.0; 100.0 1.0e-5 3.162277660168379e-7 1.1434881495 2.1373753743220803e-7 5424.0; 100.0 1.0e-5 1.0000000000000004e-7 1.386314 6.110888428059659e-8 6552.0; 100.0 1.0e-5 3.162277660168379e-8 1.6853056505 1.878257922617229e-8 8124.0; 100.0 1.0e-5 1.0000000000000004e-8 2.0423584995 5.497073277202368e-9 9942.0; 100.0 1.0e-5 3.1622776601683795e-9 2.4973586 1.8689220737186543e-9 12201.0; 100.0 1.0e-5 1.0000000000000009e-9 2.9810454005 6.044517696177685e-10 14748.0; 100.0 1.0e-5 3.1622776601683795e-10 3.9137082994999997 2.011805123563884e-10 18303.0; 100.0 1.0e-5 1.0000000000000006e-10 4.6711159005 6.748034912622587e-11 22446.0; 100.0 1.0e-5 3.1622776601683794e-11 5.54599375 2.340553277900067e-11 27168.0; 100.0 1.0e-5 1.0000000000000006e-11 6.9406785499999994 8.126628614033054e-12 32907.0; 100.0 1.0e-5 3.1622776601683794e-12 8.136130099999999 2.9436644220802246e-12 39861.0]

#Plots.plot(hcat(RK4_fix[:,6],T5_ad[:,6],T5_G_ad[:,6],T5_C_ad[:,6]),hcat(RK4_fix[:,5],T5_ad[:,5],T5_G_ad[:,5],T5_C_ad[:,5]),xscale=:log10,yscale=:log10,label=["RK4 naive" "T5 ad" "T5 G ad" "T5 C ad"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

Plots.plot(hcat(RK4_fix[:,6],RK4_ad[:,6],RK4_G_ad[:,6],RK4_C_ad_norm[:,6],T5_G_ad[:,6]),hcat(RK4_fix[:,5],RK4_ad[:,5],RK4_G_ad[:,5],RK4_C_ad_norm[:,5],T5_G_ad[:,5]),xscale=:log10,yscale=:log10,label=["RK4 naive" "RK4 ad" "RK4 G ad" "RK4 C ad" "T5 G ad"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

Plots.plot(hcat(RK4_fix[:,4],RK4_ad[:,4],RK4_G_ad[:,4],RK4_C_ad[:,4],T5_G_ad[:,4]),hcat(RK4_fix[:,5],RK4_ad[:,5],RK4_G_ad[:,5],RK4_C_ad[:,5],T5_G_ad[:,5]),xscale=:log10,yscale=:log10,label=["RK4 naive" "RK4 ad" "RK4 G ad" "RK4 C ad" "T5 G ad"],markershapes= [:circle],xlabel = "CPU time", ylabel = "eigenvalue error")

Plots.savefig("1.pdf")
println(A_RK4_T)