#RK linear interpolation
A_RK4l=[20.0 20.0 10.0 0.110958501 0.1497034803615127 720.0; 20.0 20.0 30.0 0.3356477 0.006400808495084624 2320.0; 20.0 20.0 50.0 0.582289399 0.002076203814808402 3920.0; 20.0 20.0 70.0 0.5916901 0.001028045212959581 5520.0; 20.0 20.0 100.0 0.6620699 0.0004950195208678912 7920.0; 20.0 20.0 150.0 0.9585936 0.00021759303696141502 11920.0; 20.0 20.0 250.0 1.519647 7.774851426997128e-5 19920.0; 20.0 20.0 500.0 2.9956476 1.9334654869249643e-5 39920.0; 20.0 20.0 1500.0 8.7806613 2.1313702214778167e-6 119920.0; 20.0 20.0 2500.0 16.9850497 7.591613701156019e-7 199920.0; 20.0 20.0 3500.0 23.2621645 3.813374848438197e-7 279920.0; 20.0 20.0 4500.0 28.937041901 2.2589781190886927e-7 359920.0]

A_RK3l=[20.0 20.0 10.0 0.053736701 0.3683933872456685 540.0; 20.0 20.0 30.0 0.181240499 0.016562754823629035 1740.0; 20.0 20.0 50.0 0.247550801 0.004248175087838617 2940.0; 20.0 20.0 70.0 0.3432164 0.0018097804529526476 4140.0; 20.0 20.0 100.0 0.485631 0.0007596057373673336 5940.0; 20.0 20.0 150.0 0.7470994 0.00029494643860624565 8940.0; 20.0 20.0 250.0 1.2393196 9.423782827990002e-5 14940.0; 20.0 20.0 500.0 2.5791675 2.1371675525646325e-5 29940.0; 20.0 20.0 1500.0 7.4316077 2.206112964694462e-6 89940.0; 20.0 20.0 2500.0 13.179387 7.752504929500197e-7 149940.0; 20.0 20.0 3500.0 18.252786 3.87176075063149e-7 209940.0; 20.0 20.0 4500.0 22.397835401 2.2862701211884536e-7 269940.0]

A_RK2l=[20.0 20.0 10.0 0.039490001 2.2452505647572583 360.0; 20.0 20.0 30.0 0.118820301 0.13572192567752794 1160.0; 20.0 20.0 50.0 0.1772955 0.046893136625024014 1960.0; 20.0 20.0 70.0 0.2604265 0.023581756582362504 2760.0; 20.0 20.0 100.0 0.3534018 0.011439023962923298 3960.0; 20.0 20.0 150.0 0.5396668 0.005045550497153207 5960.0; 20.0 20.0 250.0 0.89752 0.001805425949276801 9960.0; 20.0 20.0 500.0 1.730379801 0.0004491399397363598 19960.0; 20.0 20.0 1500.0 4.954689601 4.950739160153118e-5 59960.0; 20.0 20.0 2500.0 9.129880201 1.763297868149476e-5 99960.0; 20.0 20.0 3500.0 13.455762999 8.857099016156519e-6 139960.0; 20.0 20.0 4500.0 16.244942099 5.2467265698587205e-6 179960.0]

A_RK1l=[20.0 20.0 10.0 0.0316417 23.928003750514463 180.0; 20.0 20.0 30.0 0.0964844 2.510711840807373 580.0; 20.0 20.0 50.0 0.102833499 1.071663140012098 980.0; 20.0 20.0 70.0 0.1607651 0.6651274578749566 1380.0; 20.0 20.0 100.0 0.2477339 0.4202979794854556 1980.0; 20.0 20.0 150.0 0.3792049 0.2590767946252088 2980.0; 20.0 20.0 250.0 0.5486924 0.14580658731159907 4980.0; 20.0 20.0 500.0 1.091435501 0.06893763671352905 9980.0; 20.0 20.0 1500.0 3.0826311 0.021192553807768008 29980.0; 20.0 20.0 2500.0 5.576298501 0.011960903024989566 49980.0; 20.0 20.0 3500.0 8.0687035 0.00803599357678825 69980.0; 20.0 20.0 4500.0 10.6032467 0.005863600761956331 89980.0]

Plots.plot(A_RK4l[:,3],hcat(A_RK4l[:,5],A_RK3l[:,5],A_RK2l[:,5],A_RK1l[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK4 l" "RK3 l" "RK2 l" "RK1 l"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4l[:,6],A_RK3l[:,6],A_RK2l[:,6],A_RK1l[:,6]),hcat(A_RK4l[:,5],A_RK3l[:,5],A_RK2l[:,5],A_RK1l[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["RK4 l" "RK3 l" "RK2 l" "RK1 l"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")
println(tangent((A_RK1c[5,3],A_RK1c[5,5]),(A_RK1c[end,3],A_RK1c[end,5])))
#RK cubic interpolation
A_RK4c=[20.0 20.0 10.0 0.107672101 0.10027183113750548 720.0; 20.0 20.0 30.0 0.2857937 0.0010083472729172796 2320.0; 20.0 20.0 50.0 0.392074799 0.00012984816807007883 3920.0; 20.0 20.0 70.0 0.5605815 3.479249487021417e-5 5520.0; 20.0 20.0 100.0 0.777422801 8.900512725719147e-6 7920.0; 20.0 20.0 150.0 1.2031029 1.9825493429090575e-6 11920.0; 20.0 20.0 250.0 1.957201301 3.2466211793739866e-7 19920.0; 20.0 20.0 500.0 3.912417201 3.210682708399265e-8 39920.0; 20.0 20.0 1500.0 11.3524038 1.017201095891938e-9 119920.0; 20.0 20.0 2500.0 20.813408499 2.1309844056891184e-10 199920.0; 20.0 20.0 3500.0 26.0576498 7.646918512531606e-11 279920.0; 20.0 20.0 4500.0 30.705724801 3.553083746213508e-11 359920.0]

A_RK3c=[20.0 20.0 10.0 0.0923469 0.3582315293138603 540.0; 20.0 20.0 30.0 0.189434599 0.012347555961036288 1740.0; 20.0 20.0 50.0 0.285629601 0.0025693101202235356 2940.0; 20.0 20.0 70.0 0.3915312 0.0009200506163059275 4140.0; 20.0 20.0 100.0 0.590274401 0.0003113207193899941 5940.0; 20.0 20.0 150.0 0.850716199 9.125723071429004e-5 8940.0; 20.0 20.0 250.0 1.456458801 1.9541161719926182e-5 14940.0; 20.0 20.0 500.0 3.036731499 2.4266759080835723e-6 29940.0; 20.0 20.0 1500.0 7.8637124 8.944820804867315e-8 89940.0; 20.0 20.0 2500.0 14.2128975 1.927435291233604e-8 149940.0; 20.0 20.0 3500.0 19.7742631 6.997594806964571e-9 209940.0; 20.0 20.0 4500.0 24.753090099 3.2717826540165166e-9 269940.0]

A_RK2c=[20.0 20.0 10.0 0.083179401 2.2998491456498034 360.0; 20.0 20.0 30.0 0.137600099 0.1410948706055126 1160.0; 20.0 20.0 50.0 0.2525713 0.04874112549664422 1960.0; 20.0 20.0 70.0 0.301890199 0.024504337154106038 2760.0; 20.0 20.0 100.0 0.454546301 0.011883250286069211 3960.0; 20.0 20.0 150.0 0.6593372 0.005240152341652274 5960.0; 20.0 20.0 250.0 1.082129601 0.0018746321946728549 9960.0; 20.0 20.0 500.0 2.1390951 0.0004662707441436834 19960.0; 20.0 20.0 1500.0 5.9189288 5.1389067882331584e-5 59960.0; 20.0 20.0 2500.0 10.8343248 1.8302695070258757e-5 99960.0; 20.0 20.0 3500.0 14.282001301 9.193397188343469e-6 139960.0; 20.0 20.0 4500.0 19.2072671 5.4459078882694925e-6 179960.0]

A_RK1c=[20.0 20.0 10.0 0.0460801 23.930840471660037 180.0; 20.0 20.0 30.0 0.1083295 2.510711840807372 580.0; 20.0 20.0 50.0 0.161546299 1.0716631400121 980.0; 20.0 20.0 70.0 0.2215056 0.6651274578749609 1380.0; 20.0 20.0 100.0 0.299770899 0.42029797948545805 1980.0; 20.0 20.0 150.0 0.476613 0.25907679462520894 2980.0; 20.0 20.0 250.0 0.839743199 0.1458065873115988 4980.0; 20.0 20.0 500.0 1.4841139 0.06893763671353062 9980.0; 20.0 20.0 1500.0 4.339357401 0.02119255380776824 29980.0; 20.0 20.0 2500.0 6.6523335 0.011960903024988865 49980.0; 20.0 20.0 3500.0 9.806183401 0.008035993576786743 69980.0; 20.0 20.0 4500.0 12.4485154 0.0058636007619537744 89980.0]

Plots.plot(A_RK4c[:,3],hcat(A_RK4c[:,5],A_RK3c[:,5],A_RK2c[:,5],A_RK1c[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK4 c" "RK3 c" "RK2 c" "RK1 c"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4c[:,6],A_RK3c[:,6],A_RK2c[:,6],A_RK1c[:,6]),hcat(A_RK4c[:,5],A_RK3c[:,5],A_RK2c[:,5],A_RK1c[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["RK4 c" "RK3 c" "RK2 c" "RK1 c"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#LMS
A_LM4=[20.0 20.0 10.0 0.090285 42.62794797908344 120.0; 20.0 20.0 30.0 0.1016461 0.03192370963759419 520.0; 20.0 20.0 50.0 0.183637 0.0045834818955555814 920.0; 20.0 20.0 70.0 0.213839099 0.0012425752088286446 1320.0; 20.0 20.0 100.0 0.289465201 0.0003081544226185536 1920.0; 20.0 20.0 150.0 0.536366801 6.287967608989551e-5 2920.0; 20.0 20.0 250.0 0.696223799 8.529039941177841e-6 4920.0; 20.0 20.0 500.0 1.324534399 5.853191882494499e-7 9920.0; 20.0 20.0 1500.0 4.0759072 9.96718411948922e-9 29920.0; 20.0 20.0 2500.0 6.812457999 1.6769354568609021e-9 49920.0; 20.0 20.0 3500.0 10.1941425 5.395856053736599e-10 69920.0; 20.0 20.0 4500.0 12.689187199 2.351281318108765e-10 89920.0]

A_LM3=[20.0 20.0 10.0 0.0690661 25.717956436124034 140.0; 20.0 20.0 30.0 0.1044422 0.10200436990228533 540.0; 20.0 20.0 50.0 0.1516739 0.02292951616714459 940.0; 20.0 20.0 70.0 0.231912201 0.008468517695800773 1340.0; 20.0 20.0 100.0 0.2856687 0.00292818103256263 1940.0; 20.0 20.0 150.0 0.4306129 0.0008720547567996693 2940.0; 20.0 20.0 250.0 0.680435201 0.0001890103594045526 4940.0; 20.0 20.0 500.0 1.324001001 2.367889072622592e-5 9940.0; 20.0 20.0 1500.0 3.836125799 8.778305984406302e-7 29940.0; 20.0 20.0 2500.0 6.496520101 1.8937017719232573e-7 49940.0; 20.0 20.0 3500.0 10.118902301 6.878448934681533e-8 69940.0; 20.0 20.0 4500.0 12.588358499 3.216932227368008e-8 89940.0]

A_LM2=[20.0 20.0 10.0 0.045983601 28.15692464914311 160.0; 20.0 20.0 30.0 0.109121001 0.38164879509569993 560.0; 20.0 20.0 50.0 0.177633399 0.12611955729218988 960.0; 20.0 20.0 70.0 0.227097699 0.06298414536416105 1360.0; 20.0 20.0 100.0 0.311142299 0.030521364622232247 1960.0; 20.0 20.0 150.0 0.4877486 0.013483295824409014 2960.0; 20.0 20.0 250.0 0.824409699 0.0048374328424360975 4960.0; 20.0 20.0 500.0 1.484748 0.0012067568042717519 9960.0; 20.0 20.0 1500.0 4.391016601 0.00013331115601954507 29960.0; 20.0 20.0 2500.0 7.712119401 4.750369172708328e-5 49960.0; 20.0 20.0 3500.0 10.4814444 2.386612002212117e-5 69960.0; 20.0 20.0 4500.0 12.788378399 1.4139297365273201e-5 89960.0]

A_LM1=[20.0 20.0 10.0 0.033665899 23.931431167183334 180.0; 20.0 20.0 30.0 0.09902 2.5107118408073705 580.0; 20.0 20.0 50.0 0.146583601 1.0716631400121037 980.0; 20.0 20.0 70.0 0.2041575 0.6651274578749554 1380.0; 20.0 20.0 100.0 0.2933531 0.4202979794854595 1980.0; 20.0 20.0 150.0 0.4433257 0.2590767946252077 2980.0; 20.0 20.0 250.0 0.7687429 0.1458065873116002 4980.0; 20.0 20.0 500.0 1.422427499 0.06893763671353029 9980.0; 20.0 20.0 1500.0 4.206607899 0.021192553807766093 29980.0; 20.0 20.0 2500.0 7.0501296 0.01196090302498808 49980.0; 20.0 20.0 3500.0 10.210700701 0.008035993576781532 69980.0; 20.0 20.0 4500.0 11.5382823 0.005863600761952998 89980.0]

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_LM4[:,6],A_LM3[:,6],A_LM2[:,6],A_LM1[:,6]),hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#DE fix step
A_EE_J=[20.0 20.0 10.0 0.1682286 23.928033108325362 200.0; 20.0 20.0 30.0 0.2488451 2.51071184080782 600.0; 20.0 20.0 50.0 0.3352442 1.071663140012322 1000.0; 20.0 20.0 70.0 0.4254232 0.6651274578750997 1400.0; 20.0 20.0 100.0 0.557453101 0.42029797948555003 2000.0; 20.0 20.0 150.0 0.7027151 0.259076794625268 3000.0; 20.0 20.0 250.0 1.0147153 0.14580658731163165 5020.0; 20.0 20.0 500.0 1.753148601 0.06893763671354501 10000.0; 20.0 20.0 1500.0 4.906641899 0.021192553807775374 30000.0; 20.0 20.0 2500.0 8.173084099 0.011960903024992637 50020.0; 20.0 20.0 3500.0 11.7002445 0.008035993576788058 70000.0; 20.0 20.0 4500.0 14.4436388 0.005863600761958353 90000.0]

A_BS3=[20.0 20.0 10.0 0.227372499 0.3599386238952525 560.0; 20.0 20.0 30.0 0.3856738 0.012321690941519194 1760.0; 20.0 20.0 50.0 0.4953148 0.0025605577371718286 2960.0; 20.0 20.0 70.0 0.657154201 0.0009164037677570488 4160.0; 20.0 20.0 100.0 0.8879734 0.00030995741459891957 5960.0; 20.0 20.0 150.0 1.1029465 9.082827440588292e-5 8960.0; 20.0 20.0 250.0 1.817085901 1.9444304178137192e-5 15020.0; 20.0 20.0 500.0 3.4766605 2.4141834323769285e-6 29960.0; 20.0 20.0 1500.0 10.205616001 8.89764750490993e-8 89960.0; 20.0 20.0 2500.0 17.188382999 1.917229337906743e-8 150020.0; 20.0 20.0 3500.0 20.9817115 6.960591480293147e-9 209960.0; 20.0 20.0 4500.0 26.5108581 3.2545862870892753e-9 269960.0]

A_RK4_J=[20.0 20.0 10.0 0.286694699 0.09387945332211947 740.0; 20.0 20.0 30.0 0.344631301 0.0009050062196188972 2340.0; 20.0 20.0 50.0 0.619240001 0.00011059795994635626 3940.0; 20.0 20.0 70.0 0.7449073 2.8068575822645574e-5 5540.0; 20.0 20.0 100.0 1.040191299 6.611910074304533e-6 7940.0; 20.0 20.0 150.0 1.446793801 1.2867810895347057e-6 11940.0; 20.0 20.0 250.0 2.2876473 1.6479287814937063e-7 20020.0; 20.0 20.0 500.0 4.469338999 1.0208941165034366e-8 39940.0; 20.0 20.0 1500.0 12.834756699 1.2651889988813484e-10 119940.0; 20.0 20.0 2500.0 21.1061047 1.7174347245183387e-11 200020.0; 20.0 20.0 3500.0 25.0667349 5.649488703709853e-12 279940.0; 20.0 20.0 4500.0 32.6382205 3.270199624183628e-12 359940.0]

A_T5=[20.0 20.0 10.0 0.2387284 0.00389381001228147 1100.0; 20.0 20.0 30.0 0.4861095 3.977248615980822e-6 3500.0; 20.0 20.0 50.0 0.876642901 2.0121126500707344e-7 5900.0; 20.0 20.0 70.0 1.0014253 3.070338204799763e-8 8300.0; 20.0 20.0 100.0 1.5512087 4.482938781695435e-9 11900.0; 20.0 20.0 150.0 2.0120807 5.369150836863848e-10 17900.0; 20.0 20.0 250.0 3.406847899 3.929026326686574e-11 30020.0; 20.0 20.0 500.0 6.732406399 1.719626177452509e-12 59900.0; 20.0 20.0 1500.0 19.884140501 1.2687401976755569e-12 179900.0; 20.0 20.0 2500.0 32.114985 9.758454432983036e-13 300020.0; 20.0 20.0 3500.0 37.7728349 1.4511483207308356e-12 419900.0; 20.0 20.0 4500.0 48.9236591 1.739803466959458e-12 539900.0]

Plots.plot(A_EE_J[:,3],hcat(A_EE_J[:,5],A_BS3[:,5],A_RK4_J[:,5],A_T5[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["EE DE" "BS3 DE" "RK4 DE" "T5 DE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_EE_J[:,6],A_BS3[:,6],A_RK4_J[:,6],A_T5[:,6]),hcat(A_EE_J[:,5],A_BS3[:,5],A_RK4_J[:,5],A_T5[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["EE DE" "BS3 DE" "RK4 DE" "T5 DE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#DE adaptive test
A_BS3_ad=[20.0 20.0 10.0 40.687823501 2.366160228974135e-12 446391.0; 20.0 20.0 30.0 44.515936799 1.2533754294881846e-11 451416.0; 20.0 20.0 50.0 42.8635245 4.593883814926827e-12 450024.0; 20.0 20.0 70.0 47.6982938 7.398986371252379e-13 450195.0; 20.0 20.0 100.0 39.5030339 6.05715758698511e-12 449901.0; 20.0 20.0 150.0 38.6794818 7.456731360040599e-12 451119.0; 20.0 20.0 250.0 40.603949 7.621149154863324e-12 451530.0; 20.0 20.0 500.0 39.8689897 6.7765729097108634e-12 453294.0; 20.0 20.0 1500.0 41.0715018 1.0593437783740724e-11 450003.0; 20.0 20.0 2500.0 42.7322665 1.1431155478688471e-11 451161.0; 20.0 20.0 3500.0 44.0465674 1.996882982849415e-12 448776.0; 20.0 20.0 4500.0 45.042638199 2.2316625577365206e-12 449544.0]

A_RK4_J_ad=[20.0 20.0 10.0 21.1109678 2.6368012731438937e-12 264684.0; 20.0 20.0 30.0 21.0594258 1.1540043976489936e-12 272508.0; 20.0 20.0 50.0 22.001728501 4.314045108975212e-13 276240.0; 20.0 20.0 70.0 20.6280804 4.3899236226710477e-13 270324.0; 20.0 20.0 100.0 20.4981847 2.5359579037618712e-12 267378.0; 20.0 20.0 150.0 21.000792001 2.571799110853689e-12 265806.0; 20.0 20.0 250.0 20.9434244 2.3786449511606593e-12 264690.0; 20.0 20.0 500.0 21.9141775 3.442294604339981e-13 273972.0; 20.0 20.0 1500.0 23.2321351 1.9723671979389325e-12 270204.0; 20.0 20.0 2500.0 24.818035 4.602113483709253e-13 274452.0; 20.0 20.0 3500.0 27.259607601 4.96512701808014e-13 270048.0; 20.0 20.0 4500.0 30.4893182 1.3773315740207392e-12 268632.0]

A_T5_ad=[20.0 20.0 10.0 4.424231899 5.256707398128031e-12 52722.0; 20.0 20.0 30.0 4.6042297 9.097273418055976e-14 53574.0; 20.0 20.0 50.0 4.870274801 1.9465157537712886e-13 53802.0; 20.0 20.0 70.0 4.8240084 4.847433210635068e-13 53598.0; 20.0 20.0 100.0 4.8896867 2.2584741724956378e-13 53058.0; 20.0 20.0 150.0 5.415791901 2.2964189691306596e-13 53322.0; 20.0 20.0 250.0 6.0961605 1.8619261021116505e-13 53862.0; 20.0 20.0 500.0 5.5976571 9.942541535518372e-14 53916.0; 20.0 20.0 1500.0 6.159136399 1.3851156871385394e-13 53184.0; 20.0 20.0 2500.0 7.557259699 2.1367902763799686e-13 54054.0; 20.0 20.0 3500.0 8.2645493 7.751742210816065e-14 53784.0; 20.0 20.0 4500.0 9.7583945 7.988824220291035e-14 53814.0]

Plots.plot(A_BS3_ad[:,3],hcat(A_BS3_ad[:,5],A_RK4_J_ad[:,5],A_T5_ad[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 DE ad" "RK4 DE ad" "T5 DE ad"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_BS3_ad[:,6],A_RK4_J_ad[:,6],A_T5_ad[:,6]),hcat(A_BS3_ad[:,5],A_RK4_J_ad[:,5],A_T5_ad[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["EE DE" "BS3 DE" "RK4 DE" "T5 DE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#PW+QuadGK
BS3_CPW=[20.0 20.0 0.001 2.0151306995000002 0.037773377403280524 1948.5; 20.0 20.0 1.0e-5 4.1653766505 0.0008113068414521994 5898.0; 20.0 20.0 1.0e-7 7.3316911000000005 1.447810566293381e-5 20104.5; 20.0 20.0 1.0e-9 11.8842463995 2.719265612524699e-7 72256.5]

RK4_CPW=[20.0 20.0 0.001 7.0446038 0.003970165768538991 4860.0; 20.0 20.0 1.0e-5 4.683130650500001 6.131673170182824e-5 11745.0; 20.0 20.0 1.0e-7 9.2837977495 9.951152049861347e-7 27558.0; 20.0 20.0 1.0e-9 19.686411698999997 1.6143597616055933e-8 67971.0]

T5_CPW=[20.0 20.0 0.001 4.51023265 0.0026157523945654834 1791.0; 20.0 20.0 1.0e-5 3.7111475499999997 3.197101139832525e-5 4083.0; 20.0 20.0 1.0e-7 7.50244765 2.533704123847033e-7 9684.0; 20.0 20.0 1.0e-9 15.50466885 2.6291894913214005e-9 20559.0]

Plots.plot(BS3_CPW[:,3],hcat(BS3_CPW[:,5],RK4_CPW[:,5],T5_CPW[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 CPW" "RK4 CPW" "T5 CPW"],markershapes= [:circle],xlabel = "solver tolerance",ylabel = "eigenvalue error")

Plots.plot(hcat(BS3_CPW[:,6],RK4_CPW[:,6],T5_CPW[:,6]),hcat(BS3_CPW[:,5],RK4_CPW[:,5],T5_CPW[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["BS3 CPW" "RK4 CPW" "T5 CPW"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#Continous
BS3_CI=[20.0 20.0 0.001 0.877408801 0.06121366128413316 2305.5; 20.0 20.0 1.0e-5 2.7443799999999996 0.002023075369671404 7576.5; 20.0 20.0 1.0e-7 7.2250929995 4.077386166895767e-5 24832.5; 20.0 20.0 1.0e-9 27.6596714995 8.722304074039303e-7 89665.5]

RK4_CI=[20.0 20.0 0.001 9.3711363995 0.003650374308771963 6696.0; 20.0 20.0 1.0e-5 4.7455687 0.00012958420917817952 15303.0; 20.0 20.0 1.0e-7 10.296638550499999 3.095723396612297e-6 34158.0; 20.0 20.0 1.0e-9 26.45013065 5.885437228135283e-8 79878.0]

T5_CI=[20.0 20.0 0.001 4.0395677495 0.003226437031110772 2583.0; 20.0 20.0 1.0e-5 1.9281875505000001 9.75905791229906e-5 6009.0; 20.0 20.0 1.0e-7 4.150009649999999 9.679217365759696e-7 12918.0; 20.0 20.0 1.0e-9 7.8049161005 9.813095506090672e-9 25221.0]

Plots.plot(BS3_CI[:,3],hcat(BS3_CI[:,5],RK4_CI[:,5],T5_CI[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["BS3 CI" "RK4 CI" "T5 CI"],markershapes= [:circle],xlabel = "solver tolerance",ylabel = "eigenvalue error")

Plots.plot(hcat(BS3_CI[:,6],RK4_CI[:,6],T5_CI[:,6]),hcat(BS3_CI[:,5],RK4_CI[:,5],T5_CI[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["BS3 CI" "RK4 CI" "T5 CI"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")