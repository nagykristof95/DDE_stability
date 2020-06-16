using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using Plots

A_LM4=[12.0 12.0 100.0 0.15244792005000002 0.00029657392095018347 1152.0; 12.0 12.0 200.0 0.3180125602 1.9038758038258767e-5 2352.0; 12.0 12.0 300.0 0.48387656504999993 3.8023517429020925e-6 3552.0; 12.0 12.0 400.0 0.54394977005 1.2121189768469364e-6 4752.0; 12.0 12.0 500.0 0.6411981899500001 4.995390976932534e-7 5952.0; 12.0 12.0 600.0 0.7681994301 2.4222868695566743e-7 7152.0; 12.0 12.0 700.0 0.8446453148999999 1.3141993530924302e-7 8352.0; 12.0 12.0 800.0 1.0007487149999998 7.741394688808659e-8 9552.0; 12.0 12.0 900.0 1.07686980515 4.855916244370232e-8 10752.0; 12.0 12.0 1000.0 1.21031186995 3.2008058267957705e-8 11952.0; 12.0 12.0 1100.0 1.2747194348 2.1962159642496463e-8 13152.0; 12.0 12.0 1200.0 1.3429870450999999 1.557699441716101e-8 14352.0; 12.0 12.0 1300.0 1.4728817399499998 1.136004269733297e-8 15552.0; 12.0 12.0 1400.0 1.54594626515 8.483433314255473e-9 16752.0; 12.0 12.0 1500.0 1.6831084501 6.466025800898495e-9 17952.0]

A_LM3=[12.0 12.0 100.0 0.15373337504999998 0.002921893106020544 1164.0; 12.0 12.0 200.0 0.3152327299 0.00036787063194875056 2364.0; 12.0 12.0 300.0 0.44169671509999997 0.00010921233731986545 3564.0; 12.0 12.0 400.0 0.54987395525 4.611583822629367e-5 4764.0; 12.0 12.0 500.0 0.6623130101 2.362369814565118e-5 5964.0; 12.0 12.0 600.0 0.7618607254 1.3675792525273766e-5 7164.0; 12.0 12.0 700.0 0.8055690802000001 8.614229691843837e-6 8364.0; 12.0 12.0 800.0 0.90629190985 5.77188888872731e-6 9564.0; 12.0 12.0 900.0 1.03931191025 4.054338916500643e-6 10764.0; 12.0 12.0 1000.0 1.13263221495 2.9559357228816937e-6 11964.0; 12.0 12.0 1100.0 1.25825659985 2.221035557384589e-6 13164.0; 12.0 12.0 1200.0 1.3810276397999999 1.7108890131906197e-6 14364.0; 12.0 12.0 1300.0 1.51222666485 1.3457441536234498e-6 15564.0; 12.0 12.0 1400.0 1.56900216025 1.0775352977899293e-6 16764.0; 12.0 12.0 1500.0 1.6209881599 8.761162975993878e-7 17964.0]

A_LM2=[12.0 12.0 100.0 0.1351545751 0.030524424133213175 1176.0; 12.0 12.0 200.0 0.28938839005 0.007568222008120016 2376.0; 12.0 12.0 300.0 0.40193740470000006 0.0033577959967817116 3576.0; 12.0 12.0 400.0 0.5061780899 0.0018874464438511538 4776.0; 12.0 12.0 500.0 0.6013010551000001 0.0012075271092535412 5976.0; 12.0 12.0 600.0 0.6987507851 0.0008383760964680541 7176.0; 12.0 12.0 700.0 0.7920716902 0.0006158596425670004 8376.0; 12.0 12.0 800.0 0.8954970450999999 0.0004714685241843788 9576.0; 12.0 12.0 900.0 0.9991502851 0.00037248949540301266 10776.0; 12.0 12.0 1000.0 1.09305243995 0.0003016984477773633 11976.0; 12.0 12.0 1100.0 1.25940019015 0.0002493257118326122 13176.0; 12.0 12.0 1200.0 1.3723345500000002 0.00020949476134449795 14376.0; 12.0 12.0 1300.0 1.3925582551 0.0001784986989214313 15576.0; 12.0 12.0 1400.0 1.513670585 0.00015390544682470562 16776.0; 12.0 12.0 1500.0 1.7812171849999998 0.0001340656648240304 17976.0]

A_LM1=[12.0 12.0 100.0 0.14738722005 0.42198673437932827 1188.0; 12.0 12.0 200.0 0.3273479202 0.18846417193568726 2388.0; 12.0 12.0 300.0 0.42544047479999997 0.12114751285245756 3588.0; 12.0 12.0 400.0 0.55559845005 0.08924110641585149 4788.0; 12.0 12.0 500.0 0.66714805495 0.0706318979552537 5988.0; 12.0 12.0 600.0 0.79357214015 0.05844272796705494 7188.0; 12.0 12.0 700.0 0.85505921995 0.04984066207485629 8388.0; 12.0 12.0 800.0 0.9482473200999999 0.04344557808383522 9588.0; 12.0 12.0 900.0 1.0650116098 0.038504791825125256 10788.0; 12.0 12.0 1000.0 1.1247951399000002 0.03457291728573917 11988.0; 12.0 12.0 1100.0 1.5003923901499998 0.03136957708909435 13188.0; 12.0 12.0 1200.0 1.56757242005 0.028709468766856594 14388.0; 12.0 12.0 1300.0 1.52323855 0.02646521692894129 15588.0; 12.0 12.0 1400.0 1.71102701005 0.0245463795431438 16788.0; 12.0 12.0 1500.0 1.7127235600500001 0.022886966935605532 17988.0]

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_LM4[:,6],A_LM3[:,6],A_LM2[:,6],A_LM1[:,6]),hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,4],A_LM3[:,4],A_LM2[:,4],A_LM1[:,4]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

open("temp.txt", "w") do io
           writedlm(io, A_LM1)
       end

Plots.savefig("1.pdf")
