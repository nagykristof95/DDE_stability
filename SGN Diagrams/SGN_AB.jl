using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using Plots

#point A

A_LM4=[12.0 12.0 100.0 0.5624203505 0.014262304271152649 1152.0; 12.0 12.0 200.0 0.36035665 0.006950957685270864 2352.0; 12.0 12.0 300.0 0.4980449 0.004614810175545371 3552.0; 12.0 12.0 400.0 0.7282423499999999 0.0034558041623685976 4752.0; 12.0 12.0 500.0 0.8098748995 0.0027624571857962498 5952.0; 12.0 12.0 600.0 0.8679331004999999 0.0023009333452921814 7152.0; 12.0 12.0 700.0 1.001938349 0.0019715816320458536 8352.0; 12.0 12.0 800.0 1.1823405505 0.0017247239435796528 9552.0; 12.0 12.0 900.0 1.2664887005 0.0015328112736589855 10752.0; 12.0 12.0 1000.0 1.523461551 0.0013793343780429128 11952.0; 12.0 12.0 1100.0 1.6329814995 0.0012537966055333196 13152.0; 12.0 12.0 1200.0 1.61558015 0.001149204837417378 14352.0; 12.0 12.0 1300.0 1.72648565 0.0010607202050090512 15552.0; 12.0 12.0 1400.0 1.8817130995 0.0009848878246540857 16752.0; 12.0 12.0 1500.0 1.9637306005 0.000919174991935756 17952.0]

A_LM3=[12.0 12.0 100.0 0.196012049 0.01319508076897506 1164.0; 12.0 12.0 200.0 0.3382395 0.006830861546236777 2364.0; 12.0 12.0 300.0 0.47613189899999997 0.0045822000791928925 3564.0; 12.0 12.0 400.0 0.5991737500000001 0.003442800986169693 4764.0; 12.0 12.0 500.0 0.7467199005 0.0027560519732662946 5964.0; 12.0 12.0 600.0 0.8387081 0.0022973290063108883 7164.0; 12.0 12.0 700.0 0.9557902995 0.001969359440553091 8364.0; 12.0 12.0 800.0 1.1023832494999999 0.0017232597296008441 9564.0; 12.0 12.0 900.0 1.2247487989999999 0.0015317965181468464 10764.0; 12.0 12.0 1000.0 1.3394838495 0.0013786026663858508 11964.0; 12.0 12.0 1100.0 1.4453792010000002 0.0012532518590934536 13164.0; 12.0 12.0 1200.0 1.5648962995 0.0011487884815108518 14364.0; 12.0 12.0 1300.0 1.704438549 0.0010603949004541972 15564.0; 12.0 12.0 1400.0 1.8115899500000001 0.000984628866610693 16764.0; 12.0 12.0 1500.0 1.9594921495 0.0009189655113660735 17964.0]

A_LM2=[12.0 12.0 100.0 0.1751715505 0.019113796844466343 1176.0; 12.0 12.0 200.0 0.33724849999999995 0.001806798814611713 2376.0; 12.0 12.0 300.0 0.46411060000000004 0.0012429604837732354 3576.0; 12.0 12.0 400.0 0.5918292505 0.0015025374326763526 4776.0; 12.0 12.0 500.0 0.7137538999999999 0.0015021521858069205 5976.0; 12.0 12.0 600.0 0.8335907 0.0014224110960728475 7176.0; 12.0 12.0 700.0 0.9549121 0.0013247258503392226 8376.0; 12.0 12.0 800.0 1.0781949005 0.0012287679310287831 9576.0; 12.0 12.0 900.0 1.2052348995 0.00114054994460709 10776.0; 12.0 12.0 1000.0 1.3074623 0.0010613643380501328 11976.0; 12.0 12.0 1100.0 1.4128845509999999 0.0009908588149645113 13176.0; 12.0 12.0 1200.0 1.5326714495 0.0009281616877119899 14376.0; 12.0 12.0 1300.0 1.6693709005 0.0008723044770025309 15576.0; 12.0 12.0 1400.0 1.7971526 0.0008223761523514989 16776.0; 12.0 12.0 1500.0 2.0007191505 0.0007775715619919001 17976.0]

A_LM1=[12.0 12.0 100.0 0.1598856 0.42799842087521406 1188.0; 12.0 12.0 200.0 0.32320595 0.19135609450354826 2388.0; 12.0 12.0 300.0 0.4703958495 0.12304148154901659 3588.0; 12.0 12.0 400.0 0.5833893 0.09064790496957255 4788.0; 12.0 12.0 500.0 0.7316056 0.07175056036192064 5988.0; 12.0 12.0 600.0 0.8520151 0.059371114007087666 7188.0; 12.0 12.0 700.0 0.96827915 0.05063404790727091 8388.0; 12.0 12.0 800.0 1.2358378 0.044138220477127905 9588.0; 12.0 12.0 900.0 2.2308820505 0.03911938196122327 10788.0; 12.0 12.0 1000.0 1.3550785505 0.03512525861507974 11988.0; 12.0 12.0 1100.0 1.4373759504999999 0.03187111601813023 13188.0; 12.0 12.0 1200.0 1.5209255005 0.0291687611263306 14388.0; 12.0 12.0 1300.0 1.6587872 0.02688882548583083 15588.0; 12.0 12.0 1400.0 2.4991974505 0.024939448324970095 16788.0; 12.0 12.0 1500.0 1.7746502995000002 0.023253602616277968 17988.0]

Plots.plot(A_LM4[:,3],hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error",markerstrokewidth = 0,markersize = 5)

Plots.plot(hcat(A_LM4[:,6],A_LM3[:,6],A_LM2[:,6],A_LM1[:,6]),hcat(A_LM4[:,5],A_LM3[:,5],A_LM2[:,5],A_LM1[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["LM4" "LM3" "LM2" "LM1"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error",markerstrokewidth = 0,markersize = 5)

open("temp.txt", "w") do io
           writedlm(io, A_LM1)
       end

Plots.savefig("1.pdf")