using Printf
using Statistics
using LinearAlgebra
using Interpolations
using DelimitedFiles
using Plots

#constant
A_RK4_k=[12.0 12.0 100.0 0.3001772 0.020154088383299087 4752.0; 12.0 12.0 200.0 0.4540847 0.013040686850559459 9552.0; 12.0 12.0 300.0 0.6705201 0.0076486298114882435 14352.0; 12.0 12.0 400.0 0.8481294 0.005781068962432587 19152.0; 12.0 12.0 500.0 1.066977199 0.005850347629248018 23952.0; 12.0 12.0 600.0 1.320111501 0.004232921493904526 28752.0; 12.0 12.0 700.0 1.411466099 0.004225601007035612 33552.0; 12.0 12.0 800.0 1.576926899 0.0026373990237438666 38352.0; 12.0 12.0 900.0 1.7957894 0.0022955328865696175 43152.0; 12.0 12.0 1000.0 1.9996121 0.0029814687425695164 47952.0; 12.0 12.0 1100.0 2.187613899 0.0025432167123324296 52752.0; 12.0 12.0 1200.0 2.390638599 0.0019833980325676412 57552.0; 12.0 12.0 1300.0 2.5426643 0.001699752827036287 62352.0; 12.0 12.0 1400.0 2.8340848 0.0020547592327691713 67152.0; 12.0 12.0 1500.0 2.9579473 0.0015007475064067806 71952.0]

A_RK3_k=[12.0 12.0 100.0 0.219798301 0.019985050831408623 3564.0; 12.0 12.0 200.0 0.390588301 0.013018604050873468 7164.0; 12.0 12.0 300.0 0.5553554 0.007642617291532883 10764.0; 12.0 12.0 400.0 0.6803983 0.005778602800814641 14364.0; 12.0 12.0 500.0 0.878712199 0.00584891619997967 17964.0; 12.0 12.0 600.0 1.0054932 0.004232162505824564 21564.0; 12.0 12.0 700.0 1.177468699 0.0042251097084903206 25164.0; 12.0 12.0 800.0 1.364425299 0.0026371078237424116 28764.0; 12.0 12.0 900.0 1.4781928 0.002295317801526308 32364.0; 12.0 12.0 1000.0 1.641827601 0.002981289705721057 35964.0; 12.0 12.0 1100.0 1.8128214 0.0025430850553984303 39564.0; 12.0 12.0 1200.0 1.9508422 0.001983304836240532 43164.0; 12.0 12.0 1300.0 2.158966699 0.0016996829908170353 46764.0; 12.0 12.0 1400.0 2.31962 0.0020546955657387523 50364.0; 12.0 12.0 1500.0 2.3923687 0.0015007015236336281 53964.0]

A_RK2_k=[12.0 12.0 100.0 0.1789601 0.019964737879764376 2376.0; 12.0 12.0 200.0 0.3092584 0.011704731984970787 4776.0; 12.0 12.0 300.0 0.4153872 0.008161992392891459 7176.0; 12.0 12.0 400.0 0.528091699 0.00625235632528909 9576.0; 12.0 12.0 500.0 0.634553401 0.005059454087098762 11976.0; 12.0 12.0 600.0 0.7913934 0.004234121918160055 14376.0; 12.0 12.0 700.0 0.865565801 0.003636151005953395 16776.0; 12.0 12.0 800.0 0.9854371 0.003179809514721006 19176.0; 12.0 12.0 900.0 1.116183599 0.002822603853933634 21576.0; 12.0 12.0 1000.0 1.2280337 0.0025348841575051056 23976.0; 12.0 12.0 1100.0 1.353228401 0.002297043249712144 26376.0; 12.0 12.0 1200.0 1.4602377 0.0020975974000462774 28776.0; 12.0 12.0 1300.0 1.603789699 0.0019282598128148546 31176.0; 12.0 12.0 1400.0 1.701069201 0.0017828902156962548 33576.0; 12.0 12.0 1500.0 1.8562789 0.0016560373264109104 35976.0]

A_RK1_k=[12.0 12.0 100.0 0.380195801 0.423504069570898 1188.0; 12.0 12.0 200.0 0.2198575 0.19766224791363599 2388.0; 12.0 12.0 300.0 0.2938553 0.12453486612218413 3588.0; 12.0 12.0 400.0 0.368179499 0.09128533514122136 4788.0; 12.0 12.0 500.0 0.438204901 0.07580135600613912 5988.0; 12.0 12.0 600.0 0.4747741 0.060866948601977035 7188.0; 12.0 12.0 700.0 0.5771651 0.05257768151793855 8388.0; 12.0 12.0 800.0 0.6380724 0.04390274933119775 9588.0; 12.0 12.0 900.0 0.7049522 0.03953357726011462 10788.0; 12.0 12.0 1000.0 0.777064999 0.03737580040503891 11988.0; 12.0 12.0 1100.0 0.8765139 0.03361934802249949 13188.0; 12.0 12.0 1200.0 0.9296702 0.02984011837102506 14388.0; 12.0 12.0 1300.0 1.043836801 0.027031599567435526 15588.0; 12.0 12.0 1400.0 1.1031565 0.02631177642037721 16788.0; 12.0 12.0 1500.0 1.1489548 0.023481226687698926 17988.0]

Plots.plot(A_RK4_k[:,3],hcat(A_RK4_k[:,5],A_RK3_k[:,5],A_RK2_k[:,5],A_RK1_k[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4_k[:,6],A_RK3_k[:,6],A_RK2_k[:,6],A_RK1_k[:,6]),hcat(A_RK4_k[:,5],A_RK3_k[:,5],A_RK2_k[:,5],A_RK1_k[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")


#Linear
A_RK4_l=[12.0 12.0 100.0 0.307561099 0.0004950195208625873 4752.0; 12.0 12.0 200.0 0.524955 0.00012180966608616196 9552.0; 12.0 12.0 300.0 0.707139999 5.389814207118713e-5 14352.0; 12.0 12.0 400.0 0.9339089 3.025195222200019e-5 19152.0; 12.0 12.0 500.0 1.2264984 1.9334654865257502e-5 23952.0; 12.0 12.0 600.0 1.346544699 1.3413129850882598e-5 28752.0; 12.0 12.0 700.0 1.695238001 9.846179302147852e-6 33552.0; 12.0 12.0 800.0 1.8104945 7.532729271764562e-6 38352.0; 12.0 12.0 900.0 2.2945431 5.947476800804724e-6 43152.0; 12.0 12.0 1000.0 2.3971442 4.814022784624226e-6 47952.0; 12.0 12.0 1100.0 2.4475224 3.975671797438925e-6 52752.0; 12.0 12.0 1200.0 2.969103801 3.3382093628968554e-6 57552.0; 12.0 12.0 1300.0 3.0132335 2.842225405982911e-6 62352.0; 12.0 12.0 1400.0 3.2277246 2.448752331512803e-6 67152.0; 12.0 12.0 1500.0 3.411795599 2.1313702132128034e-6 71952.0]

A_RK3_l=[12.0 12.0 100.0 0.1848875 0.0007596057727981994 3564.0; 12.0 12.0 200.0 0.379962601 0.00015418487626335024 7164.0; 12.0 12.0 300.0 0.544989101 6.340513095370348e-5 10764.0; 12.0 12.0 400.0 0.7054044 3.424299045286556e-5 14364.0; 12.0 12.0 500.0 0.887602699 2.1371708118615105e-5 17964.0; 12.0 12.0 600.0 1.072424101 1.4589450724198043e-5 21564.0; 12.0 12.0 700.0 1.3037797 1.0585790742943037e-5 25164.0; 12.0 12.0 800.0 1.424649099 8.027618241820525e-6 28764.0; 12.0 12.0 900.0 1.611686099 6.294725474517406e-6 32364.0; 12.0 12.0 1000.0 1.81946 5.066974244812444e-6 35964.0; 12.0 12.0 1100.0 1.952489199 4.1655985952655445e-6 39564.0; 12.0 12.0 1200.0 2.110584201 3.484424298156023e-6 43164.0; 12.0 12.0 1300.0 2.2720814 2.957175892334238e-6 46764.0; 12.0 12.0 1400.0 2.538111701 2.5407526383281476e-6 50364.0; 12.0 12.0 1500.0 2.6356298 2.206144791150234e-6 53964.0]

A_RK2_l=[12.0 12.0 100.0 0.124328701 0.011439310404940234 2376.0; 12.0 12.0 200.0 0.256725601 0.0028277198636559304 4776.0; 12.0 12.0 300.0 0.384867399 0.0012521085374732936 7176.0; 12.0 12.0 400.0 0.5097898 0.0007030080279300505 9576.0; 12.0 12.0 500.0 0.625142801 0.00044942654676232016 11976.0; 12.0 12.0 600.0 0.781699899 0.0003118719482023613 14376.0; 12.0 12.0 700.0 0.933799199 0.00022901052803226145 16776.0; 12.0 12.0 800.0 0.981786401 0.00017526790487406763 19176.0; 12.0 12.0 900.0 1.171375 0.0001384417691463279 21576.0; 12.0 12.0 1000.0 1.2578498 0.00011211129939900819 23976.0; 12.0 12.0 1100.0 1.399598499 9.263633179209416e-5 26376.0; 12.0 12.0 1200.0 1.4846427 7.78281629750683e-5 28776.0; 12.0 12.0 1300.0 1.651582501 6.630662933871004e-5 31176.0; 12.0 12.0 1400.0 1.820846801 5.716647081862856e-5 33576.0; 12.0 12.0 1500.0 1.8603958 4.9793926106741e-5 35976.0]

A_RK1_l=[12.0 12.0 100.0 0.064360499 0.4219867429471442 1188.0; 12.0 12.0 200.0 0.169368799 0.1884641808708395 2388.0; 12.0 12.0 300.0 0.2504456 0.12114752189373915 3588.0; 12.0 12.0 400.0 0.3119551 0.08924111550745852 4788.0; 12.0 12.0 500.0 0.382674 0.07063190707621704 5988.0; 12.0 12.0 600.0 0.5038901 0.05844273710724676 7188.0; 12.0 12.0 700.0 0.5635994 0.04984067122862143 8388.0; 12.0 12.0 800.0 0.6323597 0.04344558724769039 9588.0; 12.0 12.0 900.0 0.746412499 0.038504800996777715 10788.0; 12.0 12.0 1000.0 0.7918693 0.034572926463594104 11988.0; 12.0 12.0 1100.0 0.908460199 0.03136958627200163 13188.0; 12.0 12.0 1200.0 0.9291919 0.028709477953962784 14388.0; 12.0 12.0 1300.0 1.006973899 0.02646522611958815 15588.0; 12.0 12.0 1400.0 1.1449959 0.024546388736819028 16788.0; 12.0 12.0 1500.0 1.1479198 0.022886976131897132 17988.0]

Plots.plot(A_RK4_l[:,3],hcat(A_RK4_l[:,5],A_RK3_l[:,5],A_RK2_l[:,5],A_RK1_l[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4_l[:,6],A_RK3_l[:,6],A_RK2_l[:,6],A_RK1_l[:,6]),hcat(A_RK4_l[:,5],A_RK3_l[:,5],A_RK2_l[:,5],A_RK1_l[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#Quadratic
A_RK4_q=[12.0 12.0 100.0 0.2540405 9.687573612477964e-6 4752.0; 12.0 12.0 200.0 0.501476 8.012730550722566e-7 9552.0; 12.0 12.0 300.0 0.742773401 2.0179085790054862e-7 14352.0; 12.0 12.0 400.0 0.9611522 7.823021000233317e-8 19152.0; 12.0 12.0 500.0 1.1819103 3.806132893905593e-8 23952.0; 12.0 12.0 600.0 1.6279499 2.1292681554553938e-8 28752.0; 12.0 12.0 700.0 1.951496401 1.3090675549088612e-8 33552.0; 12.0 12.0 800.0 2.368146101 8.61447323980436e-9 38352.0; 12.0 12.0 900.0 1.9812795 5.9673871874272765e-9 43152.0; 12.0 12.0 1000.0 2.227182501 4.302831423267958e-9 47952.0; 12.0 12.0 1100.0 2.371742801 3.204105872008906e-9 52752.0; 12.0 12.0 1200.0 2.850372899 2.44980444114969e-9 57552.0; 12.0 12.0 1300.0 2.858852301 1.914865065442184e-9 62352.0; 12.0 12.0 1400.0 3.285120501 1.524998696068889e-9 67152.0; 12.0 12.0 1500.0 3.394366201 1.234168322596005e-9 71952.0]

A_RK3_q=[12.0 12.0 100.0 0.1965645 0.00031210568265360516 3564.0; 12.0 12.0 200.0 0.4066025 3.8378908614216294e-5 7164.0; 12.0 12.0 300.0 0.557921299 1.1308898323325296e-5 10764.0; 12.0 12.0 400.0 0.7327941 4.7577465365745344e-6 14364.0; 12.0 12.0 500.0 0.9134612 2.4319162176269296e-6 17964.0; 12.0 12.0 600.0 1.0702027 1.4057973179959702e-6 21564.0; 12.0 12.0 700.0 1.3267672 8.845812655567742e-7 25164.0; 12.0 12.0 800.0 1.483490899 5.922476075020622e-7 28764.0; 12.0 12.0 900.0 1.6154227 4.1576174702742724e-7 32364.0; 12.0 12.0 1000.0 1.682690899 3.0297792978437483e-7 35964.0; 12.0 12.0 1100.0 1.894388699 2.2756269780858599e-7 39564.0; 12.0 12.0 1200.0 2.1021859 1.7523684513843254e-7 43164.0; 12.0 12.0 1300.0 2.2027622 1.3779896779676172e-7 46764.0; 12.0 12.0 1400.0 2.509451401 1.103092625936496e-7 50364.0; 12.0 12.0 1500.0 2.6764436 8.967119333566903e-8 53964.0]

A_RK2_q=[12.0 12.0 100.0 0.149586499 0.011882652290716953 2376.0; 12.0 12.0 200.0 0.303903001 0.00293626830447361 4776.0; 12.0 12.0 300.0 0.4351759 0.0012999925958695638 7176.0; 12.0 12.0 400.0 0.5692797 0.0007298377469827964 9576.0; 12.0 12.0 500.0 0.687481 0.00046655549359807415 11976.0; 12.0 12.0 600.0 0.8051783 0.00032374649867352584 14376.0; 12.0 12.0 700.0 0.9386754 0.00023772313928022325 16776.0; 12.0 12.0 800.0 1.0437889 0.00018193128610581946 19176.0; 12.0 12.0 900.0 1.1837501 0.00014370177139012375 21576.0; 12.0 12.0 1000.0 1.3290213 0.00011636835730762115 23976.0; 12.0 12.0 1100.0 1.388115001 9.615185414112711e-5 26376.0; 12.0 12.0 1200.0 1.559117 8.078001974444898e-5 28776.0; 12.0 12.0 1300.0 1.651180499 6.882003956087169e-5 31176.0; 12.0 12.0 1400.0 1.819090399 5.933213443598765e-5 33576.0; 12.0 12.0 1500.0 1.938710099 5.167914801824123e-5 35976.0]

A_RK1_q=[12.0 12.0 100.0 0.1015621 0.4219867343795105 1188.0; 12.0 12.0 200.0 0.2091522 0.188464171935867 2388.0; 12.0 12.0 300.0 0.313113201 0.1211475128526361 3588.0; 12.0 12.0 400.0 0.4339983 0.08924110641603036 4788.0; 12.0 12.0 500.0 0.4816458 0.07063189795543369 5988.0; 12.0 12.0 600.0 0.5716082 0.05844272796723331 7188.0; 12.0 12.0 700.0 0.701295999 0.04984066207503608 8388.0; 12.0 12.0 800.0 0.7218552 0.04344557808401611 9588.0; 12.0 12.0 900.0 0.8348313 0.03850479182530448 10788.0; 12.0 12.0 1000.0 0.9001587 0.034572917285919026 11988.0; 12.0 12.0 1100.0 0.9541737 0.03136957708927295 13188.0; 12.0 12.0 1200.0 1.0307493 0.028709468767034823 14388.0; 12.0 12.0 1300.0 1.1178865 0.026465216929119215 15588.0; 12.0 12.0 1400.0 1.2593112 0.024546379543322498 16788.0; 12.0 12.0 1500.0 1.394272401 0.022886966935783924 17988.0]

Plots.plot(A_RK4_q[:,3],hcat(A_RK4_q[:,5],A_RK3_q[:,5],A_RK2_q[:,5],A_RK1_q[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4_q[:,6],A_RK3_q[:,6],A_RK2_q[:,6],A_RK1_q[:,6]),hcat(A_RK4_q[:,5],A_RK3_q[:,5],A_RK2_q[:,5],A_RK1_q[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#Cubic
A_RK4_c=[12.0 12.0 100.0 0.312675001 8.90051272596629e-6 4752.0; 12.0 12.0 200.0 0.8531439 7.073404478427969e-7 9552.0; 12.0 12.0 300.0 0.7113148 1.7409079574517344e-7 14352.0; 12.0 12.0 400.0 0.9384229 6.657233337611436e-8 19152.0; 12.0 12.0 500.0 1.257447799 3.210682592371865e-8 23952.0; 12.0 12.0 600.0 1.3950216 1.7854953958989228e-8 28752.0; 12.0 12.0 700.0 1.6696661 1.0930625377525253e-8 33552.0; 12.0 12.0 800.0 1.8171005 7.1703481384475046e-9 38352.0; 12.0 12.0 900.0 2.137529399 4.954989869324417e-9 43152.0; 12.0 12.0 1000.0 2.3573754 3.566017484393808e-9 47952.0; 12.0 12.0 1100.0 2.755095501 2.6513498791130343e-9 52752.0; 12.0 12.0 1200.0 2.9065443 2.0246114001061092e-9 57552.0; 12.0 12.0 1300.0 2.794637099 1.5808489325497542e-9 62352.0; 12.0 12.0 1400.0 3.6175097 1.257861169767802e-9 67152.0; 12.0 12.0 1500.0 3.3577195 1.0172000727290635e-9 71952.0]

A_RK3_c=[12.0 12.0 100.0 0.1893095 0.0003113207566016173 3564.0; 12.0 12.0 200.0 0.5504046 3.829112462279342e-5 7164.0; 12.0 12.0 300.0 0.6124292 1.1283948729995178e-5 10764.0; 12.0 12.0 400.0 0.7295897 4.747447552383068e-6 14364.0; 12.0 12.0 500.0 1.4922435 2.426713239937589e-6 17964.0; 12.0 12.0 600.0 1.259444 1.4028135051542228e-6 21564.0; 12.0 12.0 700.0 1.3488605 8.827145165662846e-7 25164.0; 12.0 12.0 800.0 1.528305499 5.91003217019116e-7 28764.0; 12.0 12.0 900.0 1.938680301 4.1489116639329627e-7 32364.0; 12.0 12.0 1000.0 1.943448101 3.023452638490504e-7 35964.0; 12.0 12.0 1100.0 2.051746499 2.2708859541761336e-7 39564.0; 12.0 12.0 1200.0 2.120002701 1.7487245813067323e-7 43164.0; 12.0 12.0 1300.0 2.3492245 1.3751290259460676e-7 46764.0; 12.0 12.0 1400.0 2.538161801 1.1008059493957638e-7 50364.0; 12.0 12.0 1500.0 2.7155536 8.948554104732187e-8 53964.0]

A_RK2_c=[12.0 12.0 100.0 0.146271499 0.011883539806535831 2376.0; 12.0 12.0 200.0 0.2985744 0.0029363633591158525 4776.0; 12.0 12.0 300.0 0.464889299 0.0013000192373945412 7176.0; 12.0 12.0 400.0 0.570689101 0.0007298486688239726 9576.0; 12.0 12.0 500.0 0.7419946 0.00046656098861849537 11976.0; 12.0 12.0 600.0 0.8865098 0.000323749641390841 14376.0; 12.0 12.0 700.0 1.042985 0.00023772510163679194 16776.0; 12.0 12.0 800.0 1.1523367 0.00018193259233505212 19176.0; 12.0 12.0 900.0 1.2513982 0.000143702684212219 21576.0; 12.0 12.0 1000.0 1.333045 0.00011636902007869176 23976.0; 12.0 12.0 1100.0 1.553585699 9.615235044373296e-5 26376.0; 12.0 12.0 1200.0 1.6930002 8.078040096179931e-5 28776.0; 12.0 12.0 1300.0 1.9061578 6.882033869256197e-5 31176.0; 12.0 12.0 1400.0 2.0580764 5.933237344969e-5 33576.0; 12.0 12.0 1500.0 1.999093199 5.1679341998391315e-5 35976.0]

A_RK1_c=[12.0 12.0 100.0 0.208450899 0.4219867343794475 1188.0; 12.0 12.0 200.0 0.275412 0.18846417193580153 2388.0; 12.0 12.0 300.0 0.333938101 0.12114751285257036 3588.0; 12.0 12.0 400.0 0.446110101 0.08924110641596532 4788.0; 12.0 12.0 500.0 0.5458056 0.07063189795536637 5988.0; 12.0 12.0 600.0 0.6245838 0.0584427279671666 7188.0; 12.0 12.0 700.0 0.819312199 0.049840662074973466 8388.0; 12.0 12.0 800.0 0.844221999 0.043445578083945426 9588.0; 12.0 12.0 900.0 0.883153 0.03850479182523522 10788.0; 12.0 12.0 1000.0 0.964643099 0.03457291728585101 11988.0; 12.0 12.0 1100.0 1.026474201 0.03136957708920558 13188.0; 12.0 12.0 1200.0 1.073087001 0.028709468766968654 14388.0; 12.0 12.0 1300.0 1.1573125 0.026465216929054566 15588.0; 12.0 12.0 1400.0 1.267828601 0.024546379543256166 16788.0; 12.0 12.0 1500.0 1.3803196 0.022886966935718685 17988.0]

Plots.plot(A_RK4_c[:,3],hcat(A_RK4_c[:,5],A_RK3_c[:,5],A_RK2_c[:,5],A_RK1_c[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK4_c[:,6],A_RK3_c[:,6],A_RK2_c[:,6],A_RK1_c[:,6]),hcat(A_RK4_c[:,5],A_RK3_c[:,5],A_RK2_c[:,5],A_RK1_c[:,5]),xscale=:log10, yscale=:log10,title="Comparaison of methods A",label=["RK4" "RK3" "RK2" "EE"],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

#Comparison
Plots.plot(A_RK3_k[:,3],hcat(A_RK3_k[:,5],A_RK3_l[:,5],A_RK3_q[:,5],A_RK3_c[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK1 k" "RK3 l" "RK3 q" "RK3 c" ],markershapes= [:circle],xlabel = "number of steps",ylabel = "eigenvalue error")

Plots.plot(hcat(A_RK3_k[:,6],A_RK3_l[:,6],A_RK3_q[:,6],A_RK3_c[:,6]),hcat(A_RK3_k[:,5],A_RK3_l[:,5],A_RK3_q[:,5],A_RK3_c[:,5]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK1 k" "RK3 l" "RK3 q" "RK3 c" ],markershapes= [:circle],xlabel = "function calls",ylabel = "eigenvalue error")

Plots.plot(A_RK4_q[:,3],hcat(A_RK4_q[:,6],A_RK3_q[:,6],A_RK2_q[:,6],A_RK1_q[:,6]),xscale=:log10, yscale=:log10,title="Largest eigenvalue error A",label=["RK4 q" "RK3 q" "RK2 q" "RK1 q" ],markershapes= [:circle],xlabel = "number of steps",ylabel = "function calls")

open("temp.txt", "w") do io
           writedlm(io, A_RK4_k)
       end

Plots.savefig("1.pdf")