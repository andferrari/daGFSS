@everywhere using LightGraphs
using GraphIO
using EzXML
@everywhere using LinearAlgebra
using Plots
@everywhere using LaTeXStrings
@everywhere using DelimitedFiles
@everywhere using Statistics


@everywhere path = "/home/verduci/daGFSS/"
# load graph

@everywhere include(path*"gfss_func.jl")
include(path*"performance_func.jl")
include(path*"detection_func.jl")
@everywhere include(path*"signaux_func.jl")

g = loadgraph(path*"donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
#L = LightGraphs.laplacian_matrix(g)
L = NormalizedLaplacian(g)
A=adjacency_matrix(g)
d, v = eigen(Array(L));
λmax = maximum(d)

L= L - (λmax/2)I
node_labels = Int.(label_propagation(g, 10000)[1])


nb=10000 #nb de signaux créés
nt=512 #temps d'étude
ρ = 1.0/0.3
λ=0.01
Λ=0.1

#lancer synth_gfss_parallele pour avoir les poles et les residues
poles=[0.7344635585485741 + 0.5293681924014867im 0.7344635585485741 - 0.5293681924014867im -0.7344635202651353 + 0.5293681751822871im -0.7344635202651353 - 0.5293681751822871im]
residues=[-0.05322205451164147 - 0.08206089787078102im -0.05322205451164147 + 0.08206089787078102im 0.16044187053514897 - 0.6853621323079733im 0.16044187053514897 + 0.6853621323079733im]
φ, ψ = calcul_psi_phi(poles, residues)
c=0.6682305081233931
variancei=readdlm("variancei.csv")
σi=sqrt.(variancei)





sig1_sc=SharedArray{Float64}(250, nt*nb)
t_diaGFSS_sc=SharedArray{Float64}(250, nt*nb)
t_daGFSS_sc=SharedArray{Float64}(nb,nt)
t_carre_diaGFSS_sc=SharedArray{Float64}(250, nt*nb)
t_carre_daGFSS_sc=SharedArray{Float64}(nb,nt)
t_daGFSS_carre_sc=SharedArray{Float64}(nb,nt)
sumtot_daGFSS_sc=SharedArray{Float64}(nb,nt)



#création des signaux
@distributed for k in 1:nb
    sig1_sc[:,(k-1)*512+1:k*512]= gener_sigg(g, node_labels, 1.0, 2 , 3.0, σ2 = 7);
    t_carre_diaGFSS_sc[:,(k-1)*512+1:k*512], t_diaGFSS_sc[:,(k-1)*512+1:k*512] = diaGFSS(sig1_sc[:,(k-1)*512+1:k*512], L, ψ, φ, c; λ = 0.01, Λ=0.1)
end

t_voisin_sc=-(t_diaGFSS_sc + A*t_diaGFSS_sc)
t_voisin_carre_sc = (t_voisin_sc.^2)./σi

t_carre_voisin_sc=(t_carre_diaGFSS_sc + A*t_carre_diaGFSS_sc)./σi
t_carre_diaGFSS_sc=t_carre_diaGFSS_sc./σi

@distributed for k in 1:nb
    t_carre_daGFSS_sc[k,:] = [norm(t_carre_voisin_sc[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS_carre_sc[k,:] = [norm(t_voisin_carre_sc[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS_sc[k,:] = [norm(t_carre_diaGFSS_sc[:,i])^2 for i in (k-1)*512+1:k*512]
end

@distributed for k in 1:nb
    sumtot_daGFSS_sc[k,:]=sum(t_carre_diaGFSS_sc[:,(k-1)*512+1:k*512],dims=1)
end
sumtot_daGFSS_sc





sig1=SharedArray{Float64}(250, nt*nb)
t_diaGFSS=SharedArray{Float64}(250, nt*nb)
t_daGFSS=SharedArray{Float64}(nb,nt)
t_carre_diaGFSS=SharedArray{Float64}(250, nt*nb)
t_carre_daGFSS=SharedArray{Float64}(nb,nt)
t_daGFSS_carre=SharedArray{Float64}(nb,nt)
sumtot_daGFSS=SharedArray{Float64}(nb,nt)


#création des signaux
@distributed for k in 1:nb
    sig1[:,(k-1)*512+1:k*512]= gener_sigg(g, node_labels, 1.0, 1 , 3.0, σ2 = 7);
    t_carre_diaGFSS[:,(k-1)*512+1:k*512], t_diaGFSS[:,(k-1)*512+1:k*512] = diaGFSS(sig1[:,(k-1)*512+1:k*512], L, ψ, φ, c; λ = 0.01, Λ=0.1)
end

t_voisin=-(t_diaGFSS + A*t_diaGFSS)
t_voisin_carre = (t_voisin.^2)./σi

t_carre_voisin=(t_carre_diaGFSS + A*t_carre_diaGFSS)./σi
t_carre_diaGFSS=t_carre_diaGFSS./σi

@distributed for k in 1:nb
    t_carre_daGFSS[k,:] = [norm(t_carre_voisin[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS_carre[k,:] = [norm(t_voisin_carre[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS[k,:] = [norm(t_carre_diaGFSS[:,i])^2 for i in (k-1)*512+1:k*512]
end
t_carre_daGFSS

@distributed for k in 1:nb
    sumtot_daGFSS[k,:]=sum(t_carre_diaGFSS[:,(k-1)*512+1:k*512],dims=1)
end







#t_carre
init=300
init_sc=350
fin=512
plot(t_carre_daGFSS_sc[100,init_sc:fin],xlabel="temps", ylabel="t_carre_aGFSS")
n_rupt = 400
Δ_rupt = 30
pas=10
mini=round(minimum(t_carre_daGFSS_sc[:,init_sc:fin]))
maxi=round(maximum(t_carre_daGFSS_sc[:,init_sc:fin]))+1000
x1_sc=(mini):pas:maxi
T100=tab_threshold(mini,maxi,pas,nt) #creation d'un tableau de seuils

pfausse, tdetecttot= performance_algo(nt, t_carre_daGFSS_sc, ρ ,d ,v ,T100)


plot(t_carre_daGFSS[9000,init:fin],xlabel="temps", ylabel="t_carre_daGFSS")
mini=round(minimum(t_carre_daGFSS[:,init:fin]))
maxi=round(maximum(t_carre_daGFSS[:,init:fin]))+100000
pas=(maxi-mini)/size(T100)[1]
x1=(mini):pas:maxi
T1000=tab_threshold(mini,maxi,pas,nt) #creation d'un tableau de seuils

plot(t_carre_voisin[:,init:fin]',xlabel="temps", ylabel="t_carre_aGFSS", label="")
T20001=80*ones(250,512)

pdetect, retardtot, tdetecttot= performance_algo2(nt, sig1, t_carre_daGFSS, t_carre_voisin, ρ ,d ,v ,T1000, T20001)






#t_voisin_carre
plot(t_daGFSS_carre_sc[100,init_sc:fin],xlabel="temps", ylabel="t_daGFSS_carre")
pas=1000
mini2=round(minimum(t_daGFSS_carre_sc[:,init_sc:fin]))
maxi2=round(maximum(t_daGFSS_carre_sc[:,init_sc:fin]))+100000
x2=(mini2):pas:maxi2
T1002=tab_threshold(mini2,maxi2,pas,nt)


pfausse2, tdetecttot2= performance_algo(nt, t_daGFSS_carre_sc, ρ ,d ,v ,T1002,)


plot(t_daGFSS_carre[10,init:fin],xlabel="temps", ylabel="t_daGFSS_carre")
mini2=round(minimum(t_daGFSS_carre[:,init:fin]))
maxi2=round(maximum(t_daGFSS_carre[:,init:fin]))
pas=(maxi2-mini2)/size(T1002)[1]
x2=(mini2):pas:maxi2
T10002=tab_threshold(mini2,maxi2,pas,nt)

plot(t_voisin_carre[:,init:fin]',xlabel="temps", ylabel="t_carre_aGFSS", label="")
T2000=1000*ones(250,512)

pdetect2, retardtot2, tdetecttot2= performance_algo2(nt, sig1, t_daGFSS_carre, t_voisin_carre, ρ ,d ,v ,T10002, T2000)





#t
plot(t_daGFSS_sc[100,init_sc:fin],xlabel="temps", ylabel="t_daGFSS")
pas=0.1
mini3=round(minimum(t_daGFSS_sc[:,init_sc:fin]))
maxi3=round(maximum(t_daGFSS_sc[:,init_sc:fin]))+10
x3=(mini3):pas:maxi3
T1003=tab_threshold(mini3,maxi3,pas,nt)


pfausse3, tdetecttot3= performance_algo(nt, t_daGFSS_sc, ρ ,d ,v ,T1003)


plot(t_daGFSS[100,init:fin],xlabel="temps", ylabel="t_daGFSS")
mini3=round(minimum(t_daGFSS[:,init:fin]))
maxi3=round(maximum(t_daGFSS[:,init:fin]))
pas=(maxi3-mini3)/size(T1003)[1]
x3=(mini3):pas:maxi3
T10003=tab_threshold(mini3,maxi3,pas,nt)

plot(t_carre_diaGFSS[:,init:fin]',xlabel="temps", ylabel="t_diaGFSS", label="")
T20002=10*ones(250,512)

pdetect3, retardtot3, tdetecttot3= performance_algo2(nt, sig1, t_daGFSS, t_carre_diaGFSS, ρ ,d ,v ,T10003, T20002)






#t_sum

plot(sumtot_daGFSS_sc[100,init_sc:fin],xlabel="temps", ylabel="t_sum_daGFSS")
pas=0.1
mini4=round(minimum(sumtot_daGFSS_sc[:,init_sc:fin]))
maxi4=round(maximum(sumtot_daGFSS_sc[:,init_sc:fin]))
x4=(mini4):pas:maxi4
T1004=tab_threshold(mini4,maxi4+6,pas,nt) #creation d'un tableau de seuils


pfausse4, tdetecttot4= performance_algo(nt, sumtot_daGFSS_sc, ρ ,d ,v ,T1004)


plot(sumtot_daGFSS[100,init:fin],xlabel="temps", ylabel="t_sum_daGFSS")
mini4=round(minimum(sumtot_daGFSS[:,init:fin]))
maxi4=round(maximum(sumtot_daGFSS[:,init:fin]))
pas=(maxi4-mini4)/size(T1004)[1]
x4=(mini4):pas:maxi4
T10004=tab_threshold(mini4,maxi4,pas,nt) #creation d'un tableau de seuils


pdetect4, retardtot4, tdetecttot4= performance_algo2(nt, sig1, sumtot_daGFSS, t_carre_diaGFSS, ρ ,d ,v ,T10004, T20002)

pyplot()
plot(pfausse2, pdetect2, ratio=:equal, xlabel=L"pfa", ylabel=L"pd", label="daGFSS", w=3)
plot!(pfausse3, pdetect3, label="daGFSS Independent", w=3)
plot!(pfausse, pdetect, label="daGFSS with 2-norm", w=3)
plot!(pfausse4, pdetect4, label="daGFSS Centralized", dpi=600, w=3)


plot(pfausse2,retardtot2, xlabel="pfa", ylabel ="retard", label="daGFSS")
plot!(pfausse3, retardtot3, label="daGFSS Independent")
plot!(pfausse, retardtot, label="daGFSS with 2-norm")
plot!(pfausse4, retardtot4, label="daGFSS Centralized")
