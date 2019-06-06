using LightGraphs
using GraphIO
using EzXML
using LinearAlgebra
using Plots

using LaTeXStrings
pyplot()

#using Compose
#using Cairo
using DelimitedFiles
using Statistics
#using Graphs


# load graph

include("gfss_func.jl")
include("performance_func.jl")
include("detection_func.jl")
include("signaux_func.jl")

g = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
#L = LightGraphs.laplacian_matrix(g)
L = NormalizedLaplacian(g)
A=adjacency_matrix(g)
d, v = eigen(Array(L));
λmax = maximum(d)

L= L - (λmax/2)I
node_labels = Int.(label_propagation(g, 10000)[1])


nb=100 #nb de signaux créés
nt=512 #temps d'étude
ρ = 1.0/0.3

#lancer synth_gfss_parallele pour avoir les poles et les residues
φ, ψ = calcul_psi_phi(poles, residues)
c=addTerm

sig1=zeros(250, nt*nb)
t_diaGFSS=zeros(250, nt*nb)
t_daGFSS=zeros(nb,nt)
t_carre_diaGFSS=zeros(250, nt*nb)
t_carre_daGFSS=zeros(nb,nt)
t_daGFSS_carre=zeros(nb,nt)



#création des signaux
for k in 1:nb
    sig1[:,(k-1)*512+1:k*512]= gener_sigg(g, node_labels, 3.0, 1 , 3.0, σ2 = 7);
    t_carre_diaGFSS[:,(k-1)*512+1:k*512], t_diaGFSS[:,(k-1)*512+1:k*512] = diaGFSS(sig1[:,(k-1)*512+1:k*512], L, ψ, φ, c; λ = 0.01, Λ=0.1)
end


#writedlm( "sig1_3_7.csv",  sig1)
#writedlm( "t_diaGFSS_3_7.csv",  t_diaGFSS)
#writedlm( "t_carre_diaGFSS_3_7.csv",  t_carre_diaGFSS)

#sig1=readdlm("sig1_3_7.csv")
#t_diaGFSS = readdlm("t_diaGFSS_3_7.csv")
#t_carre_diaGFSS = readdlm("t_carre_diaGFSS_3_7.csv")


t_voisin=-(t_diaGFSS + A*t_diaGFSS)
t_voisin_carre = t_voisin.^2

t_carre_voisin=(t_carre_diaGFSS + A*t_carre_diaGFSS)

for k in 1:nb
    t_carre_daGFSS[k,:] = [norm(t_carre_voisin[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS_carre[k,:] = [norm(t_voisin_carre[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS[k,:] = [norm(t_carre_diaGFSS[:,i])^2 for i in (k-1)*512+1:k*512]
end

#writedlm( "t_daGFSS_3_7.csv",  t_daGFSS)
#writedlm( "t_daGFSS_carre_3_7.csv",  t_daGFSS_carre)
#writedlm( "t_carre_daGFSS_3_7.csv",  t_carre_daGFSS)

#t_daGFSS = readdlm("t_daGFSS_3_7.csv")
#t_daGFSS_carre=readdlm("t_daGFSS_carre_3_7.csv")
#t_carre_daGFSS = readdlm("t_carre_daGFSS_3_7.csv")



#t_carre
init=300
fin=512
plot(t_carre_daGFSS[100,init:fin],xlabel="temps", ylabel="t_carre_aGFSS")
inter=fin-init # taille de la fenêtre étudiée
n_rupt = 400
Δ_rupt = 30
mini=round(minimum(t_carre_daGFSS[:,init:fin]))
maxi=round(maximum(t_carre_daGFSS[:,init:fin]))
pas=(maxi-mini)/size(T100)[1]
x1=(mini):pas:maxi
T1000=tab_threshold(mini,maxi,pas,nt) #creation d'un tableau de seuils

plot(t_carre_voisin[:,init:fin]',xlabel="temps", ylabel="t_carre_aGFSS", label="")
T20001=80*ones(250,512)

pdetect, retardtot, tdetecttot= performance_algo2(nt, sig1, t_carre_daGFSS, t_carre_voisin, ρ ,d ,v ,T1000, T20001)



#t_voisin_carre
plot(t_daGFSS_carre[100,init:fin],xlabel="temps", ylabel="t_daGFSS_carre")
mini2=round(minimum(t_daGFSS_carre[:,init:fin]))
maxi2=round(maximum(t_daGFSS_carre[:,init:fin]))
pas=(maxi2-mini2)/size(T1002)[1]
x2=(mini2):pas:maxi2
T10002=tab_threshold(mini2,maxi2,pas,nt)

plot(t_voisin_carre[:,init:fin]',xlabel="temps", ylabel="t_carre_aGFSS", label="")
T2000=1000*ones(250,512)

pdetect2, retardtot2, tdetecttot2= performance_algo2(nt, sig1, t_daGFSS_carre, t_voisin_carre, ρ ,d ,v ,T10002, T2000)


#t
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
for k in 1:nb
    sumtot_daGFSS[k,:]=sum(t_carre_diaGFSS[:,(k-1)*512+1:k*512],dims=1)
end
sumtot_daGFSS

plot(sumtot_daGFSS[100,init:fin],xlabel="temps", ylabel="t_sum_daGFSS")
mini4=round(minimum(sumtot_daGFSS[:,init:fin]))
maxi4=round(maximum(sumtot_daGFSS[:,init:fin]))
pas=(maxi4-mini4)/size(T1004)[1]
x4=(mini4):pas:maxi4
T10004=tab_threshold(mini4,maxi4,pas,nt) #creation d'un tableau de seuils


pdetect4, retardtot4, tdetecttot4= performance_algo2(nt, sig1, sumtot_daGFSS, t_carre_diaGFSS, ρ ,d ,v ,T10004, T20002)


plot(pfausse2, pdetect2, ratio=:equal, xlabel="pfa", ylabel="pd", label="voisin carre")
plot!(pfausse3, pdetect3, label="original")
plot!(pfausse, pdetect, label="carre voisin")
plot!(pfausse4, pdetect4, label="sum")

plot(pfausse2,retardtot2, xlabel="pfa", ylabel ="retard",  label="voisin carre")
plot!(pfausse3, retardtot3, label="original")
plot!(pfausse, retardtot, label="carre voisin")
plot!(pfausse4, retardtot4, label="su")