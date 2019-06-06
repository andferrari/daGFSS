using LightGraphs
using GraphIO
using EzXML
using LinearAlgebra
using Plots
using LaTeXStrings
pyplot()
using DelimitedFiles
using Statistics



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
ρ = 1.0

#lancer synth_gfss_parallele pour avoir les poles et les residues
φ, ψ = calcul_psi_phi(poles, residues)
c=addTerm

sig1=zeros(250, nt*nb)
t_diaGFSS=zeros(250, nt*nb)
t_daGFSS=zeros(nb,nt)
t_carre_diaGFSS=zeros(250, nt*nb)
t_carre_daGFSS=zeros(nb,nt)
t_daGFSS_carre=zeros(nb,nt)
sumtot_daGFSS=zeros(nb,nt)


#création des signaux
for k in 1:nb
    sig1[:,(k-1)*512+1:k*512]= gener_sigg(g, node_labels, 3.0, 2 , 3.0, σ2 = 7);
    t_carre_diaGFSS[:,(k-1)*512+1:k*512], t_diaGFSS[:,(k-1)*512+1:k*512] = diaGFSS(sig1[:,(k-1)*512+1:k*512], L, ψ, φ, c; λ = 0.01, Λ=0.1)
end

#writedlm( "sig2_3_7.csv",  sig1)
#writedlm( "t_diaGFSS2_7.csv",  t_diaGFSS)
#writedlm( "t_carre_diaGFSS2_7.csv",  t_carre_diaGFSS)

#sig1=readdlm("sig2_3_7.csv")
#t_diaGFSS = readdlm("t_diaGFSS2_7.csv")
#t_carre_diaGFSS = readdlm("t_carre_diaGFSS2_7.csv")


t_voisin=-(t_diaGFSS + A*t_diaGFSS)
t_voisin_carre = t_voisin.^2

t_carre_voisin=(t_carre_diaGFSS + A*t_carre_diaGFSS)

for k in 1:nb
    t_carre_daGFSS[k,:] = [norm(t_carre_voisin[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS_carre[k,:] = [norm(t_voisin_carre[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS[k,:] = [norm(t_carre_diaGFSS[:,i])^2 for i in (k-1)*512+1:k*512]
end

#writedlm( "t_daGFSS2_7.csv",  t_daGFSS)
#writedlm( "t_daGFSS_carre2_7.csv",  t_daGFSS_carre)
#writedlm( "t_carre_daGFSS2_7.csv",  t_carre_daGFSS)

#t_daGFSS = readdlm("t_daGFSS2_7.csv")
#t_daGFSS_carre=readdlm("t_daGFSS_carre2_7.csv")
#t_carre_daGFSS = readdlm("t_carre_daGFSS2_7.csv")

for k in 1:nb
    sumtot_daGFSS[k,:]=sum(t_carre_diaGFSS[:,(k-1)*512+1:k*512],dims=1)
end
sumtot_daGFSS


#t_carre
init=350
fin=512
plot(t_carre_daGFSS[100,init:fin],xlabel="temps", ylabel="t_carre_aGFSS")
inter=fin-init # taille de la fenêtre étudiée
n_rupt = 400
Δ_rupt = 30
pas=100
mini=round(minimum(t_carre_daGFSS[:,init:fin]))
maxi=round(maximum(t_carre_daGFSS[:,init:fin]))
x1=(mini+50):pas:maxi
T100=tab_threshold(mini,maxi,pas,nt) #creation d'un tableau de seuils


pfausse, tdetecttot= performance_algo(nt, t_carre_daGFSS, ρ ,d ,v ,T100)





#t_voisin_carre
plot(t_daGFSS_carre[100,init:fin],xlabel="temps", ylabel="t_daGFSS_carre")
pas=10000
mini2=round(minimum(t_daGFSS_carre[:,init:fin]))
maxi2=round(maximum(t_daGFSS_carre[:,init:fin]))
x2=(mini2):pas:maxi2
T1002=tab_threshold(mini2,maxi2,pas,nt)


pfausse2, tdetecttot2= performance_algo(nt, t_daGFSS_carre, ρ ,d ,v ,T1002,)

#t
plot(t_daGFSS[100,init:fin],xlabel="temps", ylabel="t_daGFSS")
pas=1
mini3=round(minimum(t_daGFSS[:,init:fin]))
maxi3=round(maximum(t_daGFSS[:,init:fin]))
x3=(mini3):pas:maxi3
T1003=tab_threshold(mini3,maxi3,pas,nt)


pfausse3, tdetecttot3= performance_algo(nt, t_daGFSS, ρ ,d ,v ,T1003)


#t_sum

plot(sumtot_daGFSS[100,init:fin],xlabel="temps", ylabel="t_sum_daGFSS")
pas=1.5
mini4=round(minimum(sumtot_daGFSS[:,init:fin]))
maxi4=round(maximum(sumtot_daGFSS[:,init:fin]))
x4=(mini4):pas:maxi4
T1004=tab_threshold(mini4,maxi4+6,pas,nt) #creation d'un tableau de seuils


pfausse4, tdetecttot4= performance_algo(nt, sumtot_daGFSS, ρ ,d ,v ,T1004)
