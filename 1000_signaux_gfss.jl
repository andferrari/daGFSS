using LightGraphs
using GraphIO
using GraphPlot
using EzXML
using ColorSchemes
using LinearAlgebra
using Plots

using LaTeXStrings
using PyCall
pyplot()

#using Compose
#using Cairo
using DelimitedFiles
#using Graphs

# load graph
include("gfss_func.jl")
include("performance_func.jl")
include("detection_func.jl")
include("signaux_func.jl")
g = loadgraph("MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
L = LightGraphs.laplacian_matrix(g)
d, v = eigen(Array(L));

node_labels = Int.(label_propagation(g, 10000)[1])
nodefillc = get(ColorSchemes.jet, rescale(node_labels));

# save layout
locs_x, locs_y = spring_layout(g)
writedlm( "locs.csv",  [locs_x, locs_y])
locs = readdlm("locs.csv")
plt = gplot(g, locs[1,:], locs[2,:], nodefillc=nodefillc, arrowlengthfrac=0)
# draw(PDF("../paper/figs/karate.png"), plt)

# ############################################
# signal on graph

nb=1000 #nb de signaux créés
nt=512 #temps d'étude


sig1=zeros(250, nt*nb)
sig2=zeros(250, nt*nb)
t_iaGFSS =zeros(250, nt*nb)
t_aGFSS = zeros(nb,nt)
t_iaGFSS2 =zeros(250, nt*nb)
t_aGFSS2 = zeros(nb,nt)
t_diaGFSS=zeros(250, nt*nb)
t_diaGFSS2=zeros(250, nt*nb)
t_daGFSS=zeros(nb,nt)
t_daGFSS2=zeros(nb,nt)



#création des signaux avec changement de comportement "sig1 sig1_15"
for k in 1:nb
    sig1[:,(k-1)*512+1:k*512]= gener_sigg(g, node_labels, 1.0, 1 , 1.0, σ2 = 3);
    t_iaGFSS[:,(k-1)*512+1:k*512] = iaGFSS(sig1[:,(k-1)*512+1:k*512], ρ, d, v; λ = 0.01, Λ = 0.1)
    t_aGFSS[k,:] = [norm(t_iaGFSS[:,i])^2 for i in (k-1)*512+1:k*512]
    t_diaGFSS[:,(k-1)*512+1:k*512] = diaGFSS(sig1[:,(k-1)*512+1:k*512], L, 0.01, 0.01, 0.1; λ = 0.01, Λ=0.1)
    t_daGFSS[k,:] = [norm(t_diaGFSS[:,i])^2 for i in (k-1)*512+1:k*512]
end

writedlm( "sig1_3.csv",  sig1)
writedlm( "t_aGFSS_3.csv",  t_aGFSS)
writedlm( "t_diaGFSS3.csv",  t_diaGFSS)
writedlm( "t_daGFSS3.csv",  t_daGFSS)

sig1 = readdlm("sig1_15.csv")
t_aGFSS = readdlm("t_aGFSS_15.csv")
t_diaGFSS = readdlm("t_diaGFSS15.csv")
t_daGFSS = readdlm("t_daGFSS15.csv")

init=250
fin=512
plot(t_aGFSS[800,init:fin],xlabel="temps", ylabel="t_aGFSS")
inter=fin-init # taille de la fenêtre étudiée
n_rupt = 400
Δ_rupt = 30
pas=0.1
mini=round(minimum(t_aGFSS[:,init:fin]))
maxi=round(maximum(t_aGFSS[:,init:fin]))
x=(mini+pas):pas:maxi
T1000=tab_threshold(mini,maxi,pas,nt) #creation d'un tableau de seuils


# verification des performances de la methode utilisée
pdetect, pfausse, retardtot, tdetecttot= performance_algo(nt, t_aGFSS,ρ[2] ,d ,v ,T1000)

plt=same_plot(pdetect,pfausse)
plt11=plot(x,plt', xlabel="seuil", ylabel="%", label=["pd" "pfa"])
plt12=plot(pfausse, pdetect, ratio=:equal, xlabel="pfa", ylabel="pd", label="CORE")
plt13=plot(x,moytot_interpics, xlabel="seuil", ylabel="moyenne du temps entre chaque pic ")
plt14=plot(x,retardtot, xlabel="seuil", ylabel ="retard")
plt15=plot(x,tdetecttot, xlabel="seuil", ylabel="temps 1ère FA ")


#avec le nouvel algorithme
plot(t_diaGFSS[1,1:512], xlabel="temps", ylabel="t_diaGFSS")
plot(t_daGFSS[800,init:fin], xlabel="temps", ylabel="t_aGFSS")
pas=0.00001
mini=round(minimum(t_daGFSS[:,init:fin])*10000)/10000+0.003
maxi=round(maximum(t_daGFSS[:,init:fin])*10000)/10000
x=(mini+pas):pas:maxi
T1000=tab_threshold(mini,maxi,pas,nt)

pdetect3, pfausse3, retardtot3, tdetecttot3= performance_algo(nt, t_daGFSS,ρ ,d ,v ,T1000)

p20=same_plot(pdetect3,pfausse3)
plt21=plot(x,p20', xlabel="seuil", ylabel="%", label=["pd" "pfa"])
plt22=plot(pfausse3, pdetect3, ratio=:equal, xlabel="pfa", ylabel="pd", label="CORE")
plt23=plot(x,moytot_interpics3, xlabel="seuil", ylabel="moyenne du temps entre chaque pic ")
plt24=plot(x,retardtot3, xlabel="seuil", ylabel ="retard")
plt25=plot(x,tdetecttot3, xlabel="seuil", ylabel="temps 1ère FA ")
#plot(t_interpic3[800,:,10], xlabel="temps", ylabel="temps entre les différents pics")
#plot(x,nbpics3[800,:], xlabel="seuil", ylabel="nb pics (1signal)")
#plot(x,moy_interpics3[800,:], xlabel="seuil", ylabel="moyenne du temps entre chaque pic ")

#comparaison des deux méthodes
plot(plt11,plt21)
savefig("pda_pf_3.png")
plot(plt12,plt22)
savefig("core_3.png")
plot(plt14,plt24)
savefig("retard_3.png")
plot(plt15,plt25)
savefig("1FA_3.png")


#création signaux sans changement de comportemegitnt "sig2" "sig2_1_20" "sig2_5_5" "sig2_rand"
for k in 1:nb
    sig2[:,(k-1)*512+1:k*512]= gener_sigg(g, node_labels, 1.0, 2, 1.0, σ2 =10);
    t_iaGFSS2[:,(k-1)*512+1:k*512] = iaGFSS(sig2[:,(k-1)*512+1:k*512], ρ, d, v; λ = 0.01, Λ = 0.1)
    t_aGFSS2[k,:] = [norm(t_iaGFSS2[:,i])^2 for i in (k-1)*512+1:k*512]
    t_diaGFSS2[:,(k-1)*512+1:k*512] = diaGFSS(sig2[:,(k-1)*512+1:k*512], L, 0.01, 0.01, 0.1; λ = 0.01, Λ=0.1)
    t_daGFSS2[k,:] = [norm(t_diaGFSS2[:,i])^2 for i in (k-1)*512+1:k*512]
end

writedlm( "sig2_rand.csv",  sig2)
writedlm( "t_aGFSS2rand.csv",  t_aGFSS2)
writedlm( "t_diaGFSS2_15.csv",  t_diaGFSS2)
writedlm( "t_daGFSS2_15.csv",  t_daGFSS2)

sig2 = readdlm("sig2_1_20.csv")
t_aGFSS2 = readdlm("t_aGFSS21_20.csv")
t_diaGFSS2 = readdlm("t_diaGFSS2_15.csv")
t_daGFSS2 = readdlm("t_daGFSS2_15.csv")


pas=0.1
mini=round(minimum(t_aGFSS2[:,init:fin]))
maxi=round(maximum(t_aGFSS2[:,init:fin]))-60
x=(mini+pas):pas:maxi
T1000=tab_threshold(mini,maxi,pas,nt)
pdetect2, pfausse2, moytot_interpics2, retardtot2, tdetecttot2= performance_algo(nt, t_aGFSS2,ρ ,d ,v ,T1000)



Plots.plot(x,pfausse2, xlabel="seuil", ylabel="pfa")

plot(t_aGFSS2[100,init:fin], xlabel="temps", ylabel="t_aGFSS")
plot(t_interpic2[800,:,100], xlabel="temps", ylabel="temps entre les différents pics")
plot(x,nbpics2[800,:], xlabel="seuil", ylabel="nb pics (1signal)")
plot(x,moy_interpics2[800,:],xlabel="seuil", ylabel="moyenne du temps entre chaque pic ")
plot(pfausse2,moytot_interpics2, xlabel="seuil", ylabel="moyenne du temps entre chaque pic ")
plot(x,tdetecttot2, xlabel="seuil", ylabel="temps 1ère FA ")
