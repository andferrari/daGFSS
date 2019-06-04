using LightGraphs
using GraphIO
using EzXML
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
include("signaux_func.jl")
g = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
L = LightGraphs.laplacian_matrix(g)
d, v = eigen(Array(L));

node_labels = Int.(label_propagation(g, 10000)[1])


# ############################################
# signal on graph

nb=1000 #nb de signaux créés
nt=512 #temps d'étude


sig2=zeros(250, nt*nb)
t_iaGFSS2 =zeros(250, nt*nb)
t_aGFSS2 = zeros(nb,nt)
t_diaGFSS2=zeros(250, nt*nb)
t_daGFSS2=zeros(nb,nt)


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


avec le nouvel algorithme
plot(t_diaGFSS[1,1:512], xlabel="temps", ylabel="t_diaGFSS")
plot(t_daGFSS[800,init:fin], xlabel="temps", ylabel="t_aGFSS")
pas=0.00001
mini=round(minimum(t_daGFSS[:,init:fin])*10000)/10000+0.003
maxi=round(maximum(t_daGFSS[:,init:fin])*10000)/10000
x=(mini+pas):pas:maxi
T1000=tab_threshold(mini,maxi,pas,nt)

#avec le nouvel algo
pdetect4, pfausse4, retardtot4, tdetecttot4= performance_algo(nt, t_daGFSS2,ρ ,d ,v ,T1000)

p20=same_plot(pdetect4,pfausse4)
plt31=plot(x,p20', xlabel="seuil", ylabel="%", label=["pd" "pfa"])
plt32=plot(pfausse4, pdetect4, ratio=:equal, xlabel="pfa", ylabel="pd", label="CORE")
plt34=plot(x,retardtot4, xlabel="seuil", ylabel ="retard")
plt35=plot(x,tdetecttot4, xlabel="seuil", ylabel="temps 1ère FA ")


#comparaison des deux méthodes
plot(plt11,plt31)
plot(plt12,plt32)
plot(plt14,plt34)
plot(plt15,plt35)
