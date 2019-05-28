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
using Statistics
#using Graphs


# load graph
include("gfss_func.jl")
include("performance_func.jl")
include("detection_func.jl")
include("signaux_func.jl")
g = loadgraph("/Users/lverduci/Documents/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
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
ρ = 1.0



#choix des facteurs
#λ2=0.00001:0.00001:1/(1+50)
#Λ2=(1/(1+50))+0.00001:0.00001:0.9999
#λ3=zeros(length(Λ2))
#for k in 1:length(Λ2)
#    λ3[k]=solve(λ2,Λ2[k],nv(g))
#end
#λ3

#tabarg=argument_min(λ3,Λ2,nv(g))
#findmin(tabarg)
#Λ2[413]
#λ3[413]

#for k in 1:length(Λ2)
#    verif=log(Λ2[k]/λ3[k])/log((1-λ3[k])/(1-Λ2[k]))
#    println(verif)
#end

# ############################################
# affichage des clusters

σ2 = 3
init=250
fin=512
x=init:1:fin
sig1, rnd= gener_sigg(g, node_labels, 3.0,1, 3.0, σ2 = 3);
s1=sig1[findall(in(1), node_labels), x]
m1=mean(s1,dims=1)
plt1=plot(x,m1',ribbon=sqrt(3)*ones(512), color=RGB(1,136/255,5/255),label="")
s2=sig1[findall(in(2), node_labels), x]
m2=mean(s2,dims=1)
plot!(x,m2',ribbon=sqrt(3)*ones(512), color=RGB(0,196/255,1),label="")
s3=sig1[findall(in(3), node_labels), x]
m3=mean(s3,dims=1)
plot!(x,m3',ribbon=sqrt(3)*ones(512), color=RGB(211/255,179/255,176/255),label="")
s4=sig1[findall(in(4), node_labels), x]
m4=mean(s4,dims=1)
plot!(x,m4',ribbon=sqrt(3)*ones(512), color=RGB(223/255,137/255,1),label="")
s5=sig1[findall(in(5), node_labels), x]
m5=mean(s5,dims=1)
plot!(x,m5',ribbon=sqrt(3)*ones(512), color=RGB(1,85/255,132/255),label="")
s6=sig1[findall(in(6), node_labels), x]
m6=mean(s6,dims=1)
plot!(x,m6',ribbon=sqrt(3)*ones(512), color=RGB(0,189/255,148/255),label="")
s7=sig1[findall(in(7), node_labels),x]
m7=mean(s7,dims=1)
plot!(x,m7',ribbon=sqrt(3)*ones(512), color=RGB(115/255,192/255,0),label="")
s8=sig1[findall(in(8), node_labels), x]
m8=mean(s8,dims=1)
plot!(x,m8',ribbon=sqrt(3)*ones(512), color=RGB(76/255,70/255,62/255), label="")
plt9=plot!(xlab=L"time",ylab=L"signals \ on \ vertices")
vline!([300], label="")
#plot(plt9,p3, layout=grid(2,1,heights=[0.7,0.3]), dpi=300)


###########################################
#detection 1er algorithme

t_iaGFSS = iaGFSS(sig1, ρ, d, v; λ = 0.01, Λ = 0.1)
T=get_threshold(t_iaGFSS,0.003,2,10)
p1 = plot(t_iaGFSS'[:,1:10], xlabel="temps", ylabel="t_iaGFSS")
p2 = plot(T'[:,1:10], xlabel="temps", ylabel="seuil noeud")
plot(p1, p2, layout=(2,1))
#detection noeud q
q=6
plot(t_iaGFSS[q,:])
plot!(T[q,:],xlabel="temps", ylabel="seuil pour le noeud q")


#detection norme

t_aGFSS = [norm(t_iaGFSS[:,k])^2 for k in 1:512]
T2=get_threshold(t_aGFSS,0.018,1.8,1.5)
p3 = plot(x,t_aGFSS[x],label="", xlab = L"time", ylab=L"t_{aGFSS}")
plot(t_aGFSS)
plot!(T2',xlabel="temps", ylabel="seuil norme")

#verification noeud

rnd
t_change=detect_t_change(t_aGFSS ,ρ , d, v,T2; λ = 0.01, Λ=0.1 )
r=Int(floor(length(t_change)*rand(1)[1]))
n_change=detect_n_change(sig1 , t_iaGFSS,ρ , d, v,t_change[r],T; λ = 0.01, Λ=0.1 )
detect=detect_change(node_labels, sig1,t_aGFSS, t_iaGFSS, t_change ,ρ , d, v,T, T2, 50; λ = 0.01, Λ=0.1)

# ############################################

# detection algorithme 2
φ[103]
ψ[200]
c=optim_c(a, b, ρ, φ[103], ψ[200])
t_diaGFSS = diaGFSS(sig1, L, ψ[253], φ[132], c; λ = 0.01, Λ=0.1)
plot(t_diaGFSS', xlabel="temps", ylabel="t_diaGFSS", label="")
savefig("t_diaGFSS_sig1-4.pdf")
# T3=get_threshold(t_diaGFSS,0.006,1.8,8) adaptatif
x=(init+30):1:fin
T3=5*ones(250,512)
p6 = plot(x,t_diaGFSS[:,x]', xlabel="temps", ylabel="t_diaGFSS",label="")
p7 = plot!(x,T3[:,x]', xlabel="temps", ylabel="seuil", label="")
savefig("diaGFSS_σ_7_sig1-4.pdf")

q=170
p8=same_plot(t_diaGFSS[q,x],T3[q,x])
plot(p8', xlabel="temps", ylabel="seuil noeud q")
savefig("diaGFSS_σ_7_1noeud_sig1-4.pdf")




t_daGFSS = [norm(t_diaGFSS[:,k])^2 for k in 1:512]
#T4=get_threshold(t_daGFSS[init:fin],0.018,1.8,10)
T4=1000*ones(1,263)
p8 = plot(t_daGFSS[x], xlabel="temps", ylabel="t_daGFSS", label="")
p9 = plot!(T4',xlabel="temps", ylabel="seuil norme")
savefig("daGFSS_σ_3_sig1-4.pdf")

rnd
t_change2=detect_t_change(t_daGFSS[x] ,ρ , d, v,T4; λ = 0.01, Λ=0.1 ).+(init+30)
r2=Int(floor(length(t_change2)*rand(1)[1]))
n_change2=detect_n_change(sig1 , t_diaGFSS,ρ , d, v,t_change2[r2],T3; λ = 0.01, Λ=0.1 )
detect=detect_change(node_labels,sig1,t_daGFSS,t_diaGFSS, t_change2 ,ρ , d, v,T3, T4, 30; λ = 0.01, Λ=0.1)


#savefig(plt, "res.pdf")

#for k in 1:512
#    p_sig = gplot(g, locs[1,:], locs[2,:], nodefillc=get(ColorSchemes.jet, rescale(t_iaGFSS[:,k])))
#    draw(PNG("/Users/lverduci/Documents/graph_$k.png"), p_sig)
#end

#; ffmpeg -i Animate/graph_%d.png -vf scale=320:240 -pix_fmt yuv420p Animate/out.mp4
#; open Animate/out.mp4
