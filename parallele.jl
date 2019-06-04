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
d, v = eigen(Array(L));
λmax = maximum(d)

L= L - (λmax/2)I
node_labels = Int.(label_propagation(g, 10000)[1])

# ############################################
# signal on graph

nb=1000 #nb de signaux créés
nt=512 #temps d'étude
ρ = 1.0
init=250
fin=512

φ, ψ = calcul_psi_phi(poles, residues)
c=addTerm

sig1 = gener_sigg(g, node_labels, 3.0,1, 3.0, σ2 = 1);
t_diaGFSS = diaGFSS(sig1, L, ψ, φ, c; λ = 0.01, Λ=0.1)
plot(t_diaGFSS[:,:]', xlabel="temps", ylabel="t_diaGFSS", label="")

x=(init+30):1:fin
T3=5*ones(250,512)
p6 = plot(x,t_diaGFSS[:,x]', xlabel="temps", ylabel="t_diaGFSS",label="")
p7 = plot!(x,T3[:,x]', xlabel="temps", ylabel="seuil", label="")

q=240
p8=same_plot(t_diaGFSS[q,x],T3[q,x])
plot(p8', xlabel="temps", ylabel="seuil noeud q")

t_daGFSS = [norm(t_diaGFSS[:,k])^2 for k in 1:512]
T4=1000*ones(1,263)
p8 = plot(t_daGFSS[x], xlabel="temps", ylabel="t_daGFSS", label="")
p9 = plot!(T4',xlabel="temps", ylabel="seuil norme")

t_change2=detect_t_change(t_daGFSS[x] ,ρ , d, v,T4; λ = 0.01, Λ=0.1 ).+(init+30)
r2=Int(floor(length(t_change2)*rand(1)[1]))
n_change2=detect_n_change(sig1 , t_diaGFSS,ρ , d, v,t_change2[r2],T3; λ = 0.01, Λ=0.1 )
detect=detect_change(node_labels,sig1,t_daGFSS,t_diaGFSS, t_change2 ,ρ , d, v,T3, T4, 30; λ = 0.01, Λ=0.1)
