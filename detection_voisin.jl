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
L = NormalizedLaplacian(g)
A = adjacency_matrix(g)
d, v = eigen(Array(L));
λmax = maximum(d)
d
L = L - (λmax/2)I
node_labels = Int.(label_propagation(g, 10000)[1])



nt = 512 #temps d'étude
ρ = 1.0/0.3
init = 300
fin = 512
λ = 0.01
Λ = 0.1
σ2 = 7

poles = [0.7344635585485741 + 0.5293681924014867im 0.7344635585485741 - 0.5293681924014867im -0.7344635202651353 + 0.5293681751822871im -0.7344635202651353 - 0.5293681751822871im]
residues = [-0.05322205451164147 - 0.08206089787078102im -0.05322205451164147 + 0.08206089787078102im 0.16044187053514897 - 0.6853621323079733im 0.16044187053514897 + 0.6853621323079733im]
φ, ψ = calcul_psi_phi(poles, residues)
c = 0.6682305081233931

variancei = readdlm("donnees/variancei.csv")
σi = sqrt.(variancei)


sig1 = gener_sigg(g, node_labels, 1.0,1, 0.5, σ2 = 7);
t_carre_diaGFSS, t_diaGFSS  = diaGFSS(sig1, L, ψ, φ, c; λ = 0.01, Λ=0.1)
t_voisin = -(t_diaGFSS + A*t_diaGFSS)
t_voisin_carre = (t_voisin.^2)./σi
maxt = maximum(t_voisin_carre)

x = init:1:fin
T3 = 20*ones(250,512)
q = 90
p8 = same_plot(t_voisin_carre[q,x],T3[q,x])
plot(p8', xlabel="temps", ylabel="seuil noeud q")

t_daGFSS = [norm(t_voisin_carre[:,k])^2 for k in 1:512]
T4 = 50000*ones(1,263)
p8 = plot(t_daGFSS[x], xlabel="temps", ylabel="t_daGFSS", label="")
p9 = plot!(T4',xlabel="temps", ylabel="seuil norme")

t_change2 = detect_t_change(t_daGFSS[x] ,ρ , d, v,T4; λ = 0.01, Λ=0.1 ).+(init)
r2 = Int(floor(length(t_change2)*rand(1)[1]))
n_change2 = detect_n_change(sig1 , t_voisin_carre,ρ , d, v,t_change2[r2],T3; λ = 0.01, Λ=0.1 )
detect = detect_change(node_labels,sig1, t_daGFSS,t_voisin_carre, t_change2 ,ρ , d, v,T3, T4, 30; λ = 0.01, Λ=0.1)
