using LightGraphs
using GraphIO
using EzXML
using LinearAlgebra
using Plots
using LaTeXStrings
pyplot()

using DelimitedFiles
using Statistics



include("gfss_func.jl")
include("performance_func.jl")
include("detection_func.jl")
include("signaux_func.jl")

g = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
L = NormalizedLaplacian(g)
A = adjacency_matrix(g)
d, v = eigen(Array(L));
λmax = 2
L = L - (λmax/2)I
node_labels = Int.(label_propagation(g, 10000)[1])



nt = 512 #temps d'étude
ρ = 1.0/0.3
init = 300
fin = 512
λ = 0.01
Λ = 0.1
σ2 = 7

poles = [0.9284586365913845 + 0.6691948262233165im 0.9284586365913845 - 0.6691948262233165im -0.9284586223955065 + 0.6691948202913867im -0.9284586223955065 - 0.6691948202913867im]
residues = [-0.09550841212039587 - 0.10204555134224505im -0.09550841212039587 + 0.10204555134224504im -0.023277450874456127 - 0.8479373939514138im  -0.023277450874456127 + 0.8479373939514138im]
φ, ψ = calcul_psi_phi(poles, residues)
c = 0.7257546820209667

variancei = readdlm("donnees/variancei.csv")
σi = sqrt.(variancei)


sig1 = gener_sigg(g, node_labels, 1.0,1, 3, σ2 = 3);
t_square_diaGFSS, t_diaGFSS  = diaGFSS(sig1, L, ψ, φ, c; λ = 0.01, Λ=0.1)
t_neigh = -(t_diaGFSS + A*t_diaGFSS)
t_neigh_square = (t_neigh.^2)./σi

x = init:fin
T3 = 55*ones(300,512)
q = 90
plot(t_neigh_square[:,x]', label="")
plot!(T3[q,x], xlabel="temps", ylabel="seuil noeud q")
savefig("GFSS_v2.pdf")

t_daGFSS = [norm(t_neigh_square[:,k])^2 for k in 1:512]
T4 = 120000*ones(1,213)
p8 = plot(t_daGFSS[x], xlabel="temps", ylabel="t_daGFSS", label="")
p9 = plot!(T4',xlabel="temps", ylabel="seuil norme")

t_change2 = detect_t_change(t_daGFSS[x] ,ρ , d, v,T4; λ = 0.01, Λ=0.1 ).+(init-1)
r2 = Int(floor(length(t_change2)*rand(1)[1]))
n_change2 = detect_n_change(sig1 , t_neigh_square,ρ , d, v,t_change2[r2],T3; λ = 0.01, Λ=0.1 )
detect = detect_change(node_labels,sig1, t_daGFSS,t_neigh_square, t_change2 ,ρ , d, v,T3, T4, 30; λ = 0.01, Λ=0.1)
