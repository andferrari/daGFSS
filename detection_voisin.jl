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

# ############################################
# signal on graph

nb=1000 #nb de signaux créés
nt=512 #temps d'étude
ρ = 1.0/0.3
init=250
fin=512

#lancer synth_gfss_parallele pour avoir les poles et les residues
φ, ψ = calcul_psi_phi(poles, residues)
c=addTerm

sig1 = gener_sigg(g, node_labels, 3.0,1, 3.0, σ2 = 7);
t_carre_diaGFSS, t_diaGFSS  = diaGFSS(sig1, L, ψ, φ, c; λ = 0.01, Λ=0.1)

t_voisin=-(t_diaGFSS + A*t_diaGFSS)
t_voisin_carre = t_voisin.^2
plot(t_voisin_carre',xlabel="temps", ylabel="t_voisin", label="")

t_carre_voisin=(t_carre_diaGFSS + A*t_carre_diaGFSS)
plot(t_carre_voisin',xlabel="temps", ylabel="t_voisin", label="")

x=(init+80):1:fin
