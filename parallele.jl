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

φ, ψ = calcul_psi_phi(poles, residues)
c=addTerm

sig1 = gener_sigg(g, node_labels, 3.0,1, 3.0, σ2 = 1);
t_diaGFSS = diaGFSS(sig1, L, ψ, φ, c; λ = 0.01, Λ=0.1)
plot(t_diaGFSS[:,:]', xlabel="temps", ylabel="t_diaGFSS", label="")
