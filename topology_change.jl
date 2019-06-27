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
include("gfss_topo_func.jl")

poles = [0.7344635585485741 + 0.5293681924014867im 0.7344635585485741 - 0.5293681924014867im -0.7344635202651353 + 0.5293681751822871im -0.7344635202651353 - 0.5293681751822871im]
residues = [-0.05322205451164147 - 0.08206089787078102im -0.05322205451164147 + 0.08206089787078102im 0.16044187053514897 - 0.6853621323079733im 0.16044187053514897 + 0.6853621323079733im]
φ, ψ = calcul_psi_phi(poles, residues)
c = 0.6682305081233931

variancei = readdlm("donnees/variancei.csv")
σi = sqrt.(variancei)

ntot = 512
ρ = 1.0/0.3
λ = 0.01
Λ = 0.1
σ2 = 7


g = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
L1 = NormalizedLaplacian(g)
d, v = eigen(Array(L1));
λmax = 2
L1 = L1 - (λmax/2)I
node_labels = Int.(label_propagation(g, 10000)[1])
sig1 = gener_sigg(g, node_labels, 1.0,1, 1, nt=512, n_rupt = 400, σ2 = 7);


init=301
x = init:512
T3 = 50*ones(250,ntot)
T4 = 200000*ones(1,263)

proba=10
L = zeros(250,250,ntot)
L[:,:,1] = L1
for k in 2:ntot
    println(k)
    global g
    g1 = copy(g)
    e = collect(edges(g1))
    L2= zeros(250,250)
    for i in 1:length(e)
        rnd = round(1000*rand(1)[1])
        if rnd<=proba
            rem_edge!(g1,e[i])
        end
    end
    L[:,:,k] = NormalizedLaplacian(g1) - (λmax/2)I
end
L
t_square_diaGFSS, t_diaGFSS  = diaGFSS_topo(sig1, L, ψ, φ, c; λ = 0.01, Λ=0.1)

A = adjacency_matrix(g)
t_neigh = -(t_diaGFSS + A*t_diaGFSS)
t_neigh_square = (t_neigh.^2)./σi


t_daGFSS = [norm(t_neigh_square[:,k])^2 for k in 1:ntot]
plot!(t_daGFSS[x])

maxn = maximum(t_neigh_square[:,x])
maxnorme = maximum(t_daGFSS[x])

t_change2 = detect_t_change(t_daGFSS[x] ,ρ , d, v,T4; λ = 0.01, Λ=0.1 ).+(init-1)
retar = t_change2[1]-400
detect = detect_change(node_labels,sig1, t_daGFSS,t_neigh_square, t_change2 ,ρ , d, v,T3, T4, 300; λ = 0.01, Λ=0.1)
compt = detect[4]/1400
