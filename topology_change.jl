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
sig1 = gener_sigg(g, node_labels, 1.0,1, 1, nt=512, n_rupt = 400, σ2 = 7);
ntot = 512
ρ = 1.0/0.3
λ = 0.01
Λ = 0.1
σ2 = 7
proba=1:18
retar=zeros(length(proba))
compt=zeros(length(proba))

for it in 1:length(proba)
    g = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
    L1 = NormalizedLaplacian(g)
    d, v = eigen(Array(L1));
    λmax = maximum(d)
    L1 = L1 - (λmax/2)I
    node_labels = Int.(label_propagation(g, 10000)[1])

    L = zeros(250,250,ntot)
    L[:,:,1] = L1

    A = adjacency_matrix(g)
    println(it)
    global A
    for k in 2:ntot
        for i in 1:size(L)[1]
            for j in 1:size(L)[2]
                rnd = round(10000*rand(1)[1])
                if rnd<=proba[it] && k%2==0 && A[i,j]==1
                    add_edge!(g,i,j)
                elseif rnd<=proba[it] && k%2==1 && A[i,j]==1
                    rem_edge!(g,i,j)
                end
            end
        end
        L1 = NormalizedLaplacian(g)
        d2 = real(eigvals(Array(L1)))
        L1 = L1 - (λmax/2)I
        L[:,:,k] = L1
    end

    A = adjacency_matrix(g)
    t_square_diaGFSS, t_diaGFSS  = diaGFSS_topo(sig1, L, ψ, φ, c; λ = 0.01, Λ=0.1)


    t_neigh = -(t_diaGFSS + A*t_diaGFSS)
    t_neigh_square = (t_neigh.^2)./σi

    init=301
    x = init:512
    T3 = 8*ones(250,ntot)
    t_daGFSS = [norm(t_neigh_square[:,k])^2 for k in 1:ntot]

    t_change2 = detect_t_change(t_daGFSS[x] ,ρ , d, v,T4; λ = 0.01, Λ=0.1 ).+(init)
    retar[it] = t_change2[1]-400
    detect = detect_change(node_labels,sig1, t_daGFSS,t_neigh_square, t_change2 ,ρ , d, v,T3, T4, 30; λ = 0.01, Λ=0.1)
    compt[it] = detect[4]/11400
end
