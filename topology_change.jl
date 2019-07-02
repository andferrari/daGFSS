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

poles = [0.9284586365913845 + 0.6691948262233165im 0.9284586365913845 - 0.6691948262233165im -0.9284586223955065 + 0.6691948202913867im -0.9284586223955065 - 0.6691948202913867im]
residues = [-0.09550841212039587 - 0.10204555134224505im -0.09550841212039587 + 0.10204555134224504im -0.023277450874456127 - 0.8479373939514138im  -0.023277450874456127 + 0.8479373939514138im]
φ, ψ = calcul_psi_phi(poles, residues)
c = 0.7257546820209667

variancei = readdlm("donnees/variancei.csv")
σi = sqrt.(variancei)

ntot = 512
ρ = 1.0/0.3
λ = 0.01
Λ = 0.1
σ2 = 7

proba = 1:2:20
retar = zeros(length(proba))
compt = zeros(length(proba))
gg = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
node_labelsg = Int.(label_propagation(gg, 10000)[1])
sig1 = gener_sigg(gg, node_labelsg, 1.0,1, 1, σ2 = 7);

for it in 1:length(proba)
    g = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
    L1 = NormalizedLaplacian(g)
    d, v = eigen(Array(L1));
    λmax = 2
    L1 = L1 - (λmax/2)I
    node_labels = Int.(label_propagation(g, 10000)[1])


    init=300
    x = init:512
    T3 = 10*ones(250,ntot)
    T4 = 150000*ones(1,263)

    L = zeros(250,250,ntot)
    L[:,:,1] = L1
    for k in 2:ntot
        g1 = copy(g)
        e = collect(edges(g1))
        for i in 1:250
            for j in 1:250
                if has_edge(g1, i, j)==false
                    rnd = round(1000*rand(1)[1])
                    if rnd<=proba[it]
                        add_edge!(g1,i,j)
                    end
                end
            end
        end
        L[:,:,k] = NormalizedLaplacian(g1) - (λmax/2)I
    end

    A = adjacency_matrix(g)
    t_square_diaGFSS, t_diaGFSS  = diaGFSS_topo(sig1, L, ψ, φ, c; λ = 0.01, Λ=0.1)


    t_neigh = -(t_diaGFSS + A*t_diaGFSS)
    t_neigh_square = (t_neigh.^2)./σi


    t_daGFSS = [norm(t_neigh_square[:,k])^2 for k in 1:ntot]

    println(maximum(t_neigh_square[:,x]))
    println(maximum(t_daGFSS[x]))

    t_change2 = detect_t_change(t_daGFSS[x] ,ρ , d, v,T4; λ = 0.01, Λ=0.1 ).+(init-1)
    println(t_change2)
    retar[it] = t_change2[1]-400
    detect = detect_change(node_labels,sig1, t_daGFSS,t_neigh_square, t_change2 ,ρ , d, v,T3, T4, 300; λ = 0.01, Λ=0.1)
    println(detect)
    compt[it] = detect[4]/1400
end

plot(proba./1000, retar)

plot(compt)




proba = 1:2:50
retard = zeros(length(proba))
compt = zeros(length(proba))

for it in 1:length(proba)
    println(it)
    g = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
    L1 = NormalizedLaplacian(g)
    d, v = eigen(Array(L1));
    λmax = 2
    L1 = L1 - (λmax/2)I
    node_labels = Int.(label_propagation(g, 10000)[1])


    init=300
    x = init:512
    T3 = 20*ones(250,ntot)
    T4 = 150000*ones(1,263)


    L = zeros(250,250,ntot)
    L[:,:,1] = L1
    retard[it], compt[it] = test_proba(g, L, ψ, φ, c, σi, L1, node_labels, T3, T4, proba[it])

end

plot(proba./1000, retard)

plot(proba./1000,compt)
