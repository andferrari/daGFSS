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
init = 300
ρ = 1.0/0.3
λ = 0.01
Λ = 0.1
σ2 = 7
λmax = 2

gg = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
node_labelsg = Int.(label_propagation(gg, 10000)[1])
sig1 = gener_sigg(gg, node_labelsg, 1.0,1, 1, σ2 = 7);

p = 0.05
Ltest = NormalizedLaplacian(gg)
eigvals(Array(Ltest))
eigvals(Array((1-p)*Ltest + p*ones(250,250)))

proba = 0:0.1:50
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
    retard[it], compt[it] = test_proba(g, L, d, v, ψ, φ, c, σi, L1, node_labels, T3, T4, proba[it])

end

plot(proba[1:414]./1000, retard[1:414], seriestype=:scatter)

plot(proba[1:414]./1000,compt[1:414], seriestype=:scatter)
