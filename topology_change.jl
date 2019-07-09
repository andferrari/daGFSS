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

gg = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
node_labelsg = Int.(label_propagation(gg, 10000)[1])
sig1 = gener_sigg(gg, node_labelsg, 1.0,1, 1, σ2 = 7);

p = 0.01
Lt = NormalizedLaplacian(gg)
Lt1 = eigvals(Array(Lt))
Ltest = laplacian_matrix(gg)
v1 = eigvals(Array(Ltest))
v1./diag(Ltest)
plot(v1./diag(Ltest))
plot!(Lt1)
tab_v = eigvals(Array((1-p)*Ltest - p*ones(250,250) + 250*p*Matrix(I,250,250)))
(1-p).*v1 .+250*p
plot(tab_v[2:250])
plot!((1-p)*v1[2:250].+250*p)

proba = 0:1:50
retard = zeros(length(proba))
compt = zeros(length(proba))
g = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
L1 = NormalizedLaplacian(g)
L1 = L1 - I
L1 = laplacian_matrix(g)
ones(250,250) + L1 -250*Matrix(I,250,250)
d, v = eigen(Array(L1));
node_labels = Int.(label_propagation(g, 10000)[1])


init=300
x = init:512
T3 = 20*ones(250,ntot)
T4 = 200000*ones(1,263)


L = zeros(250,250,ntot)
L[:,:,1] = L1

for it in 1:length(proba)
    println(it)

    retard[it], compt[it] = test_proba(g, L, d, v, ψ, φ, c, σi, L1, node_labels, T3, T4, proba[it])

end

plot(proba./1000, retard, seriestype=:scatter)

plot(proba./1000,compt, seriestype=:scatter)
