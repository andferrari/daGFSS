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
include("performances_func.jl")
include("signals_func.jl")

g = loadgraph("data/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
L = NormalizedLaplacian(g)
A = adjacency_matrix(g)
λmax = 2.15
L = L - (λmax/2)I
d, v = eigen(Array(L));
node_labels = Int.(label_propagation(g, 10000)[1])


nt = 512
ρ = 1.0/0.3
init = 300
fin = 512
λ = 0.01
Λ = 0.1
σ2 = 7
nb = 100


poles = [0.9284586365913845 + 0.6691948262233165im 0.9284586365913845 - 0.6691948262233165im -0.9284586223955065 + 0.6691948202913867im -0.9284586223955065 - 0.6691948202913867im]
residues = [-0.09550841212039587 - 0.10204555134224505im -0.09550841212039587 + 0.10204555134224504im -0.023277450874456127 - 0.8479373939514138im  -0.023277450874456127 + 0.8479373939514138im]
φ, ψ = calcul_psi_phi(poles, residues)
c = 0.7257546820209667

R = real(variance_t(φ, ψ, c, L; σ2 =7))
R2 = (I+A)*R*(I+A')
variance_ti = ones(length(d)).*diag(R2)




sig1_sc = zeros(250, nt*nb)
t_diaGFSS_sc = zeros(250, nt*nb)
t_square_diaGFSS_sc = zeros(250, nt*nb)

for k in 1:nb
    sig1_sc[:,(k-1)*512+1:k*512] = gener_sigg(g, node_labels, 1.0, 2 , 3.0, σ2 = 7);
    t_square_diaGFSS_sc[:,(k-1)*512+1:k*512], t_diaGFSS_sc[:,(k-1)*512+1:k*512] = diaGFSS(sig1_sc[:,(k-1)*512+1:k*512], L, ψ, φ, c; λ = 0.01, Λ=0.1)
end

t_neigh_sc = -(t_diaGFSS_sc + A*t_diaGFSS_sc)

var_neigh_sc = zeros(250,nb)
for k in 1:nb
    println(k)
    var_neigh_sc[:,k] = var(t_neigh_sc[:,(k-1)*512+351:k*512], dims=2)
end

variance_ti2 = sum(var_neigh_sc,dims=2)./nb


plot(variance_ti)
plot!(variance_ti2)
