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
include("signaux_func.jl")

g = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
L = NormalizedLaplacian(g)
A = adjacency_matrix(g)
λmax = maximum(d)
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


poles = [0.7344635585485741 + 0.5293681924014867im 0.7344635585485741 - 0.5293681924014867im -0.7344635202651353 + 0.5293681751822871im -0.7344635202651353 - 0.5293681751822871im]
residues = [-0.05322205451164147 - 0.08206089787078102im -0.05322205451164147 + 0.08206089787078102im 0.16044187053514897 - 0.6853621323079733im 0.16044187053514897 + 0.6853621323079733im]
φ, ψ = calcul_psi_phi(poles, residues)
c = 0.6682305081233931


R = real(variance_t(φ, ψ, c, L; σ2 =7))
R=diag(R).*Matrix{Float64}(I, 250, 250)
R2 = (I+A)*R*(I+A')
variance_ti = ones(length(d)).*diag(R2)





sig1_sc=zeros(250, nt*nb)
t_diaGFSS_sc=zeros(250, nt*nb)
t_carre_diaGFSS_sc=zeros(250, nt*nb)




#création des signaux
for k in 1:nb
    sig1_sc[:,(k-1)*512+1:k*512] = gener_sigg(g, node_labels, 1.0, 2 , 3.0, σ2 = 7);
    t_carre_diaGFSS_sc[:,(k-1)*512+1:k*512], t_diaGFSS_sc[:,(k-1)*512+1:k*512] = diaGFSS(sig1_sc[:,(k-1)*512+1:k*512], L, ψ, φ, c; λ = 0.01, Λ=0.1)
end

t_voisin_sc = -(t_diaGFSS_sc + A*t_diaGFSS_sc)

var_voisin_sc=zeros(250,nb)
for k in 1:nb
    println(k)
    var_voisin_sc[:,k] = var(t_voisin_sc[:,(k-1)*512+351:k*512], dims=2)
end

variance_ti2 = sum(var_voisin_sc,dims=2)./nb


plot(variancei)
plot!(variance_ti)

variancei = readdlm("donnees/variancei.csv")
σi = sqrt.(variancei)
