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
include("detection_func.jl")
include("signaux_func.jl")

g = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
#L = LightGraphs.laplacian_matrix(g)
L = NormalizedLaplacian(g)
A=adjacency_matrix(g)
d, v = eigen(Array(L));
λmax = maximum(d)
d
L= L - (λmax/2)I
node_labels = Int.(label_propagation(g, 10000)[1])


# ############################################
# signal on graph

nt=512 #temps d'étude
ρ = 1.0/0.3
init=250
fin=512
λ=0.01
Λ=0.1
σ=7
nb=10000

#lancer synth_gfss_parallele pour avoir les poles et les residues
poles=[0.7344635585485741 + 0.5293681924014867im 0.7344635585485741 - 0.5293681924014867im -0.7344635202651353 + 0.5293681751822871im -0.7344635202651353 - 0.5293681751822871im]
residues=[-0.05322205451164147 - 0.08206089787078102im -0.05322205451164147 + 0.08206089787078102im 0.16044187053514897 - 0.6853621323079733im 0.16044187053514897 + 0.6853621323079733im]
φ, ψ = calcul_psi_phi(poles, residues)
c=0.6682305081233931

R = variance_t(φ, ψ, c, L; σ2 =7)
R2= (I+A)*R*(I+A')
variance_ti=ones(length(d)).*diag(R2)





sig1_sc=zeros(250, nt*nb)
t_diaGFSS_sc=zeros(250, nt*nb)
t_carre_diaGFSS_sc=zeros(250, nt*nb)




#création des signaux
for k in 1:nb
    sig1_sc[:,(k-1)*512+1:k*512]= gener_sigg(g, node_labels, 1.0, 2 , 3.0, σ2 = 7);
    t_carre_diaGFSS_sc[:,(k-1)*512+1:k*512], t_diaGFSS_sc[:,(k-1)*512+1:k*512] = diaGFSS(sig1_sc[:,(k-1)*512+1:k*512], L, ψ, φ, c; λ = 0.01, Λ=0.1)
end

t_voisin_sc=-(t_diaGFSS_sc + A*t_diaGFSS_sc)

var_voisin_sc=zeros(250,nb)
for k in 1:nb
    println(k)
    var_voisin_sc[:,k]=var(t_voisin_sc[:,(k-1)*512+351:k*512], dims=2)
end

variance_ti2=sum(var_voisin_sc,dims=2)./nb


plot(variance_ti)
plot!(variance_ti2)
