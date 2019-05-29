using LightGraphs
using GraphIO
using EzXML
using LinearAlgebra

using Plots
using LaTeXStrings
pyplot()

using MATLAB

using Convex
using Gurobi

using Polynomials

function NormalizedLaplacian(g)
    adjmat = LightGraphs.LinAlg.CombinatorialAdjacency(adjacency_matrix(g))
    Ln = I - Diagonal(adjmat.D.^(-1/2))*(adjmat.A)*Diagonal(adjmat.D.^(-1/2))
end

g = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())

Ln = NormalizedLaplacian(g)
d, v = eigen(Array(Ln));
λmax = maximum(d)#d[nv(g)]

# center normalized Laplacian => |eigen(Lc)| < 1
Lc = Ln - (λmax/2)*I

# filter response

h(μ) = min.(ρ, sqrt.(ρ./(μ.+λmax/2)))
ρ = 3

n_μ = 300
μ = range(-λmax/2, stop = λmax/2, length = n_μ)
hμ = h(μ)

# approximation

kb = 4  # num. order   \sum_{i=0}^kb b_i μ^i
ka = 4  # denum. order 1 + \sum_{i=1}^ka a_i μ^i
radius = 0.1

vb = zeros(n_μ, kb+1)
for k in 1:kb+1
    vb[:,k] = μ.^(k-1)
end

va = zeros(n_μ, ka)
for k in 1:ka
    va[:,k] = μ.^k
end
hμ_va = Diagonal(hμ)*va

# 1. Simple Least squares => unstable
#
# H = hcat(vb, -hμ_va)
# ls = H\hμ
#
# pb_ls = Poly(ls[1:kb+1])
# pa_ls = Poly([1; ls[kb+2:end]])
#
# arma_ls = [polyval(pb_ls, μ)/polyval(pa_ls, μ) for μ in μ]
# println("modules poles LS : ",abs.(roots(pa_ls)), ", λmax/2 : ", λmax/2)
# plot!(μ , arma_ls, w=2, label="ARMA LS")

# 2. Constrained least squares with CVX

a_cvx = Variable(ka)
b_cvx = Variable(kb+1)

constraints = [minimum(va*a_cvx) >= -radius, maximum(va*a_cvx) <= radius]
problem = minimize( norm(vb*b_cvx - hμ_va*a_cvx - hμ), constraints)
solve!(problem, GurobiSolver())

pb_cvx = Poly(dropdims(b_cvx.value, dims=2))
pa_cvx = Poly([1; dropdims(a_cvx.value, dims=2)])
arma_cvx = [polyval(pb_cvx, μ)/polyval(pa_cvx, μ) for μ in μ]

# println("modules poles CVX : ",abs.(roots(pa_cvx)), ", λmax/2 : ", λmax/2)

maximum(abs.(roots(pa_cvx))) < λmax/2 ? error("unstable graph filter") :

# 3. Refined => useless

# a_refined = dropdims(a_cvx.value, dims=2)
# vb_ref = zeros(n_μ, ka+1)
# for k in 1:ka+1
#     vb_ref[:,k] = μ.^(k-1)
# end
# b_refined = vb_ref\(hμ_va*a_refined + hμ)
#
# pb_refined = Poly(b_refined)
# pa_refined = pa_cvx
# arma_refined = [polyval(pb_refined, μ)/polyval(pa_refined, μ) for μ in μ]
# plot!(μ , arma_refined, label="ARMA refined")

# residue with matlab

a_final = reverse([1; dropdims(a_cvx.value, dims=2)])
b_final = reverse(dropdims(b_cvx.value, dims=2))

mat"[$residues, $poles, $addTerm] = residue($b_final, $a_final)"

arma_parallel = addTerm*ones(Complex, size(μ))
for (pole, residue) in zip(poles, residues)
    arma_parallel .+= residue./ (μ .- pole)
end
arma_parallel = real.(arma_parallel)

plot(μ, hμ, label=L"h^\ast(\mu)")
plot!(μ , arma_parallel, w=2, label="order 4 parallel ARMA")
