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

# Constrained least squares with CVX

a_cvx = Variable(ka)
b_cvx = Variable(kb+1)

constraints = [minimum(va*a_cvx) >= -radius, maximum(va*a_cvx) <= radius]
problem = minimize( norm(vb*b_cvx - hμ_va*a_cvx - hμ), constraints)
solve!(problem, GurobiSolver())

pb_cvx = Poly(dropdims(b_cvx.value, dims=2))
pa_cvx = Poly([1; dropdims(a_cvx.value, dims=2)])
arma_cvx = [polyval(pb_cvx, μ)/polyval(pa_cvx, μ) for μ in μ]

maximum(abs.(roots(pa_cvx))) < λmax/2 ? error("unstable graph filter") :

# compute residue with matlab

a_final = reverse([1; dropdims(a_cvx.value, dims=2)])
b_final = reverse(dropdims(b_cvx.value, dims=2))

mat"[$residues, $poles, $addTerm] = residue($b_final, $a_final)"

arma_parallel = addTerm*ones(Complex, size(μ))
for (pole, residue) in zip(poles, residues)
    arma_parallel .+= residue./ (μ .- pole)
end
arma_parallel = real.(arma_parallel)

plot(μ .+ λmax/2, hμ, w=3, label=L"GFSS \ filter: \ h^\ast(\mu)", dpi=600)
plot!(μ .+ λmax/2, arma_parallel, w=3, label= "Order 4 parallel ARMA GFSS", xlabel = L"\mu", dpi=600)


# julia> poles
# 4-element Array{Complex{Float64},1}:
#   0.8748969843858525 + 0.5747054405842672im
#   0.8748969843858525 - 0.5747054405842672im
#  -0.7285516058955004 + 0.5746831597316937im
#  -0.7285516058955004 - 0.5746831597316937im
#
#  julia> residues
# 4-element Array{Complex{Float64},1}:
#  -0.11419392350542883 - 0.10630527603739182im
#  -0.11419392350542877 + 0.10630527603739177im
#   0.17839845033985413 - 0.6978803839459782im
#   0.17839845033985413 + 0.6978803839459782im
#
#  julia> addTerm
#   0.41950040212415846
