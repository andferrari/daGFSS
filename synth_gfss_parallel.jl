using LightGraphs
using GraphIO
using EzXML
using LinearAlgebra

using Plots
pyplot()

using MATLAB

using Convex
using Gurobi

using Polynomials

function NormalizedLaplacian(g)
    adjmat = LightGraphs.LinAlg.CombinatorialAdjacency(adjacency_matrix(g))
    I - Diagonal(adjmat.D.^(-1/2))*(adjmat.A)*Diagonal(adjmat.D.^(-1/2))
end

g = loadgraph("data/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())

Ln = NormalizedLaplacian(g)
d, v = eigen(Array(Ln));
λmax = 2

# center normalized Laplacian => |eigen(Lc)| < 1
Lc = Ln - (λmax/2)*I

# filter response

ρ = 1/0.3
h(μ) = min.(ρ, sqrt.(ρ./(μ.+λmax/2)))


n_μ = 400
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

minimum(abs.(roots(pa_cvx))) < λmax/2 ? error("unstable graph filter") : nothing

# compute residue with matlab

a_final = reverse([1; dropdims(a_cvx.value, dims=2)])
b_final = reverse(dropdims(b_cvx.value, dims=2))

mat"[$residues, $poles, $addTerm] = residue($b_final, $a_final)"

arma_parallel = addTerm*ones(Complex, size(μ))
for (pole, residue) in zip(poles, residues)
    arma_parallel .+= residue./ (μ .- pole)
end
arma_parallel = real.(arma_parallel)

Plots.reset_defaults()
Plots.scalefontsizes(1.5)
plot(μ .+ λmax/2, hμ, w=3, label="GFSS filter \$ h^\\ast(\\mu) \$", dpi=600)
plot!(μ .+ λmax/2, arma_parallel, w=3, label= "ARMA\$_4\$ GFSS", xlabel = "\$\\mu\$", dpi=600)

# sol_ls =  hcat(vb, - hμ_va)\hμ
# pb_ls = Poly(sol_ls[1:kb+1])
# pa_ls = Poly([1; sol_ls[kb+2:end]])
# arma_ls = [polyval(pb_ls, μ)/polyval(pa_ls, μ) for μ in μ]
# plot!(μ .+ λmax/2, arma_ls, w=3, label= "ARMA\$_4\$  unconstraint", xlabel = "\\mu", dpi=600)

# savefig("../paper/figs/approx_filt.png")
