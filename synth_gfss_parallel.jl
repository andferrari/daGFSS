using Plots
pyplot()

using Convex
using Gurobi
using Polynomials


"""
     r, p, k = my_residue(b, a)

     naïve partial fraction expansion of
     Σ_{k=0}^M b[k]x^k / Σ_{k=0}^M a[k]x^k)
"""
function my_residue(b, a)

    p = roots(Poly(a))
    c = b[end]/a[end]

    poly_b = Poly(b)
    r = zeros(Complex, length(p))
    for n in 1:length(p)
        p_m = copy(p)
        deleteat!(p_m, n)
        r[n] = polyval(poly_b, p[n])/polyval(poly(p_m), p[n])
    end

    return r/a[end], p, c
end


# filter response

λmax = 2.0
ρ = 1/0.3
h(μ) = min.(ρ, sqrt.(ρ./(μ.+λmax/2)))


n_μ = 400
μ = range(-λmax/2, stop = λmax/2, length = n_μ)
hμ = h(μ)

# approximation

k_ord = 4  # filter order   \sum_{i=0}^k_ord b_i μ^i/1 + \sum_{i=1}^k_ord a_i μ^i
radius = 0.1

vb = zeros(n_μ, k_ord +1)
for k in 1:k_ord +1
    vb[:,k] = μ.^(k-1)
end

va = zeros(n_μ, k_ord )
for k in 1:k_ord
    va[:,k] = μ.^k
end
hμ_va = Diagonal(hμ)*va

# Constrained least squares with CVX

a_cvx = Variable(k_ord )
b_cvx = Variable(k_ord +1)

constraints = [minimum(va*a_cvx) >= -radius, maximum(va*a_cvx) <= radius]
problem = minimize( norm(vb*b_cvx - hμ_va*a_cvx - hμ), constraints)
solve!(problem, GurobiSolver())

pb_cvx = Poly(dropdims(b_cvx.value, dims=2))
pa_cvx = Poly([1; dropdims(a_cvx.value, dims=2)])
arma_cvx = [polyval(pb_cvx, μ)/polyval(pa_cvx, μ) for μ in μ]

minimum(abs.(roots(pa_cvx))) < λmax/2 ? error("unstable graph filter") : nothing

# partial fraction expansion

a_final = [1; dropdims(a_cvx.value, dims=2)]
b_final = dropdims(b_cvx.value, dims=2)

residues, poles, addTerm = my_residue(b_final, a_final)

arma_parallel = addTerm*ones(Complex, size(μ))
for (pole, residue) in zip(poles, residues)
    arma_parallel .+= residue./ (μ .- pole)
end
arma_parallel = real.(arma_parallel)

Plots.reset_defaults()
Plots.scalefontsizes(1.5)
plot(μ .+ λmax/2, hμ, w=3, label="GFSS filter \$ h^\\ast(\\mu) \$", dpi=600)
plot!(μ .+ λmax/2, arma_parallel, w=3, label= "ARMA\$_4\$ GFSS", xlabel = "\$\\mu\$", dpi=600)
