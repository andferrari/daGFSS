function rescale(x)
    if norm(x) == 0
        x_resc = x
    elseif minimum(x) == maximum(x)
        x_resc = normalize(x)
    else
        x_resc = (x .- minimum(x))./(maximum(x)-minimum(x))
    end
    x_resc
end


function filt_detec_graph(x, ρ, μ, Φ)
    g = zeros(size(x))
    for k in 2:length(μ)
        g += min(1.0, sqrt(ρ/μ[k]))*dot(v[:,k], x)*Φ[:,k]
    end
    return g
end

function GFSS(x, ρ, μ, Φ)

    (nv,nt) = size(x)
    t = zeros(nt)

    for k in 1:nt
        t[k] = norm(filt_detec_graph(x[:,k], ρ, μ, Φ))^2
    end
    return t

end

function iaGFSS(x, ρ, μ, Φ; λ = 0.001, Λ=0.01)

    (nv,nt) = size(x)
    t = zeros(nv, nt)

    g = filt_detec_graph(x[:,1], ρ, μ, Φ)
    g_slow = g
    g_fast = g
    t[:,1] = (g_fast - g_slow).^2

    for k in 2:nt
        g = filt_detec_graph(x[:,k], ρ, μ, Φ)

        g_slow = (1 - λ)*g_slow + λ*g
        g_fast = (1 - Λ)*g_fast + Λ*g

        t[:,k]  = (g_fast - g_slow).^2
    end
    return t

end

function diaGFSS(x, L, ψ, φ, c; λ = 0.01, Λ=0.1)

    (n_v,n_t) = size(x)
    t = zeros(n_v, n_t)
    t_voisin = zeros(n_v, n_t)
    g_slow = zeros(n_v)
    g_fast = zeros(n_v)
    v = zeros(length(ψ),n_v)
    for i in 1:length(ψ)
        v[i,:] = φ[i].re*x[:,1]
    end
    g = sum(v,dims=1)[1,:] + c*x[:,1]
    g_slow = g
    g_fast = g
    t[:,1] = (g_fast - g_slow).^2

    Rev = v
    Imv = zeros(length(ψ), n_v)
    Rev2 = zeros(length(ψ), n_v)
    Imv2 = zeros(length(ψ), n_v)

    for k in 2:n_t
        for j in 1:length(ψ)
            Rev2[j,:] = (ψ[j].re)*L*Rev[j,:] - (ψ[j].im)*L*Imv[j,:] + (φ[j].re)*x[:,k]
            Imv2[j,:] = (ψ[j].im)*L*Rev[j,:] + (ψ[j].re)*L*Imv[j,:] + (φ[j].im)*x[:,k]
        end

        g = sum(Rev2, dims=1)[1,:] + c*x[:,k]

        g_slow = (1 - λ)*g_slow + λ*g
        g_fast = (1 - Λ)*g_fast + Λ*g

        t_voisin[:,k]  = (g_fast - g_slow)
        t[:,k]  = (g_fast - g_slow).^2
        Rev = copy(Rev2)
        Imv = copy(Imv2)

    end
    return t, t_voisin

end

function calcul_psi_phi( p, r)
    ψ = ones(length(p))'./p
    φ = r.*ψ
    return φ, ψ
end



function colorbar(cmap, min, mid, max)
    n = length(cmap)
    compose(context(0.0, 0.0, 0.1, 1.0),
           (context(), Compose.rectangle([0.0], range(0.0, stop=1.0, length=n), [0.5], [1.0/n]), fill(reverse(cmap))),
           (context(), Compose.text([0.6], [0.0, 0.5, 1.0], ["$max", "$mid", "$min"], [hleft], [vtop, vcenter, vtop]))
    )
end;

function plot_tests(sig ,ρ , d, v; λ = 0.01, Λ=0.1)

    (nv, nt) = size(sig)

    t_GFSS = GFSS(sig,ρ , d, v)
    p1 = Plots.plot(t_GFSS, label="")
    xlabel!(L"t")
    ylabel!(L"t_{\rm{GFSS}}(y_t)")

    t_iaGFSS = iaGFSS(sig, ρ, d, v; λ = λ, Λ = Λ)
    t_aGFSS = [norm(t_iaGFSS[:,k])^2 for k in 1:nt]
    p2 = Plots.plot(t_aGFSS, label = "")
    xlabel!(L"t")
    ylabel!(L"t_{\rm{aGFSS}}(y_t)")

    p3 = Plots.plot(t_iaGFSS', label="")
    xlabel!(L"t")
    ylabel!(L"t_{\rm{iaGFSS}}(y_t)")

    p_sig = Plots.plot(p1, p2, p3, layout=(3,1))
end


function make_my_graph()
    gsp = pyimport("pygsp.graphs")
    G = gsp.Community(N=250, comm_density=0.3, seed=42)
    #G = gsp.Community(N=250, Nc=3, comm_sizes=[50, 120, 80], seed=42)
    g = SimpleDiGraph(G.A.todense())
    savegraph("MyGraph.graphml", g, GraphIO.GraphML.GraphMLFormat())
end

function NormalizedLaplacian(g)
    adjmat = LightGraphs.LinAlg.CombinatorialAdjacency(adjacency_matrix(g))
    I - Diagonal(adjmat.D.^(-1/2))*(adjmat.A)*Diagonal(adjmat.D.^(-1/2))
end
