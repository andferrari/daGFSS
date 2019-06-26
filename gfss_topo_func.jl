function diaGFSS_topo(x, L, ψ, φ, c; λ = 0.01, Λ=0.1)

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
            Rev2[j,:] = (ψ[j].re)*L[:,:,k]*Rev[j,:] - (ψ[j].im)*L[:,:,k]*Imv[j,:] + (φ[j].re)*x[:,k]
            Imv2[j,:] = (ψ[j].im)*L[:,:,k]*Rev[j,:] + (ψ[j].re)*L[:,:,k]*Imv[j,:] + (φ[j].im)*x[:,k]
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
