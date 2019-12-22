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

function gener_sigg(g,node_labels, α, choice, pertu;nt = 512, n_rupt = 400, Δ_rupt = 112, σ2 = 1)

    # choice 1 : change in cluster
    # choice 2 : no change
    # choixe 3 : change in 2 clusters
    # choice 4 : change in diffferent clusters

    sig1 = zeros(nv(g), nt)
    rnd=zeros(Int64,26,maximum(node_labels))

    if choice ==1
        for k in 1:length(node_labels)
            sig1[findall(in(k), node_labels),:] .= α*(k-1)
        end
        sig1[findall(in(4), node_labels), n_rupt:n_rupt+Δ_rupt] .+= pertu
    end
    if choice ==2
        for k in 1:length(node_labels)
            sig1[findall(in(k), node_labels),:] .= α*(k-1)
        end
    end


    if choice ==3
        for k in 1:length(node_labels)
            sig1[findall(in(k), node_labels),:] .= α*(k-1)
        end
        sig1[findall(in(5), node_labels), n_rupt-300:n_rupt+Δ_rupt-301] .+= pertu
        sig1[findall(in(5), node_labels), n_rupt:n_rupt+Δ_rupt-1] .+= pertu
        sig1[findall(in(3), node_labels), n_rupt-300:n_rupt+Δ_rupt-301] .+= pertu
        sig1[findall(in(3), node_labels), n_rupt:n_rupt+Δ_rupt-1] .+= pertu
    end

    if choice==4
        for k in 1:4
            index_nodes = findall(in(k), node_labels)
            sig1[index_nodes,:] .= α*(k-1)
            rnd[:,k]=rand(index_nodes,min(4,length(index_nodes)))
            sig1[rnd[:,k], n_rupt:n_rupt+Δ_rupt-1] .+= pertu
        end
    end
    sig1 += sqrt(σ2) * randn(nv(g), nt)
    return sig1
end

#detection du changement temporel avec la norme
function detect_t_change(t_aGFSS ,ρ , d, v,T; λ = 0.01, Λ=0.1 )
    j = 1
    for k in 1:length(t_aGFSS)
        if t_aGFSS[k] > T[1,k]
            j = j+1
        end
    end
    t_change = zeros(j-1)
    v = 1
    for k in 1:length(t_aGFSS)
        if t_aGFSS[k] > T[1,k]
            t_change[v] = k
            v = v+1
        end
    end
    return t_change
end

#detection du changement dans les noeuds
function detect_n_change(sig, t_iaGFSS, ρ , d, v, t, T; λ = 0.01, Λ=0.1)
    (nv, nt) = size(sig)
    j = 1
    t = Int(t)
    for k in 1:size(t_iaGFSS,1)
        if t_iaGFSS[k,t] > T[k,t]
            j = j+1
        end
    end
    n_change = zeros(j-1)
    v = 1
    for k in 1:size(t_iaGFSS,1)
        if t_iaGFSS[k,t] > T[k,t]
            n_change[v] = k
            v = v+1
        end
    end
    return n_change
end


#detection du changement dans les clusters
function detect_change(node_labels,sig,t_aGFSS, t_iaGFSS, t_change, ρ , d, v,T, T2, s; λ = 0.01, Λ=0.1)
    detect = zeros(maximum(node_labels))
    for k in 1:length(t_change)
        n_change = detect_n_change(sig ,t_iaGFSS ,ρ , d, v,t_change[k],T; λ = 0.01, Λ=0.1 )
        for j in 1:maximum(node_labels)
            for u in 1:length(n_change)
                for i in 1:length(findall(in(j), node_labels))
                    if Int(n_change[u]) == findall(in(j), node_labels)[i]
                        detect[j] += 1
                    end
                end
            end
        end
    end
    for k in 1:length(detect)
        if detect[k] >= s
            println("cluster $k has changed")
        end
    end
    return detect
end

#tableau de seuils
function tab_threshold(min, max, step, nt)
    nbvalue = Int(round((max-min)/step))
    T1000 = zeros(nbvalue,nt)
    for j in 1:nbvalue
        T1000[j,:] = j*step*ones(512)+min*ones(512)
    end
    return T1000
end

function compute_pfa(nt, t_aGFSS, ρ, d, v, T1000)
    inter = fin-init
    Texp = zeros(1,inter)
    nsigdetect = zeros(size(T1000)[1])
    tdetect = zeros(nb,size(T1000)[1])
    tdetectfin = zeros(size(T1000)[1])
    e2 = zeros(size(T1000)[1])


    for k in 1:nb
        for i in 1:size(T1000)[1]
            Texp[1,:] = T1000[i,init+1:fin]
            t_changeexp = detect_t_change(t_aGFSS[k,init+1:fin], ρ, d, v, Texp; λ = 0.01, Λ=0.1)
            if t_changeexp != zeros(0)
                    e2[i] = e2[i]+1 #false alarm
                    nsigdetect[i] += 1
                    tdetect[k,i] = t_changeexp[1] #first detection time
            end
        end
    end
    pfa = e2/nb

    for i in 1:size(T1000)[1] #mean for every threshold
        tdetectfin[i] = sum(tdetect[:,i])/nsigdetect[i]
    end
    return pfa, tdetectfin
end

function compute_pd(nt, sig1, t_aGFSS, t_diaGFSS, ρ, d, v, T1000)
    inter = fin-init
    Texp = zeros(1,inter)
    nsigdetect = zeros(size(T1000)[1])
    tdetect = zeros(nb,size(T1000)[1])
    tdetectfin = zeros(size(T1000)[1])
    e = zeros(size(T1000)[1])

    for k in 1:nb
        for i in 1:size(T1000)[1]
            Texp[1,:] = T1000[i,init+1:fin]
            t_changeexp = detect_t_change(t_aGFSS[k,init+1:fin],ρ , d, v,Texp; λ = 0.01, Λ=0.1 )
            if t_changeexp != zeros(0)
                e[i] = e[i]+1 #nb signals detected
            end
        end
    end
    pdetect = e/nb
    return pdetect
end

function compute_delay(nt, sig1, t_aGFSS, t_diaGFSS, ρ, d, v, T1000)
    inter = fin-init
    Texp = zeros(1,inter)
    delay = zeros(nb,size(T1000)[1])
    delaytot = zeros(size(T1000)[1])
    nsigdetect = zeros(size(T1000)[1])

    for k in 1:nb
        for i in 1:size(T1000)[1]
            Texp[1,:] = T1000[i,init+1:fin]
            t_changeexp = detect_t_change(t_aGFSS[k,n_rupt:fin],ρ , d, v,Texp; λ = 0.01, Λ=0.1 )
            if t_changeexp != zeros(0)
                delay[k,i] = t_changeexp[1]
                nsigdetect[i] += 1
            end
        end
    end

    for i in 1:size(T1000)[1]
        delaytot[i] = sum(delay[:,i])/nsigdetect[i]
    end
    return delaytot
end


function variance_t(φ, ψ, c, L; σ2 =7)
    Q=zeros(size(L))
    v=real(sum(φ))
    r=zeros(size(L))
    for i in 1:length(φ)
        for j in 1:length(φ)
            inver=(I-ψ[i]*conj(ψ[j])*(L^2))
            r=inv(Matrix(inver))
            Q+=φ[i]*conj(φ[j])*(σ2)*r
        end
    end
    Q+=c*(σ2)*(c+2*v)*I
    R = (λ/(2-λ) + Λ/(2-Λ) - 2*λ*Λ/(λ+Λ-λ*Λ))*Q
    return R
end

