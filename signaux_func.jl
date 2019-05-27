function gener_sigg(g,node_labels, α, choice, pertu;nt = 512, n_rupt = 400, Δ_rupt = 30, σ2 = 1)

    # choice 1 : changement du cluster 4
    # choice 2 : pas de changement
    # choixe 3 : changement de deux clusters
    # choice 4 : changementdes signaux sur plusieurs cluster en gardant la même puisssance du pic

    sig1 = zeros(nv(g), nt)
    rnd=zeros(Int64,4,maximum(node_labels))

    if choice ==1
        for k in 1:length(node_labels)
            sig1[findall(in(k), node_labels),:] .= α*(k-1)
        end
        sig1[findall(in(1), node_labels), n_rupt:n_rupt+Δ_rupt] .+= pertu

        sig1 += sqrt(σ2) * randn(nv(g), nt)
        return sig1
    end
    if choice ==2
        for k in 1:length(node_labels)
            sig1[findall(in(k), node_labels),:] .= α*(k-1)
        end

        sig1 += sqrt(σ2) * randn(nv(g), nt)
        return sig1
    end


    if choice ==3
        for k in 1:length(node_labels)
            sig1[findall(in(k), node_labels),:] .= α*(k-1)
        end
        sig1[findall(in(5), node_labels), n_rupt-300:n_rupt+Δ_rupt-301] .+= pertu
        sig1[findall(in(5), node_labels), n_rupt:n_rupt+Δ_rupt-1] .+= pertu
        sig1[findall(in(3), node_labels), n_rupt-300:n_rupt+Δ_rupt-301] .+= pertu
        sig1[findall(in(3), node_labels), n_rupt:n_rupt+Δ_rupt-1] .+= pertu

        sig1 += sqrt(σ2) * randn(nv(g), nt)
        return sig1
    end

    if choice==4
        for k in 1:4
            index_nodes = findall(in(k), node_labels)
            sig1[index_nodes,:] .= α*(k-1)
            rnd[:,k]=rand(index_nodes,min(4,length(index_nodes)))
            sig1[rnd[:,k], n_rupt:n_rupt+Δ_rupt-1] .+= pertu
        end
        sig1 += sqrt(σ2) * randn(nv(g), nt)
        return sig1, rnd
    end

    #p_sig1_h0 = gplot(g, locs[1,:], locs[2,:], nodefillc=get(ColorSchemes.jet, rescale(sig1[:,1])))
    #p_sig1_h1 = gplot(g, locs[1,:], locs[2,:], nodefillc=get(ColorSchemes.jet, rescale(sig1[:,n_rupt])))
    #p_sig2_h0 = gplot(g, locs[1,:], locs[2,:], nodefillc=get(ColorSchemes.jet, rescale(sig2[:,1])))
    #p_sig2_h1 = gplot(g, locs[1,:], locs[2,:], nodefillc=get(ColorSchemes.jet, rescale(sig2[:,n_rupt])))

    #mymap = get(ColorSchemes.jet, range(0,1,length=100))
    #cbar = colorbar(mymap, 0.0, 0.1, 1.0)
    #p = gridstack([p_sig1_h0 p_sig1_h1; p_sig2_h0 p_sig2_h1])


end
