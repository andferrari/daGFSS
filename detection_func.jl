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
