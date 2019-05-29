#detection du changement temporel avec la norme
function detect_t_change(t_aGFSS ,ρ , d, v,T; λ = 0.01, Λ=0.1 )

    j=1
    for k in 1:length(t_aGFSS)
        if t_aGFSS[k]>T[1,k]
            j=j+1
        end
    end
    t_change=zeros(j-1)
    v=1
    for k in 1:length(t_aGFSS)
        if t_aGFSS[k]>T[1,k]
            t_change[v]=k
            v=v+1
        end
    end
    return t_change
end

#detection du changement dans les noeuds
function detect_n_change(sig, t_iaGFSS, ρ , d, v, t, T; λ = 0.01, Λ=0.1)
    (nv, nt) = size(sig)

    j=1
    t=Int(t)
    for k in 1:size(t_iaGFSS,1)
        if t_iaGFSS[k,t]>T[k,t]
            j=j+1
        end
    end
    n_change=zeros(j-1)
    v=1
    for k in 1:size(t_iaGFSS,1)
        if t_iaGFSS[k,t]>T[k,t]
            n_change[v]=k
            v=v+1
        end
    end
    return n_change
end

#filtrage adaptatif
function get_threshold(GFSS, n, a, init)
    if length(size(GFSS))==2
        (nv,nt) = size(GFSS)
        T=zeros(nv,nt)
        for i in 1:nv
            S1=GFSS[i,1]
            S2=GFSS[i,1]+maximum(GFSS)/init
            T[i,1]=a*S2
            for k in 2:nt
                S1=GFSS[i,k]
                S2=(1-n)*S2+n*S1
                T[i,k]=a*S2
            end
        end
    else
        nt=length(GFSS)
        T=zeros(1,nt)
        S1=GFSS[1]
        S2=GFSS[1]+maximum(GFSS)/init
        T[1,1]=a*S2
        for k in 2:nt
            S1=GFSS[k]
            S2=(1-n)*S2+n*S1
            T[1,k]=a*S2
        end
    end


    return T
end

#detection du changement dans les clusters
function detect_change(node_labels,sig,t_aGFSS, t_iaGFSS, t_change, ρ , d, v,T, T2, s; λ = 0.01, Λ=0.1)
    detect=zeros(maximum(node_labels))
    for k in 1:length(t_change)
        n_change=detect_n_change(sig ,t_iaGFSS ,ρ , d, v,t_change[k],T; λ = 0.01, Λ=0.1 )
        for j in 1:maximum(node_labels)
            for u in 1:length(n_change)
                for i in 1:length(findall(in(j), node_labels))
                    if Int(n_change[u]) == findall(in(j), node_labels)[i]
                        detect[j]+=1
                    end
                end
            end
        end
    end
    for k in 1:length(detect)
        if detect[k]>=s
            println("le groupe $k a changé de comportement")
        end
    end
    return detect
end

#tableau de seuils
function tab_threshold(min, max, pas, nt)
    nbvaleur=Int(round((max-min)/pas))
    T1000=zeros(nbvaleur,nt)
    for j in 1:nbvaleur
        T1000[j,:]=j*pas*ones(512)+min*ones(512)
    end
    return T1000
end
