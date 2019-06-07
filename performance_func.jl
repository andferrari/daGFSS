function performance_algo(nt, t_aGFSS,ρ ,d ,v ,T1000)
    n_rupt = 400
    Δ_rupt = 30
    inter=fin-init
    Texp=zeros(1,inter)
    nsignauxdetect=zeros(size(T1000)[1])
    tdetect=zeros(nb,size(T1000)[1])
    tdetecttot=zeros(size(T1000)[1])
    e2=zeros(size(T1000)[1])


    for k in 1:nb
        for i in 1:size(T1000)[1] #test de tous les seuils
            Texp[1,:]=T1000[i,init+1:fin]
            t_changeexp=detect_t_change(t_aGFSS[k,init+1:fin],ρ , d, v,Texp; λ = 0.01, Λ=0.1 ) #on detect si il y a changement pour le suil étudié
            if t_changeexp!=zeros(0) #changement
                    e2[i]=e2[i]+1 #fausse alarme
                    nsignauxdetect[i]+=1 #donc le signal a été detecté
                    tdetect[k,i]=t_changeexp[1] #premier temps de detection
            end
        end
    end
    pfausse=e2/nb

    for i in 1:size(T1000)[1] #moyenne pour chaque seuil de tous les temps
        tdetecttot[i]=sum(tdetect[:,i])/nsignauxdetect[i]
    end


    return pfausse, tdetecttot
end

function performance_algo2(nt,sig1, t_aGFSS, t_diaGFSS, ρ ,d ,v ,T1000, T2000)
    inter=fin-init
    Texp=zeros(1,inter)
    retard=zeros(nb,size(T1000)[1])
    retardtot=zeros(size(T1000)[1])
    nsignauxdetect=zeros(size(T1000)[1])
    tdetect=zeros(nb,size(T1000)[1])
    tdetecttot=zeros(size(T1000)[1])
    e=zeros(size(T1000)[1])


    for k in 1:nb
        println(k)
        for i in 1:size(T1000)[1] #test de tous les seuils
            a=0
            Texp[1,:]=T1000[i,init+1:fin]
            t_changeexp=detect_t_change(t_aGFSS[k,init+1:fin],ρ , d, v,Texp; λ = 0.01, Λ=0.1 ) #on detect si il y a changement pour le suil étudié
            if t_changeexp!=zeros(0) #changement
                nsignauxdetect[i]+=1 #donc le signal a été detecté
                e[i]=e[i]+1 #nb de signaux détectés
                tdetect[k,i]=t_changeexp[1] #premier temps de detection
                retard[k,i]=t_changeexp[1]+init-n_rupt #si bonne detection calcul du retard
            end
        end
    end
    pdetect=e/nb

    for i in 1:size(T1000)[1] #moyenne pour chaque seuil de tous les temps
        retardtot[i]=sum(retard[:,i])/nsignauxdetect[i]
        tdetecttot[i]=sum(tdetect[:,i])/nsignauxdetect[i]
    end


    return pdetect, retardtot, tdetecttot
end

function performance_algo3(nt,sig1, t_aGFSS, t_diaGFSS, ρ ,d ,v ,T1000, T2000)
    Texp=zeros(1,inter)

    retard=zeros(nb,size(T1000)[1])
    retardtot=zeros(size(T1000)[1])
    nsignauxdetect=zeros(size(T1000)[1])
    tdetect=zeros(nb,size(T1000)[1])
    tdetecttot=zeros(size(T1000)[1])
    e=zeros(size(T1000)[1])
    e2=zeros(size(T1000)[1])
    a=0


    for k in 1:nb
        println(k)
        for i in 1:size(T1000)[1] #test de tous les seuils
            a=0
            Texp[1,:]=T1000[i,init+1:fin]
            t_changeexp=detect_t_change(t_aGFSS[k,init+1:fin],ρ , d, v,Texp; λ = 0.01, Λ=0.1 ) #on detect si il y a changement pour le suil étudié
            if t_changeexp!=zeros(0) #changement
                nsignauxdetect[i]+=1 #donc le signal a été detecté
                e[i]=e[i]+1 #nb de signaux détectés
                tdetect[k,i]=t_changeexp[1] #premier temps de detection
                if t_changeexp[1]<n_rupt-init+1 || t_changeexp[length(t_changeexp)]>=n_rupt-init+1+Δ_rupt
                    e2[i]=e2[i]+1 #fausse alarme
                else
                    retard[k,i]=t_changeexp[1]+init-n_rupt #si bonne detection calcul du retard
                    for j in 1:length(t_changeexp)
                        n_changeexp=detect_n_change(sig1, t_diaGFSS[:,250:512] ,ρ , d, v,t_changeexp[j],T2000; λ = 0.01, Λ=0.1 )
                        if n_changeexp[1]<minimum(findall(in(8), node_labels)) || n_changeexp[length(n_changeexp)]>maximum(findall(in(8), node_labels))
                            a=1 #fausse alarme
                        end
                    end
                    if a==1
                        e2[i]=e2[i]+1 #fausse alarme
                    end
                end
            end
        end
    end
    pdetect=e/nb
    pfausse=e2/nb

    for i in 1:size(T1000)[1] #moyenne pour chaque seuil de tous les temps
        retardtot[i]=sum(retard[:,i])/nsignauxdetect[i]
        tdetecttot[i]=sum(tdetect[:,i])/nsignauxdetect[i]
    end


    return pdetect, retardtot, tdetecttot
end
