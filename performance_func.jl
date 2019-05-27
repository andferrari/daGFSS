function performance_algo(nt, t_aGFSS,ρ ,d ,v ,T1000)
    init=250
    fin=512
    nb=1000
    inter=fin-init # taille de la fenêtre étudiée
    n_rupt = 400
    Δ_rupt = 30


    nbpics=zeros(nb,size(T1000)[1])
    moy_interpics=zeros(nb,size(T1000)[1])
    moytot_interpics=zeros(size(T1000)[1])
    retard=zeros(nb,size(T1000)[1])
    retardtot=zeros(size(T1000)[1])
    nsignauxdetect=zeros(size(T1000)[1])
    nsignaux2pics=zeros(size(T1000)[1])
    Texp=zeros(1,inter)
    tdetect=zeros(nb,size(T1000)[1])
    tdetecttot=zeros(size(T1000)[1])
    e=zeros(size(T1000)[1])
    e2=zeros(size(T1000)[1])

    t_interpic=zeros(nb,nt,size(T1000)[1])

    for k in 1:nb
        for i in 1:size(T1000)[1] #test de tous les seuils
            Texp[1,:]=T1000[i,init+1:fin]
            t_changeexp=detect_t_change(t_aGFSS[k,init+1:fin],ρ , d, v,Texp; λ = 0.01, Λ=0.1 ) #on detect si il y a changement pour le suil étudié
            if t_changeexp!=zeros(0) #changement
                nbpics[k,i]+=1  #donc detection de pics
                nsignauxdetect[i]+=1 #donc le signal a été detecté
                e[i]=e[i]+1 #nb de signaux détectés
                tdetect[k,i]=t_changeexp[1] #premier temps de detection
                if t_changeexp[1]<n_rupt-init+1 || t_changeexp[length(t_changeexp)]>=n_rupt-init+1+Δ_rupt
                    e2[i]=e2[i]+1 #fausse alarme
                else
                    retard[k,i]=t_changeexp[1]+init-n_rupt #si bonne detection calcul du retard
                end
                if length(t_changeexp)>1
                    for j in 2:length(t_changeexp)
                        if t_changeexp[j]-t_changeexp[j-1]>1
                            nbpics[k,i]+=1
                            t_interpic[k,Int(t_changeexp[j-1])+init-1,i]+=t_changeexp[j]-t_changeexp[j-1]
                        end
                    end
                    if  nbpics[k,i]>1
                        nsignaux2pics[i]+=1
                    end
                end
            end
        end
    end
    pdetect=e/nb
    pfausse=e2/nb

    for k in 1:nb
        for i in 1:size(T1000)[1]
            if nbpics[k,i]>1
                moy_interpics[k,i]=sum(t_interpic[k,:,i])/nbpics[k,i]
            end
        end
    end

    for i in 1:size(T1000)[1] #moyenne pour chaque seuil de tous les temps
        moytot_interpics[i]=sum(moy_interpics[:,i])/nsignaux2pics[i]
        retardtot[i]=sum(retard[:,i])/nsignauxdetect[i]
        tdetecttot[i]=sum(tdetect[:,i])/nsignauxdetect[i]
    end


    return pdetect, pfausse, moytot_interpics, retardtot, tdetecttot
end
