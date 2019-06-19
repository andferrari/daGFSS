function calculate_pfa(nt, t_aGFSS,ρ ,d ,v ,T1000)
    n_rupt = 400
    Δ_rupt = 30
    inter=fin-init
    Texp=zeros(1,inter)
    nsigdetect=zeros(size(T1000)[1])
    tdetect=zeros(nb,size(T1000)[1])
    tdetectfin=zeros(size(T1000)[1])
    e2=zeros(size(T1000)[1])


    for k in 1:nb
        for i in 1:size(T1000)[1] #test de tous les seuils
            Texp[1,:]=T1000[i,init+1:fin]
            t_changeexp=detect_t_change(t_aGFSS[k,init+1:fin],ρ , d, v,Texp; λ = 0.01, Λ=0.1 ) #on detect si il y a changement pour le suil étudié
            if t_changeexp!=zeros(0) #changement
                    e2[i]=e2[i]+1 #fausse alarme
                    nsigdetect[i]+=1 #donc le signal a été detecté
                    tdetect[k,i]=t_changeexp[1] #premier temps de detection
            end
        end
    end
    pfa=e2/nb

    for i in 1:size(T1000)[1] #moyenne pour chaque seuil de tous les temps
        tdetectfin[i]=sum(tdetect[:,i])/nsigdetect[i]
    end


    return pfa, tdetectfin
end

function calculate_pd(nt,sig1, t_aGFSS, t_diaGFSS, ρ ,d ,v ,T1000)
    inter=fin-init
    Texp=zeros(1,inter)
    nsigdetect=zeros(size(T1000)[1])
    tdetect=zeros(nb,size(T1000)[1])
    tdetectfin=zeros(size(T1000)[1])
    e=zeros(size(T1000)[1])


    for k in 1:nb
        println(k)
        for i in 1:size(T1000)[1] #test de tous les seuils
            a=0
            Texp[1,:]=T1000[i,init+1:fin]
            t_changeexp=detect_t_change(t_aGFSS[k,init+1:fin],ρ , d, v,Texp; λ = 0.01, Λ=0.1 ) #on detect si il y a changement pour le suil étudié
            if t_changeexp!=zeros(0) #changement
                e[i]=e[i]+1 #nb de signaux détectés
                tdetect[k,i]=t_changeexp[1] #premier temps de detection
                nsigdetect[i]+=1
            end
        end
    end
    pdetect=e/nb

    for i in 1:size(T1000)[1] #moyenne pour chaque seuil de tous les temps
        tdetectfin[i]=sum(tdetect[:,i])/nsigdetect[i]
    end


    return pdetect
end

function calculate_delay(nt,sig1, t_aGFSS, t_diaGFSS, ρ ,d ,v ,T1000)
    inter=fin-init
    Texp=zeros(1,inter)
    delay=zeros(nb,size(T1000)[1])
    delaytot=zeros(size(T1000)[1])
    nsigdetect=zeros(size(T1000)[1])


    for k in 1:nb
        println(k)
        for i in 1:size(T1000)[1] #test de tous les seuils
            a=0
            Texp[1,:]=T1000[i,init+1:fin]
            t_changeexp=detect_t_change(t_aGFSS[k,n_rupt:fin],ρ , d, v,Texp; λ = 0.01, Λ=0.1 )
            if t_changeexp!=zeros(0) #changement
                delay[k,i]=t_changeexp[1] #si bonne detection calcul du delay
                nsigdetect[i]+=1
            end
        end
    end

    for i in 1:size(T1000)[1] #moyenne pour chaque seuil de tous les temps
        delaytot[i]=sum(delay[:,i])/nsigdetect[i]
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
    R = (λ/(2-λ) + Λ/(2-Λ) + 2*λ*Λ/(λ+Λ-λ*Λ))*Q
    return R
end
