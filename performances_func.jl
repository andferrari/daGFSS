function calculate_pfa(nt, t_aGFSS, ρ, d, v, T1000)
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

function calculate_pd(nt, sig1, t_aGFSS, t_diaGFSS, ρ, d, v, T1000)
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

function calculate_delay(nt, sig1, t_aGFSS, t_diaGFSS, ρ, d, v, T1000)
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
