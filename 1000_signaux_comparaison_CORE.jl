@everywhere using LightGraphs
using GraphIO
using EzXML
@everywhere using LinearAlgebra
using Plots
@everywhere using LaTeXStrings
@everywhere using DelimitedFiles
@everywhere using Statistics


@everywhere path = "/home/verduci/daGFSS/"
# load graph

@everywhere include(path*"gfss_func.jl")
include(path*"performance_func.jl")
include(path*"detection_func.jl")
@everywhere include(path*"signaux_func.jl")

g = loadgraph(path*"donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
#L = LightGraphs.laplacian_matrix(g)
L = NormalizedLaplacian(g)
A=adjacency_matrix(g)
d, v = eigen(Array(L));
λmax = maximum(d)

L= L - (λmax/2)I
node_labels = Int.(label_propagation(g, 10000)[1])


nb=10000 #nb de signaux créés
nt=512 #temps d'étude
ρ = 1.0/0.3
λ=0.01
Λ=0.1

#lancer synth_gfss_parallele pour avoir les poles et les residues
poles=[0.7344635585485741 + 0.5293681924014867im 0.7344635585485741 - 0.5293681924014867im -0.7344635202651353 + 0.5293681751822871im -0.7344635202651353 - 0.5293681751822871im]
residues=[-0.05322205451164147 - 0.08206089787078102im -0.05322205451164147 + 0.08206089787078102im 0.16044187053514897 - 0.6853621323079733im 0.16044187053514897 + 0.6853621323079733im]
φ, ψ = calcul_psi_phi(poles, residues)
c=0.6682305081233931
variancei=readdlm("variancei.csv")
σi=sqrt.(variancei)





sig1_sc=SharedArray{Float64}(250, nt*nb)
t_diaGFSS_sc=SharedArray{Float64}(250, nt*nb)
t_daGFSS_sc=SharedArray{Float64}(nb,nt)
t_carre_diaGFSS_sc=SharedArray{Float64}(250, nt*nb)
t_carre_daGFSS_sc=SharedArray{Float64}(nb,nt)
t_daGFSS_carre_sc=SharedArray{Float64}(nb,nt)
sumtot_daGFSS_sc=SharedArray{Float64}(nb,nt)



#création des signaux
@distributed for k in 1:nb
    sig1_sc[:,(k-1)*512+1:k*512]= gener_sigg(g, node_labels, 1.0, 2 , 0.5, σ2 = 7);
    t_carre_diaGFSS_sc[:,(k-1)*512+1:k*512], t_diaGFSS_sc[:,(k-1)*512+1:k*512] = diaGFSS(sig1_sc[:,(k-1)*512+1:k*512], L, ψ, φ, c; λ = 0.01, Λ=0.1)
end

t_voisin_sc=-(t_diaGFSS_sc + A*t_diaGFSS_sc)
t_voisin_carre_sc = (t_voisin_sc.^2)./σi

t_carre_voisin_sc=(t_carre_diaGFSS_sc + A*t_carre_diaGFSS_sc)./σi
t_carre_diaGFSS_sc=t_carre_diaGFSS_sc./σi

@distributed for k in 1:nb
    t_carre_daGFSS_sc[k,:] = [norm(t_carre_voisin_sc[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS_carre_sc[k,:] = [norm(t_voisin_carre_sc[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS_sc[k,:] = [norm(t_carre_diaGFSS_sc[:,i])^2 for i in (k-1)*512+1:k*512]
end
t_carre_daGFSS_sc

@distributed for k in 1:nb
    sumtot_daGFSS_sc[k,:]=sum(t_carre_diaGFSS_sc[:,(k-1)*512+1:k*512],dims=1)
end
sumtot_daGFSS_sc





sig1=SharedArray{Float64}(250, nt*nb)
t_diaGFSS=SharedArray{Float64}(250, nt*nb)
t_daGFSS=SharedArray{Float64}(nb,nt)
t_carre_diaGFSS=SharedArray{Float64}(250, nt*nb)
t_carre_daGFSS=SharedArray{Float64}(nb,nt)
t_daGFSS_carre=SharedArray{Float64}(nb,nt)
sumtot_daGFSS=SharedArray{Float64}(nb,nt)


#création des signaux
@distributed for k in 1:nb
    sig1[:,(k-1)*512+1:k*512]= gener_sigg(g, node_labels, 1.0, 1 , 0.5, σ2 = 7);
    t_carre_diaGFSS[:,(k-1)*512+1:k*512], t_diaGFSS[:,(k-1)*512+1:k*512] = diaGFSS(sig1[:,(k-1)*512+1:k*512], L, ψ, φ, c; λ = 0.01, Λ=0.1)
end

t_voisin=-(t_diaGFSS + A*t_diaGFSS)
t_voisin_carre = (t_voisin.^2)./σi

t_carre_voisin=(t_carre_diaGFSS + A*t_carre_diaGFSS)./σi
t_carre_diaGFSS=t_carre_diaGFSS./σi

@distributed for k in 1:nb
    t_carre_daGFSS[k,:] = [norm(t_carre_voisin[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS_carre[k,:] = [norm(t_voisin_carre[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS[k,:] = [norm(t_carre_diaGFSS[:,i])^2 for i in (k-1)*512+1:k*512]
end
t_carre_daGFSS

@distributed for k in 1:nb
    sumtot_daGFSS[k,:]=sum(t_carre_diaGFSS[:,(k-1)*512+1:k*512],dims=1)
end
sumtot_daGFSS






#t_carre
init=300
init_sc=300
fin=512
n_rupt = 400
Δ_rupt = 30


plot(t_carre_daGFSS[100,init:fin],xlabel="temps", ylabel="t_carre_daGFSS")
mini=round(minimum(t_carre_daGFSS[:,init:fin]))
maxi=round(maximum(t_carre_daGFSS[:,init:fin]))
pas=10
x1=(mini):pas:maxi
T1000=tab_threshold(mini,maxi,pas,nt) #creation d'un tableau de seuils

plot(t_carre_voisin[:,init:fin]',xlabel="temps", ylabel="t_carre_aGFSS", label="")
T20001=80*ones(250,512)

pdetect = performance_algo2(nt, sig1, t_carre_daGFSS, t_carre_voisin, ρ ,d ,v ,T1000, T20001)
retard = calcul_retard(nt, sig1, t_carre_daGFSS, t_carre_voisin, ρ ,d ,v ,T1000, T20001)

plot(t_carre_daGFSS_sc[100,init_sc:fin],xlabel="temps", ylabel="t_carre_aGFSS")

pfausse, tdetecttot= performance_algo(nt, t_carre_daGFSS_sc, ρ ,d ,v ,T1000)

plot(pdetect, retardtot)

#t_voisin_carre
plot(t_daGFSS_carre[100,init:fin],xlabel="temps", ylabel="t_daGFSS_carre")
mini2=round(minimum(t_daGFSS_carre[:,init:fin]))
maxi2=round(maximum(t_daGFSS_carre[:,init:fin]))
pas=(maxi2-mini2)/size(T1000)[1]
x2=(mini2):pas:maxi2
T10002=tab_threshold(mini2,maxi2,pas,nt)

plot(t_voisin_carre[:,init:fin]',xlabel="temps", ylabel="t_carre_aGFSS", label="")
T2000=1000*ones(250,512)

pdetect2 = performance_algo2(nt, sig1, t_daGFSS_carre, t_voisin_carre, ρ ,d ,v ,T10002, T2000)
retard2 = calcul_retard(nt, sig1, t_daGFSS_carre, t_voisin_carre, ρ ,d ,v ,T10002, T2000)

plot(t_daGFSS_carre_sc[100,init_sc:fin],xlabel="temps", ylabel="t_daGFSS_carre")

pfausse2, tdetecttot2= performance_algo(nt, t_daGFSS_carre_sc, ρ ,d ,v ,T10002,)

plot(pdetect2, retardtot2)


#t
plot(t_daGFSS[100,init:fin],xlabel="temps", ylabel="t_daGFSS")
mini3=round(minimum(t_daGFSS[:,init:fin]))
maxi3=round(maximum(t_daGFSS[:,init:fin]))
pas=(maxi3-mini3)/size(T1000)[1]
x3=(mini3):pas:maxi3
T10003=tab_threshold(mini3,maxi3,pas,nt)

plot(t_carre_diaGFSS[:,init:fin]',xlabel="temps", ylabel="t_diaGFSS", label="")
T20002=10*ones(250,512)

pdetect3 = performance_algo2(nt, sig1, t_daGFSS, t_carre_diaGFSS, ρ ,d ,v ,T10003, T20002)
retard3 = calcul_retard(nt, sig1, t_daGFSS, t_carre_diaGFSS, ρ ,d ,v ,T10003, T20002)


plot(t_daGFSS_sc[100,init_sc:fin],xlabel="temps", ylabel="t_daGFSS")
pfausse3, tdetecttot3= performance_algo(nt, t_daGFSS_sc, ρ ,d ,v ,T10003)

plot(pdetect3, retardtot3)

#t_sum
plot(sumtot_daGFSS[100,init:fin],xlabel="temps", ylabel="t_sum_daGFSS")
mini4=round(minimum(sumtot_daGFSS[:,init:fin]))
maxi4=round(maximum(sumtot_daGFSS[:,init:fin]))
pas=(maxi4-mini4)/size(T1000)[1]
x4=(mini4):pas:maxi4
T10004=tab_threshold(mini4,maxi4,pas,nt) #creation d'un tableau de seuils


pdetect4 = performance_algo2(nt, sig1, sumtot_daGFSS, t_carre_diaGFSS, ρ ,d ,v ,T10004, T20002)
retard4 = calcul_retard(nt, sig1, sumtot_daGFSS, t_carre_diaGFSS, ρ ,d ,v ,T10004, T20004)


plot(sumtot_daGFSS_sc[100,init_sc:fin],xlabel="temps", ylabel="t_sum_daGFSS")


pfausse4, tdetecttot4= performance_algo(nt, sumtot_daGFSS_sc, ρ ,d ,v ,T10004)

plot(pdetect4, retardtot4)

pyplot()
#plot(pfausse2, pdetect2, ratio=:equal, xlabel=L"Pfa", ylabel=L"pd", label="daGFSS", w=2)
#plot!(pfausse3, pdetect3, label="daGFSS Independent", w=2)
#plot!(pfausse, pdetect, label="daGFSS with 2-norm", w=2)
#plot!(pfausse4, pdetect4, label="daGFSS Centralized", dpi=600, w=2)


plot(pdetect2,retardtot2, xlabel="Pd", ylabel ="Mean detecion delay", label="daGFSS")
plot!(pdetect3, retardtot3, label="daGFSS Independent")
plot!(pdetect, retardtot, label="daGFSS with 2-norm")
plot!(pdetect4, retardtot4, label="daGFSS Centralized",xlims=(0.04,1))


writedlm("pdetectret.csv",pdetect)
writedlm("pdetect2ret.csv",pdetect2)
writedlm("pdetect3ret.csv",pdetect3)
writedlm("pdetect4ret.csv",pdetect4)

writedlm("retardvar.csv",retardtot)
writedlm("retard2var.csv",retardtot2)
writedlm("retard3var.csv",retardtot3)
writedlm("retard4var.csv",retardtot4)
