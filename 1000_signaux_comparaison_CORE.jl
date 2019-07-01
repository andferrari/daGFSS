@everywhere using LightGraphs
using GraphIO
using EzXML
@everywhere using LinearAlgebra
using Plots
@everywhere using LaTeXStrings
@everywhere using DelimitedFiles
@everywhere using Statistics


@everywhere path = "/home/verduci/daGFSS/"
@everywhere include(path*"gfss_func.jl")
include(path*"performance_func.jl")
include(path*"detection_func.jl")
@everywhere include(path*"signaux_func.jl")

g = loadgraph(path*"donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
L = NormalizedLaplacian(g)
A=adjacency_matrix(g)
d, v = eigen(Array(L));
λmax = maximum(d)
L = L - (λmax/2)I
node_labels = Int.(label_propagation(g, 10000)[1])


nb = 10000 #
nt = 512
ρ = 1.0/0.3
λ = 0.01
Λ = 0.1


poles = [0.7344635585485741 + 0.5293681924014867im 0.7344635585485741 - 0.5293681924014867im -0.7344635202651353 + 0.5293681751822871im -0.7344635202651353 - 0.5293681751822871im]
residues = [-0.05322205451164147 - 0.08206089787078102im -0.05322205451164147 + 0.08206089787078102im 0.16044187053514897 - 0.6853621323079733im 0.16044187053514897 + 0.6853621323079733im]
φ, ψ = calcul_psi_phi(poles, residues)
c = 0.6682305081233931

variancei = readdlm("donnees/variancei.csv")
σi = sqrt.(variancei)



sig1_sc = SharedArray{Float64}(250, nt*nb)
t_diaGFSS_sc = SharedArray{Float64}(250, nt*nb)
t_daGFSS_sc = SharedArray{Float64}(nb,nt)
t_square_diaGFSS_sc = SharedArray{Float64}(250, nt*nb)
t_square_daGFSS_sc = SharedArray{Float64}(nb,nt)
t_daGFSS_square_sc = SharedArray{Float64}(nb,nt)
sumtot_daGFSS_sc = SharedArray{Float64}(nb,nt)


@distributed for k in 1:nb
    sig1_sc[:,(k-1)*512+1:k*512] = gener_sigg(g, node_labels, 1.0, 2 , 0.5, σ2 = 7);
    t_square_diaGFSS_sc[:,(k-1)*512+1:k*512], t_diaGFSS_sc[:,(k-1)*512+1:k*512] = diaGFSS(sig1_sc[:,(k-1)*512+1:k*512], L, ψ, φ, c; λ = 0.01, Λ=0.1)
end

t_neigh_sc = -(t_diaGFSS_sc + A*t_diaGFSS_sc)
t_neigh_square_sc = (t_neigh_sc.^2)./σi

t_square_neigh_sc = (t_square_diaGFSS_sc + A*t_square_diaGFSS_sc)./σi
t_square_diaGFSS_sc = t_square_diaGFSS_sc./σi


@distributed for k in 1:nb
    t_square_daGFSS_sc[k,:] = [norm(t_square_neigh_sc[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS_square_sc[k,:] = [norm(t_neigh_square_sc[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS_sc[k,:] = [norm(t_square_diaGFSS_sc[:,i])^2 for i in (k-1)*512+1:k*512]
end

@distributed for k in 1:nb
    sumtot_daGFSS_sc[k,:] = sum(t_square_diaGFSS_sc[:,(k-1)*512+1:k*512],dims=1)
end



sig1 = SharedArray{Float64}(250, nt*nb)
t_diaGFSS = SharedArray{Float64}(250, nt*nb)
t_daGFSS = SharedArray{Float64}(nb,nt)
t_square_diaGFSS = SharedArray{Float64}(250, nt*nb)
t_square_daGFSS = SharedArray{Float64}(nb,nt)
t_daGFSS_square = SharedArray{Float64}(nb,nt)
sumtot_daGFSS = SharedArray{Float64}(nb,nt)



@distributed for k in 1:nb
    sig1[:,(k-1)*512+1:k*512] = gener_sigg(g, node_labels, 1.0, 1 , 0.5, σ2 = 7);
    t_square_diaGFSS[:,(k-1)*512+1:k*512], t_diaGFSS[:,(k-1)*512+1:k*512] = diaGFSS(sig1[:,(k-1)*512+1:k*512], L, ψ, φ, c; λ = 0.01, Λ=0.1)
end

t_neigh = -(t_diaGFSS + A*t_diaGFSS)
t_neigh_square = (t_neigh.^2)./σi

t_square_neigh = (t_square_diaGFSS + A*t_square_diaGFSS)./σi
t_square_diaGFSS = t_square_diaGFSS./σi


@distributed for k in 1:nb
    t_square_daGFSS[k,:] = [norm(t_square_neigh[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS_square[k,:] = [norm(t_neigh_square[:,i])^2 for i in (k-1)*512+1:k*512]
    t_daGFSS[k,:] = [norm(t_square_diaGFSS[:,i])^2 for i in (k-1)*512+1:k*512]
end

@distributed for k in 1:nb
    sumtot_daGFSS[k,:] = sum(t_square_diaGFSS[:,(k-1)*512+1:k*512],dims=1)
end




init=300
init_sc=300
fin=512
n_rupt = 400
Δ_rupt = 112


plot(t_square_daGFSS[400,init:fin],xlabel="temps", ylabel="t_square_daGFSS")
mini = round(minimum(t_square_daGFSS[:,init:fin]))
maxi = round(maximum(t_square_daGFSS[:,init:fin]))
pas = 30
x1 = mini:pas:maxi
T1000 = tab_threshold(mini,maxi,pas,nt) #creation d'un tableau de seuils

pdetect = calculate_pd(nt, sig1, t_square_daGFSS, t_square_neigh, ρ ,d ,v ,T1000)
delay = calculate_delay(nt, sig1, t_square_daGFSS, t_square_neigh, ρ ,d ,v ,T1000)

plot(t_square_daGFSS_sc[100,init_sc:fin],xlabel="temps", ylabel="t_square_aGFSS")
pfa, tdetecttot = calculate_pfa(nt, t_square_daGFSS_sc, ρ ,d ,v ,T1000)



plot(t_daGFSS_square[100,init:fin],xlabel="temps", ylabel="t_daGFSS_square")
mini2 = round(minimum(t_daGFSS_square[:,init:fin]))
maxi2 = round(maximum(t_daGFSS_square[:,init:fin]))
pas = (maxi2-mini2)/size(T1000)[1]
x2 = mini2:pas:maxi2
T10002 = tab_threshold(mini2,maxi2,pas,nt)

pdetect2 = calculate_pd(nt, sig1, t_daGFSS_square, t_neigh_square, ρ ,d ,v ,T10002)
delay2 = calculate_delay(nt, sig1, t_daGFSS_square, t_neigh_square, ρ ,d ,v ,T10002)

plot(t_daGFSS_square_sc[100,init_sc:fin],xlabel="temps", ylabel="t_daGFSS_square")
pfa2, tdetecttot2 = calculate_pfa(nt, t_daGFSS_square_sc, ρ ,d ,v ,T10002)




plot(t_daGFSS[100,init:fin],xlabel="temps", ylabel="t_daGFSS")
mini3 = round(minimum(t_daGFSS[:,init:fin]))
maxi3 = round(maximum(t_daGFSS[:,init:fin]))
pas = (maxi3-mini3)/size(T1000)[1]
x3 = mini3:pas:maxi3
T10003 = tab_threshold(mini3,maxi3,pas,nt)


pdetect3 = calculate_pd(nt, sig1, t_daGFSS, t_square_diaGFSS, ρ ,d ,v ,T10003)
delay3 = calculate_delay(nt, sig1, t_daGFSS, t_square_diaGFSS, ρ ,d ,v ,T10003)


plot(t_daGFSS_sc[100,init_sc:fin],xlabel="temps", ylabel="t_daGFSS")
pfa3, tdetecttot3 = calculate_pfa(nt, t_daGFSS_sc, ρ ,d ,v ,T10003)




plot(sumtot_daGFSS[100,init:fin],xlabel="temps", ylabel="t_sum_daGFSS")
mini4 = round(minimum(sumtot_daGFSS[:,init:fin]))
maxi4 = round(maximum(sumtot_daGFSS[:,init:fin]))
pas = (maxi4-mini4)/size(T1000)[1]
x4 = mini4:pas:maxi4
T10004 = tab_threshold(mini4,maxi4,pas,nt) #creation d'un tableau de seuils


pdetect4 = calculate_pd(nt, sig1, sumtot_daGFSS, t_square_diaGFSS, ρ ,d ,v ,T10004)
delay4 = calculate_delay(nt, sig1, sumtot_daGFSS, t_square_diaGFSS, ρ ,d ,v ,T10004)


plot(sumtot_daGFSS_sc[100,init_sc:fin],xlabel="temps", ylabel="t_sum_daGFSS")
pfa4, tdetecttot4= calculate_pfa(nt, sumtot_daGFSS_sc, ρ ,d ,v ,T10004)




pyplot()

Plots.reset_defaults()
Plots.scalefontsizes(1.5)
plot(pfa2,pdetect2, xlabel="Pfa", ylabel ="Mean detecion delay", label="daGFSS",w=2)
plot!(pfa3, pdetect3, label="daGFSS Independent", w=2)
plot!(pfa, pdetect, label="daGFSS with 2-norm",w=2)
plot!(pfa4, pdetect4, label="daGFSS Centralized", aspect_ratio=:equal, w=2)


writedlm("donnees/pfausseret.csv", pfa)
writedlm("donnees/pfausse2ret.csv", pfa2)
writedlm("donnees/pfausse3ret.csv", pfa3)
writedlm("donnees/pfausse4ret.csv", pfa4)

writedlm("donnees/pdetectret.csv", pdetect)
writedlm("donnees/pdetect2ret.csv", pdetect2)
writedlm("donnees/pdetect3ret.csv", pdetect3)
writedlm("donnees/pdetect4ret.csv", pdetect4)
