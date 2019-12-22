using Distributed
# don't forget addprocs( ? ) 

@everywhere using LightGraphs
using GraphIO
using EzXML

@everywhere using LinearAlgebra
@everywhere using Statistics
@everywhere include("gfss_func.jl")

using DelimitedFiles
# uncomment the  end to write the results in data/

g = loadgraph(path*"data/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
L = NormalizedLaplacian(g)
A=adjacency_matrix(g)
d, v = eigen(Array(L));
L = L - I
node_labels = Int.(label_propagation(g, 10000)[1])


nb = 10000 
nt = 512
ρ = 1.0/0.3
λ = 0.01
Λ = 0.1


poles = [ 0.9284586365913845 + 0.6691948262233165im 0.9284586365913845 - 0.6691948262233165im -0.9284586223955065 + 0.6691948202913867im -0.9284586223955065 - 0.6691948202913867im]
residues = [-0.09550841212039587 - 0.10204555134224505im -0.09550841212039587 + 0.10204555134224504im -0.023277450874456127 - 0.8479373939514138im  -0.023277450874456127 + 0.8479373939514138im ]
φ, ψ = calcul_psi_phi(poles, residues)
c = 0.6682305081233931

variancei = readdlm("data/variancei.csv")
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
fin=512


mini = round(minimum(t_square_daGFSS[:,init:fin]))
maxi = round(maximum(t_square_daGFSS[:,init:fin]))
step = 8
x1 = (mini):step:maxi
T1000 = tab_threshold(mini,maxi,step,nt) #creation d'un tableau de seuils
pdetect = compute_pd(nt, sig1, t_square_daGFSS, t_square_neigh, ρ ,d ,v ,T1000)
delay = compute_delay(nt, sig1, t_square_daGFSS, t_square_neigh, ρ ,d ,v ,T1000)
pfa, tdetecttot = compute_pfa(nt, t_square_daGFSS_sc, ρ ,d ,v ,T1000)


mini2 = round(minimum(t_daGFSS_square[:,init:fin]))
maxi2 = round(maximum(t_daGFSS_square[:,init:fin]))
step = (maxi2-mini2)/size(T1000)[1]
x2 = mini2:step:maxi2
T10002 = tab_threshold(mini2,maxi2,step,nt)
pdetect2 = compute_pd(nt, sig1, t_daGFSS_square, t_neigh_square, ρ ,d ,v ,T10002)
delay2 = compute_delay(nt, sig1, t_daGFSS_square, t_neigh_square, ρ ,d ,v ,T10002)
pfa2, tdetecttot2 = compute_pfa(nt, t_daGFSS_square_sc, ρ ,d ,v ,T10002)


mini3 = round(minimum(t_daGFSS[:,init:fin]))
maxi3 = round(maximum(t_daGFSS[:,init:fin]))
step = (maxi3-mini3)/size(T1000)[1]
x3 = mini3:step:maxi3
T10003 = tab_threshold(mini3,maxi3,step,nt)
pdetect3 = compute_pd(nt, sig1, t_daGFSS, t_square_diaGFSS, ρ ,d ,v ,T10003)
delay3 = compute_delay(nt, sig1, t_daGFSS, t_square_diaGFSS, ρ ,d ,v ,T10003)
pfa3, tdetecttot3 = compute_pfa(nt, t_daGFSS_sc, ρ ,d ,v ,T10003)


mini4 = round(minimum(sumtot_daGFSS[:,init:fin]))
maxi4 = round(maximum(sumtot_daGFSS[:,init:fin]))
step = (maxi4-mini4)/size(T1000)[1]
x4 = mini4:step:maxi4
T10004 = tab_threshold(mini4,maxi4,step,nt) #creation d'un tableau de seuils
pdetect4 = compute_pd(nt, sig1, sumtot_daGFSS, t_square_diaGFSS, ρ ,d ,v ,T10004)
delay4 = compute_delay(nt, sig1, sumtot_daGFSS, t_square_diaGFSS, ρ ,d ,v ,T10004)
pfa4, tdetecttot4= compute_pfa(nt, sumtot_daGFSS_sc, ρ ,d ,v ,T10004)


# writedlm(pfa, "data/pfausse_delay.csv")
# writedlm(pfa2, "data/pfausse2_delay.csv")
# writedlm(pfa3, "data/pfausse3_delay.csv")
# writedlm(pfa4, "data/pfausse4_delay.csv")

# writedlm(pdetect, "data/pdetect_delay.csv")
# writedlm(pdetect2, "data/pdetect2_delay.csv")
# writedlm(pdetect3, "data/pdetect3_delay.csv")
# writedlm(pdetect4, "data/pdetect4_delay.csv")

# writedlm(delay, "data/delay_var.csv")
# writedlm(delay2, "data/delay2_var.csv")
# writedlm(delay3, "data/delay3_var.csv")
# writedlm(delay4, "data/delay4_var.csv")

