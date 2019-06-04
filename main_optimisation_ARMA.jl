using Optim
using LightGraphs
using GraphIO
using GraphPlot
using EzXML
using LinearAlgebra
using Plots
pyplot()
using LaTeXStrings

include("optimisation_func.jl")

# calcul des valeurs propres
g = loadgraph("donnees/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
L = LightGraphs.laplacian_matrix(g)
d, v = eigen(Array(L));
λmax=maximum(d)
λmin=minimum(d)
a=0.4
b=λmax -1

#choix des paramètres
ρ = 1.0
n_ψ = 256
ψ= range(-1/λmax, stop = 1/λmax, length = n_ψ)
n_φ = 256
φ = range(-0.8, stop = 0.8, length = n_φ)
resul=zeros(n_φ,n_ψ)

#calcul des paramètres dans une fenêtre delimitée
for i in 1:n_φ
    for j in 1:n_ψ
        c=optim_c(a, b, ρ, φ[i], ψ[j])
            resul[i,j]=Integrale1(ρ, c, ρ, φ[i], ψ[j])-Integrale1(a,c, ρ, φ[i], ψ[j])+Integrale2(b, c, ρ , φ[i], ψ[j])-Integrale2(ρ, c, ρ , φ[i], ψ[j])
    end
end
resul

heatmap(ψ,φ, (resul).^0.001, xlabel = "ψ", ylabel = "φ")
plot(ψ,φ, (resul).^0.001, st=:surface, xlabel = "ψ", ylabel = "φ")

#valeurs optimales
minimum(resul)
argmin(resul)
φ[103]
ψ[200]
c=optim_c(a, b, ρ, φ[103], ψ[200])

#comparaison des deux filtres
x=0:0.01:30
opti=φ[103]./(ones(length(x))-ψ[200]*x)+ones(length(x))*c #filtre obtenu
plot(x,opti,label=L"h(\mu)")

filtre=zeros(length(x))
for k in 1:length(x)
    if x[k]<ρ
        filtre[k]=1
    else
        filtre[k]=sqrt(ρ/x[k])
    end
end
filtre #filtre optimal
plot!(x,filtre, xlab = L"μ", label=L"h^\ast(\mu)")
