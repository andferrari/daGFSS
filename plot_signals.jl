using LightGraphs
using GraphIO
using EzXML
using ColorSchemes
using LinearAlgebra
using Plots

using LaTeXStrings
pyplot()
using Statistics


include("signals_func.jl")
include("gfss_func.jl")

g = loadgraph("data/MyGraph.graphml", GraphIO.GraphML.GraphMLFormat())
L = NormalizedLaplacian(g)
d, v = eigen(Array(L));
L = L - I

node_labels = Int.(label_propagation(g, 10000)[1])


nb = 1000
nt = 512
ρ = 1.0/0.3
σ2 = 7
init = 300
fin = 512
x = init:1:fin


sig1 = gener_sigg(g, node_labels, 2.4,1, 2.0, Δ_rupt = 112, σ2 = 7);
s1 = sig1[findall(in(1), node_labels), x]
m1 = mean(s1,dims=1)
plt1 = plot(x,m1',ribbon=sqrt(3)*ones(512), color=RGB(1,136/255,5/255),label="")

s2 = sig1[findall(in(2), node_labels), x]
m2 = mean(s2,dims=1)
plot!(x,m2',ribbon=sqrt(3)*ones(512), color=RGB(0,196/255,1),label="")

s3 = sig1[findall(in(3), node_labels), x]
m3 = mean(s3,dims=1)
plot!(x,m3',ribbon=sqrt(3)*ones(512), color=RGB(211/255,179/255,176/255),label="")

s4 = sig1[findall(in(4), node_labels), x]
m4 = mean(s4,dims=1)
plot!(x,m4',ribbon=sqrt(3)*ones(512), color=RGB(223/255,137/255,1),label="")

s5 = sig1[findall(in(5), node_labels), x]
m5 = mean(s5,dims=1)
plot!(x,m5',ribbon=sqrt(3)*ones(512), color=RGB(1,85/255,132/255),label="")

s6 = sig1[findall(in(6), node_labels), x]
m6 = mean(s6,dims=1)
plot!(x,m6',ribbon=sqrt(3)*ones(512), color=RGB(0,189/255,148/255),label="")

s7 = sig1[findall(in(7), node_labels),x]
m7 = mean(s7,dims=1)
plot!(x,m7',ribbon=sqrt(3)*ones(512), color=RGB(115/255,192/255,0),label="")

s8 = sig1[findall(in(8), node_labels), x]
m8 = mean(s8,dims=1)
plot!(x,m8',ribbon=sqrt(3)*ones(512), color=RGB(76/255,70/255,62/255), label="")

plt9 = plot!(xlab=L"time",ylab=L"signals \ on \ vertices")
vline!([400], w=3, label="")



t_iaGFSS = iaGFSS(sig1, ρ, d, v; λ = 0.01, Λ = 0.1)
t_aGFSS = [norm(t_iaGFSS[:,k])^2 for k in 1:512]
p3 = plot(x,t_aGFSS[x],label="", xlab = L"time", ylab=L"t_{aGFSS}")

plot(plt9,p3, layout=Plots.grid(2,1,heights=[0.7,0.3]), dpi=300)
