using Plots
using DelimitedFiles

pfausse = readdlm("pfausse.csv")
pfausse2 = readdlm("pfausse2.csv")
pfausse3 = readdlm("pfausse3.csv")
pfausse4 = readdlm("pfausse4.csv")

pdetect = readdlm("pdetect.csv")
pdetect2 = readdlm("pdetect2.csv")
pdetect3 = readdlm("pdetect3.csv")
pdetect4 = readdlm("pdetect4.csv")

pyplot()
plot(pfausse2, pdetect2, ratio=:equal, xlabel=L"pfa", ylabel=L"pd", label=L"Centralized \ detection \ \mathbf{d}_t ", w=3)
plot!(pfausse3, pdetect3, label=L"\mathbf{t}_{daGFSS}(i)^2", w=3)
plot!(pfausse, pdetect, label=L"independent \ detection \ \mathbf{d}_t(i)^2", w=3)
plot!(pfausse4, pdetect4, label=L"\sum_{k\in N_\mathcal{G}(i)} \ \mathbf{d}_t(k)^2", dpi=600, w=3)
