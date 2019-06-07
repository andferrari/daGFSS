using Plots
using DelimitedFiles

pfausse = readdlm("donnees/pfausse.csv")
pfausse2 = readdlm("donnees/pfausse2.csv")
pfausse3 = readdlm("donnees/pfausse3.csv")
pfausse4 = readdlm("donnees/pfausse4.csv")

pdetect = readdlm("donnees/pdetect.csv")
pdetect2 = readdlm("donnees/pdetect2.csv")
pdetect3 = readdlm("donnees/pdetect3.csv")
pdetect4 = readdlm("donnees/pdetect4.csv")

pyplot()

Plots.reset_defaults()
Plots.scalefontsizes(1.5)
#plot(pfausse4, pdetect4, ratio=:equal, xlabel="Pfa", ylabel="Pd", label="Centralized detection \$\\|\\mathbf{d}_t\\|\$", w=3)
#plot!(pfausse2, pdetect2, label=L"\mathbf{t}_{daGFSS}(i)", w=3)
#plot!(pfausse3, pdetect3, label="Independent detection \$|\\mathbf{d}_t(i)|\$", w=3)
#plot!(pfausse, pdetect, label=L"\sum_{k\in N_\mathcal{G}(i)} \ \mathbf{d}_t(k)^2", w=3)

plot(pfausse2, pdetect2, label="daGFSS", w=3, dpi=300)
plot!(pfausse, pdetect, label="daGFSS with 2-norm", w=3, dpi=300)
plot!(pfausse4, pdetect4, ratio=:equal, xlabel="Pfa", ylabel="Pd", label="daGFSS Centralized", w=3, dpi=300)
plot!(pfausse3, pdetect3, label="daGFSS Independent", w=3, dpi=300)
