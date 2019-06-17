using Plots
using DelimitedFiles

pfausse = readdlm("donnees/pfausseret.csv")
pfausse2 = readdlm("donnees/pfausse2ret.csv")
pfausse3 = readdlm("donnees/pfausse3ret.csv")
pfausse4 = readdlm("donnees/pfausse4ret.csv")

pdetect = readdlm("donnees/pdetectret.csv")
pdetect2 = readdlm("donnees/pdetect2ret.csv")
pdetect3 = readdlm("donnees/pdetect3ret.csv")
pdetect4 = readdlm("donnees/pdetect4ret.csv")

retard=readdlm("donnees/retardvar.csv")
retard2=readdlm("donnees/retard2var.csv")
retard3=readdlm("donnees/retard3var.csv")
retard4=readdlm("donnees/retard4var.csv")


pyplot()

Plots.reset_defaults()
Plots.scalefontsizes(1.5)
#plot(pfausse4, pdetect4, ratio=:equal, xlabel="Pfa", ylabel="Pd", label="Centralized detection \$\\|\\mathbf{d}_t\\|\$", w=3)
#plot!(pfausse2, pdetect2, label=L"\mathbf{t}_{daGFSS}(i)", w=3)
#plot!(pfausse3, pdetect3, label="Independent detection \$|\\mathbf{d}_t(i)|\$", w=3)
#plot!(pfausse, pdetect, label=L"\sum_{k\in N_\mathcal{G}(i)} \ \mathbf{d}_t(k)^2", w=3)

plot(pfausse2, pdetect2, label="daGFSS", w=2, dpi=300)
plot!(pfausse, pdetect, label="daGFSS with 2-norm", w=2, dpi=300)
plot!(pfausse4, pdetect4, ratio=:equal, xlabel="Pfa", ylabel="Pd", label="daGFSS Centralized", w=2, dpi=300)
plot!(pfausse3, pdetect3, label="daGFSS Independent", w=2, dpi=300)
savefig("newCore.pdf")

Plots.reset_defaults()
Plots.scalefontsizes(1.5)

plot(pfausse2,retard2, xlabel="Pfa", ylabel ="Detect Delay",  label="daGFSS", w=2, dpi=300)
plot!(pfausse3, retard3, label="daGFSS Independent", w=2, dpi=300)
plot!(pfausse, retard, label="daGFSS with 2-norm", w=2, dpi=300)
plot!(pfausse4, retard4, label="daGFSS Centralized", w=2, dpi=300)
