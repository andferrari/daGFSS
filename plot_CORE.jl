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

retard=readdlm("donnees/retard.csv")
retard2=readdlm("donnees/retard2.csv")
retard3=readdlm("donnees/retard3.csv")
retard4=readdlm("donnees/retard4.csv")


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

Plots.reset_defaults()
Plots.scalefontsizes(1.5)

plot(pfausse2[1:745],retard2, xlabel="pfa", ylabel ="retard",  label="daGFSS", w=2, dpi=300)
plot!(pfausse3, retard3[57:1123], label="daGFSS Independent", w=2, dpi=300)
plot!(pfausse[20:729], retard[20:729], label="daGFSS with 2-norm", w=2, dpi=300)
plot!(pfausse4[25:211], retard4[25:211], label="daGFSS Centralized", w=2, dpi=300)
