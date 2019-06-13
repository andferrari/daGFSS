using Plots
using DelimitedFiles

pfausse = readdlm("donnees/pfaussevar.csv")
pfausse2 = readdlm("donnees/pfausse2var.csv")
pfausse3 = readdlm("donnees/pfausse3var.csv")
pfausse4 = readdlm("donnees/pfausse4var.csv")

pdetect = readdlm("donnees/pdetectvar.csv")
pdetect2 = readdlm("donnees/pdetect2var.csv")
pdetect3 = readdlm("donnees/pdetect3var.csv")
pdetect4 = readdlm("donnees/pdetect4var.csv")

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

plot(pfausse2, pdetect2, label="daGFSS", w=1, dpi=300)
plot!(pfausse, pdetect, label="daGFSS with 2-norm", w=1, dpi=300)
plot!(pfausse4, pdetect4, ratio=:equal, xlabel="Pfa", ylabel="Pd", label="daGFSS Centralized", w=1, dpi=300)
plot!(pfausse3, pdetect3, label="daGFSS Independent", w=1, dpi=300)
xlims!((0,0.5))
ylims!((0.5,1))

Plots.reset_defaults()
Plots.scalefontsizes(1.5)

plot(pfausse2[1:415],retard2[1:415], xlabel="pfa", ylabel ="retard",  label="daGFSS", w=2, dpi=300)
plot!(pfausse3[10:400], retard3[10:400], label="daGFSS Independent", w=2, dpi=300)
plot!(pfausse[1:378], retard[1:378], label="daGFSS with 2-norm", w=2, dpi=300)
plot!(pfausse4[60:480], retard4[60:480], label="daGFSS Centralized", w=2, dpi=300)
