using Plots
using DelimitedFiles

pfa = readdlm("donnees/pfausseret.csv")
pfa2 = readdlm("donnees/pfausse2ret.csv")
pfa3 = readdlm("donnees/pfausse3ret.csv")
pfa4 = readdlm("donnees/pfausse4ret.csv")

pdetect = readdlm("donnees/pdetectret.csv")
pdetect2 = readdlm("donnees/pdetect2ret.csv")
pdetect3 = readdlm("donnees/pdetect3ret.csv")
pdetect4 = readdlm("donnees/pdetect4ret.csv")

delay = readdlm("donnees/retardvar.csv")
delay2 = readdlm("donnees/retard2var.csv")
delay3 = readdlm("donnees/retard3var.csv")
delay4 = readdlm("donnees/retard4var.csv")


pyplot()

Plots.reset_defaults()
Plots.scalefontsizes(1.5)

plot(pfa2, pdetect2, label="daGFSS", w=2, dpi=300)
plot!(pfa, pdetect, label="daGFSS with 2-norm", w=2, dpi=300)
plot!(pfa4, pdetect4, ratio=:equal, xlabel="Pfa", ylabel="Pd", label="daGFSS Centralized", w=2, dpi=300)
plot!(pfa3, pdetect3, label="daGFSS Independent", w=2, dpi=300)
savefig("newCore.pdf")

Plots.reset_defaults()
Plots.scalefontsizes(1.5)

plot(pfa2,delay2, xlabel="Pfa", ylabel ="Mean detection delay",  label="daGFSS", w=2, dpi=300)
plot!(pfa3, delay3, label="daGFSS Independent", w=2, dpi=300)
plot!(pfa, delay, label="daGFSS with 2-norm", w=2, dpi=300)
plot!(pfa4[1:800], delay4[1:800], label="daGFSS Centralized", w=2, dpi=300, ylims=(7,23))

savefig("mean_detection_delay.pdf")
