using Plots
using DelimitedFiles

pfa = readdlm("data/pfausse_delay.csv")
pfa2 = readdlm("data/pfausse2_delay.csv")
pfa3 = readdlm("data/pfausse3_delay.csv")
pfa4 = readdlm("data/pfausse4_delay.csv")

pdetect = readdlm("data/pdetect_delay.csv")
pdetect2 = readdlm("data/pdetect2_delay.csv")
pdetect3 = readdlm("data/pdetect3_delay.csv")
pdetect4 = readdlm("data/pdetect4_delay.csv")

delay = readdlm("data/delay_var.csv")
delay2 = readdlm("data/delay2_var.csv")
delay3 = readdlm("data/delay3_var.csv")
delay4 = readdlm("data/delay4_var.csv")


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
