# daGFSS

- `synth_gfss_parallel.jl`: plot approximated filter response
- `affichage_signaux.jl`: plot signals on different clusters and detection's norm
- `synth_gfss_parallel.jl`: calculate poles and residues
- `1000_signaux_voisins_sans_changement.jl`: To be throw after getting poles and residues, calculate different pfa
- `1000_signaux_voisins.jl`: To be throw after "1000_signaux_voisins_sans_changement" calculate pd and plot CORE
- `detection_func.jl`: give functions that allow to detect t_change, n_change and cluster_change
- `gfss_func.jl`: give functions that calculate GFSS
- `detection_voisin.jl`: Test GSFF algo with 1 signal
- `performance_func.jl`: give functions that allow to calculate pfa, pd
- `signaux_func.jl`: give the function that create the signal we want
