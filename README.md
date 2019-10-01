`julia` code for:
> "Distributed Change Detection in Streaming Graph Signals", <br />
> A. Ferrari, C. Richard and L. Verduci, <br />
> *2019 IEEE International Workshop on Computational Advances in Multi-Sensor Adaptive Processing (CAMSAP)*.


- `synth_gfss_parallel.jl`: plot approximated filter response
- `affichage_signaux.jl`: plot signals on different clusters and detection's norm
- `synth_gfss_parallel.jl`: calculate poles and residues
- `1000_signaux_comparaison_CORE.jl`: calculate pfa, pd and delay for every algo with 1000 signals
- `detection_func.jl`: give functions that allow to detect t_change, n_change and cluster_change
- `gfss_func.jl`: give functions that calculate GFSS
- `detection_voisin.jl`: Test GSFF algo with 1 signal
- `performance_func.jl`: give functions that allow to calculate pfa, pd
- `signaux_func.jl`: give the function that create the signal we want
- `add_sig_gephi.jl`: create a csv to be imported in gephi
- `plot_CORE.jl`: plot pd and mean delay depending on pfa
