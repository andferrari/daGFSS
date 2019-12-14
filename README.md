`julia` code for:
> "Distributed Change Detection in Streaming Graph Signals", <br />
> A. Ferrari, C. Richard and L. Verduci, <br />
> *2019 IEEE International Workshop on Computational Advances in Multi-Sensor Adaptive Processing (CAMSAP)*.


- `synth_gfss_parallel.jl`: compute approximated filter response
- `affichage_signaux.jl`: plot signals on different clusters and detectors outputs
- `1000_signaux_comparaison_CORE.jl`: compute pfa, pd and delay for each algorithm from 1000 signals
- `detection_func.jl`: functions to detect and compute t_change, n_change and cluster_change
- `gfss_func.jl`: functions 
- `detection_voisin.jl`: test daGFSS algo with 1 signal
- `performance_func.jl`: compute pfa, pd
- `signaux_func.jl`: generate signals
- `add_sig_gephi.jl`: create the csv to be imported in gephi
- `plot_CORE.jl`: plot pd and mean delay as a function of pfa
