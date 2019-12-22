`julia` code for:
> "Distributed Change Detection in Streaming Graph Signals", <br />
> A. Ferrari, C. Richard and L. Verduci, <br />
> *2019 IEEE International Workshop on Computational Advances in Multi-Sensor Adaptive Processing (CAMSAP)*.


- Figure 2: `plot_signals.jl` plots the signals on the different clusters and aGFSS detector output.

- Figure 3: `synth_gfss_parallel.jl` computes and plots the approximated filter response. 
- Figures 4 and 5:  `Monte_Carlo_ROC.jl` computes and saves the probabilities of false alarm, probabilities of detection and mean delay for each algorithm. Monte-Carlo simulation uses `julia` distributed facilities. `plot_ROC.jl`: plots the ROCs and the mean detection delays.
- Figure 6: `add_sig_gephi.jl` creates the .csv to be imported in gephi. 
- `test_variance.jl`: check the analytical expression of the variance in Eq. (12) using Monte-Carlo simulations.

All the functions are included in `gfss_func.jl`.
