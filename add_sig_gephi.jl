using DelimitedFiles
using ColorSchemes
using LinearAlgebra
using Colors

function rescale(x)
    if norm(x) == 0
        x_resc = x
    elseif minimum(x) == maximum(x)
        x_resc = normalize(x)
    else
        x_resc = (x .- minimum(x))./(maximum(x)-minimum(x))
    end
    x_resc
end

function write_nodes_sig(vals, n_vertex, path)

    palette = ColorSchemes.viridis
    vals_n = rescale(vals)

    nodes = ["n"*string(i-1) for i in 1:n_vertex]
    colors_sig = ["#"*hex((get(palette, vals_n[i]))) for i in 1:n_vertex]

    writedlm(path, ["Id" "sig" "Color"; hcat(nodes,vals,colors_sig)], ',')

end

n_vertex = 250
vals = readdlm("donnees/t_voisin_carre_420.csv")
write_nodes_sig(vals, n_vertex, "donnees/plot_t_420.csv")
