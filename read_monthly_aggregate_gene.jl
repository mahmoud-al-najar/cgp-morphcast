using NPZ
using CSV
using JSON
using Plots
using DataFrames
using Statistics
using CartesianGeneticProgramming

include("utils.jl")

struct Month
    topo0
    topo1
    tp
    hs
    direction
end

gene_json = JSON.parsefile("/home/mn/JuliaProjects/cgp-morphcast/gens/2021-05-19T15:31:14.968/23000/0010.dna")
fitness = gene_json["fitness"]
functions = gene_json["fs"]
chromosome = convert(Array{Float64,1}, gene_json["chromosome"])
xs = gene_json["xs"]
ys = gene_json["ys"]

ind = CGPInd(get_config("topo.yaml"), chromosome)

##########################################

dataset_path = get_dataset_path()
n_months = 236
months = []

for i in 1:n_months
    topo_path = string(dataset_path, '/', i, "_topo.npy")
    topo = npzread(topo_path)
    topo0 = norm_to_range(topo[1, :], -10.0, 10.0, -1, 1)
    topo1 = norm_to_range(topo[2, :], -10.0, 10.0, -1, 1)
    d_topo = topo1 - topo0

    @assert topo0 + d_topo == topo1

    waves_path = string(dataset_path, '/', i, "_waves.csv")
    waves_csv = CSV.read(waves_path, DataFrame)

    tp_month = norm_to_range(waves_csv.tp, 0, 25, -1, 1)
    hs_month = norm_to_range(waves_csv.hs, 0, 10, -1, 1)
    direction_month = norm_to_range(waves_csv.direction, 0, 360, -1, 1)

    input_row = vcat(topo0, 
        mean(tp_month), std(tp_month), maximum(tp_month), median(tp_month), 
        mean(hs_month), std(hs_month), maximum(hs_month), median(hs_month), 
        mean(direction_month), std(direction_month), maximum(direction_month), median(direction_month))
    pred_delta_topo = process(ind, input_row)

    println("final plot:")

    plot_d_topo = plot(d_topo, label="Target Δtopo")
    plot!(plot_d_topo, pred_delta_topo, label="Predicted Δtopo")

    plot_target = plot(norm_to_range(topo0 + pred_delta_topo, -1, 1, -10.0, 10.0))

    display(plot(plot_d_topo, plot_target, layout=(2, 1)))
    readline()
end
