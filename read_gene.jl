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

gene_json = JSON.parsefile("path/to/gene")
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

    waves_path = string(dataset_path, '/', i, "_waves.csv")
    waves_csv = CSV.read(waves_path, DataFrame)

    tp_month = norm_to_range(waves_csv.tp, 0, 25, -1, 1)
    hs_month = norm_to_range(waves_csv.hs, 0, 10, -1, 1)
    direction_month = norm_to_range(waves_csv.direction, 0, 360, -1, 1)

    tp_weeks = []
    hs_weeks = []
    direction_weeks = []

    for i in range(0, stop=27)
        tp = tp_month[(i*8) + 1 : (i*8) + 8]
        hs = hs_month[(i*8) + 1 : (i*8) + 8]
        direction = direction_month[(i*8) + 1 : (i*8) + 8]
        push!(tp_weeks, tp)
        push!(hs_weeks, hs)
        push!(direction_weeks, direction)
    end
    month = Month(topo0, topo1, tp_weeks, hs_weeks, direction_weeks)
    global months
    push!(months, month)
end

month = months[1]

topo0 = deepcopy(month.topo0)
for w_i in 1:size(month.tp)[1]
    global topo0
    
    x_row = vcat(topo0, month.tp[w_i], month.hs[w_i], month.direction[w_i])
    topo0 = process(ind, x_row)
end
println("final plot:")

display(plot(norm_to_range(topo0, -1, 1, -10.0, 10.0)))
readline()
