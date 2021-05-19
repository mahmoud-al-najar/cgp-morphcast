using NPZ
using CSV
using Dates
using DataFrames
using Statistics
using CartesianGeneticProgramming
using Cambrian
import Cambrian.mutate
using Plots
using StatsPlots
include("utils.jl")


struct Month
    topo0
    topo1
    d_topo
    
    tp
    tp_mean
    tp_std
    tp_max
    tp_med
    
    hs
    hs_mean
    hs_std
    hs_max
    hs_med
    
    direction
    dir_mean
    dir_std
    dir_max
    dir_med

    input_row
end

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

    # t = min_max_norm(waves_csv.t)
    tp_month = norm_to_range(waves_csv.tp, 0, 25, -1, 1)
    hs_month = norm_to_range(waves_csv.hs, 0, 10, -1, 1)
    direction_month = norm_to_range(waves_csv.direction, 0, 360, -1, 1)

    # plot_tp = density(tp_month, label="Tp")
    # plot_hs = density(hs_month, label="Hs")
    # plot_dir = density(direction_month, label="Direction")
    # display(plot(plot_tp, plot_hs, plot_dir, layout = (3, 1)))
    # readline()
    
    input_row = vcat(topo0, 
        mean(tp_month), std(tp_month), maximum(tp_month), median(tp_month), 
        mean(hs_month), std(hs_month), maximum(hs_month), median(hs_month), 
        mean(direction_month), std(direction_month), maximum(direction_month), median(direction_month))

    month = Month(topo0, topo1, topo1 - topo0, 
        tp_month, mean(tp_month), std(tp_month), maximum(tp_month), median(tp_month), 
        hs_month, mean(hs_month), std(hs_month), maximum(hs_month), median(hs_month),
        direction_month, mean(direction_month), std(direction_month), maximum(direction_month), median(direction_month),
        input_row)
    push!(months, month)
end

function evaluate(ind::CGPInd, months::Vector{Any})    
    errors = Vector{Float64}(undef, size(months)[1])

    for m_i in 1:size(months)[1]
        month = months[m_i]        
        d_topo = process(ind, month.input_row)
        errors[m_i] = mean(abs.(d_topo - month.d_topo))
    end
    [mean(errors) * -1]
end

cfg = get_config("topo.yaml")
fit(i::CGPInd) = evaluate(i, months)
mutate(i::CGPInd) = goldman_mutate(cfg, i)
e = CGPEvolution(cfg, fit)
run!(e)
