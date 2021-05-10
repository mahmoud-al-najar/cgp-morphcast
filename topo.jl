using NPZ
using CSV
using DataFrames
using Statistics
using CartesianGeneticProgramming
using Cambrian
import Cambrian.mutate
using Plots
include("utils.jl")


struct Month
    topo0
    topo1
    tp
    hs
    direction
end

dataset_path = get_dataset_path()
n_months = 236
months = []

for i in 1:n_months
    topo_path = string(dataset_path, '/', i, "_topo.npy")
    topo = npzread(topo_path)
    topo0 = norm_to_range(topo[1, :], -8.0, 8.0, -1, 1)
    topo1 = norm_to_range(topo[2, :], -8.0, 8.0, -1, 1)

    waves_path = string(dataset_path, '/', i, "_waves.csv")
    waves_csv = CSV.read(waves_path, DataFrame)

    # t = min_max_norm(waves_csv.t)
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

function evaluate(ind::CGPInd, months::Vector{Any})    
    errors = Vector{Float64}(undef, size(months)[1])

    for m_i in 1:size(months)[1]
        month = months[m_i]
        topo0 = deepcopy(month.topo0)
        for w_i in 1:size(month.tp)[1]
            x_row = vcat(topo0, month.tp[w_i], month.hs[w_i], month.direction[w_i])
            topo0 = process(ind, x_row)
        end
        errors[m_i] = mean(abs.(month.topo1 - topo0))
    end
    [mean(errors) * -1]
end

cfg = get_config("topo.yaml")
fit(i::CGPInd) = evaluate(i, months)
mutate(i::CGPInd) = goldman_mutate(cfg, i)
e = CGPEvolution(cfg, fit)
run!(e)
