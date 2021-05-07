using NPZ
using CSV
using DataFrames
using Statistics
using CartesianGeneticProgramming
using Cambrian
import Cambrian.mutate
using Plots
include("utils.jl")


dataset_path = get_dataset_path()
dataset_size = 236

X = zeros(748, dataset_size)
Y = zeros(76, dataset_size)

for i in 1:dataset_size
    global X
    topo_path = string(dataset_path, '/', i, "_topo.npy")
    topo = npzread(topo_path)
    topo0 = norm(topo[1, :], -8.0, 8.0)
    topo1 = norm(topo[2, :], -8.0, 8.0)

    waves_path = string(dataset_path, '/', i, "_waves.csv")
    waves_csv = CSV.read(waves_path, DataFrame)

    # t = min_max_norm(waves_csv.t)
    tp = norm(waves_csv.tp, 0, 25)
    hs = norm(waves_csv.hs, 0, 10)
    direction = norm(waves_csv.direction, 0, 360)

    x_row = vcat(topo0, tp, hs, direction)
    y_row = topo1

    X[:, i] = x_row
    Y[:, i] = y_row
end



function pass_month(ind, month_x)
    topo0 = month_x[1:76]
    
    shift = 76
    tp_month = month_x[shift:shift+224]

    shift = shift + 224
    hs_month = month_x[shift:shift+224]

    shift = shift + 224
    direction_month = month_x[shift:shift+224]

    out = []

    for i in range(0, stop=27)
        tp = tp_month[(i*8) + 1 : (i*8) + 8 + 1]
        hs = hs_month[(i*8) + 1 : (i*8) + 8 + 1]
        direction = direction_month[(i*8) + 1 : (i*8) + 8 + 1]

        x_row = vcat(topo0, tp, hs, direction)
        out = process(ind, x_row)
    end
    out
end

function evaluate(ind::CGPInd, X::AbstractArray, Y::AbstractArray)
    out = pass_month(ind, X[:, 1])
    out = denorm(out, -8, +8)
    
    errors = out - Y[:, 1]
    errors = errors .^ 2
    mse = mean(errors)
    mean_mse = mse

    for i in 2:size(X, 2)
        out = pass_month(ind, X[:, i])
        out = denorm(out, -8, +8)
        y = denorm(Y[:, i], -8, +8)
        errors = out - y
        errors = errors .^ 2
        mse = mean(errors)
        if abs(mse) == Inf || isnan(mse)
            mse = 500
        end
        mean_mse = (mean_mse + mse) / 2.0
    end
    neg_rmse = sqrt(mean_mse) * -1
    [neg_rmse]
end

cfg = get_config("topo.yaml")
fit(i::CGPInd) = evaluate(i, X, Y)
mutate(i::CGPInd) = goldman_mutate(cfg, i)
e = CGPEvolution(cfg, fit)
run!(e)
