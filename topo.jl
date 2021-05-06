using NPZ
using CSV
using DataFrames
using Statistics
using CartesianGeneticProgramming
using Cambrian
import Cambrian.mutate
include("utils.jl")


dataset_path = get_dataset_path()

X = zeros(100, 7661)
Y = zeros(76, 7661)

for i in 1:7661
    global X
    topo_path = string(dataset_path, '/', i, "_topo.npy")
    topo = npzread(topo_path)
    topo0 = min_max_norm(topo[1, :], -8.0, 8.0)
    topo1 = min_max_norm(topo[2, :], -8.0, 8.0)

    waves_path = string(dataset_path, '/', i, "_waves.csv")
    waves_csv = CSV.read(waves_path, DataFrame)

    # t = min_max_norm(waves_csv.t)
    tp = min_max_norm(waves_csv.tp, 0, 25)
    hs = min_max_norm(waves_csv.hs, 0, 10)
    direction = min_max_norm(waves_csv.direction, 0, 360)

    x_row = vcat(topo0, tp, hs, direction)
    y_row = topo1

    X[:, i] = x_row
    Y[:, i] = y_row
end

function evaluate(ind::CGPInd, X::AbstractArray, Y::AbstractArray)
    mse_list = zeros(0)

    out = process(ind, X[:, 1])
    errors = out - Y[:, 1]
    errors = errors .^ 2
    mse = mean(errors)
    mean_mse = mse

    for i in 2:500  # size(X, 2)
        out = process(ind, X[:, i])
        errors = out - Y[:, i]
        errors = errors .^ 2
        mse = mean(errors)
        if mse == -Inf || isnan(mse)
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
