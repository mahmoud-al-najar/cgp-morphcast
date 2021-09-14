using JuMP
using GLPK
using Plots  #; plotlyjs()
using Measures
using DataFrames
using NPZ
using CSV
using Dates
using DataFrames
using Statistics
using CartesianGeneticProgramming
using Cambrian
import Cambrian.mutate
using Plots, Measures
using StatsPlots
include("../utils.jl")
using Dierckx
using BSplineKit
using Distances


struct Month
    topo0
    topo1

    tp
    hs
    direction

    shoreline_change
    shoreline_change_class  # (landward, seaward)
end

function fit_poly_inf_norm(x,y,n)
	m = Model(GLPK.Optimizer)

    @variable(m, alpha[1:n+1])
	@variable(m, t >= 0.)
	@constraint(m, [i=1:length(y)], sum(alpha[j]*x[i]^(j-1) for j in 1:n+1) - y[i] <= t);
	@constraint(m, [i=1:length(y)], -t <= sum(alpha[j]*x[i]^(j-1) for j in 1:n+1) - y[i]);
	@objective(m, Min, t)

	optimize!(m)

    return JuMP.value.(alpha)
end


f2(x, args) = (args[2] * x) .+ args[1]

dataset_path = get_dataset_path(;dataset_name="series_dataset_path")
n_series = 6
series = []
count_months = 0
xs_csv_path = get_dataset_path(;dataset_name="xs_csv_path")
xs_csv = CSV.read(xs_csv_path, DataFrame)
xs = xs_csv.x

changes = []

for s in 1:n_series
    push!(series, [])
    series_dir = joinpath(dataset_path, string(s))

    files = readdir(series_dir)
    topo_files = filter(f-> !occursin("_waves.csv", f), files)
    topo_files = map(f->replace(f, "_topo.npy" => ""), topo_files)

    month_numbers = sort(map(f->parse(Int64, f), topo_files))
    n_months = length(month_numbers)
    global count_months += n_months

    println(s, " - ", n_months)

    for mi in 1:n_months
        i = month_numbers[mi]
        topo_path = string(dataset_path, '/', s, '/', i, "_topo.npy")
        if i == 312
            println("TEST")
            exit()
            continue
        end
        topo = npzread(topo_path)

        topo0 = topo[1, :]
        topo1 = topo[2, :]

        t0_index_last_positive = findlast(topo0.>0)
        t0_index_first_negative = findfirst(topo0.<0)

        t0_x_last_positive = xs[t0_index_last_positive]
        t0_x_first_negative = xs[t0_index_first_negative]
        t0_value_last_positive = topo0[t0_index_last_positive]
        t0_value_first_negative = topo0[t0_index_first_negative]

        t0_coeffs = fit_poly_inf_norm([t0_value_last_positive, t0_value_first_negative], [t0_x_last_positive, t0_x_first_negative], 1)
        t0_x_at_zero = f2([0], t0_coeffs)[1]
        
        t1_index_last_positive = findlast(topo1.>0)
        t1_index_first_negative = findfirst(topo1.<0)

        t1_x_last_positive = xs[t1_index_last_positive]
        t1_x_first_negative = xs[t1_index_first_negative]
        t1_value_last_positive = topo1[t1_index_last_positive]
        t1_value_first_negative = topo1[t1_index_first_negative]

        t1_coeffs = fit_poly_inf_norm([t1_value_last_positive, t1_value_first_negative], [t1_x_last_positive, t1_x_first_negative], 1)
        t1_x_at_zero = f2([0], t1_coeffs)[1]

        shoreline_change = t1_x_at_zero - t0_x_at_zero
        shoreline_change = norm_to_range(shoreline_change, -15.0, 15.0, -1.0, 1.0)
        if shoreline_change < 0
            class = (0, 1)  # (landward, seaward)
        else
            class = (1, 0)  # (landward, seaward)
        end

        waves_path = string(dataset_path, '/', s, '/', i, "_waves.csv")
        waves_csv = CSV.read(waves_path, DataFrame)

        tp_month = norm_to_range(waves_csv.tp, 0.0, 25.0, -1.0, 1.0)
        hs_month = norm_to_range(waves_csv.hs, 0.0, 10.0, -1.0, 1.0)
        direction_month = norm_to_range(waves_csv.direction, 0.0, 360.0, -1.0, 1.0)

        month = Month(topo0, topo1, tp_month, hs_month, direction_month, shoreline_change, class)
        push!(series[end], month)
    end
end

real_class = []
real_change = []
tp_min = []
tp_max = []
tp_avg = []
hs_min = []
hs_max = []
hs_avg = []
dir_min = []
dir_max = []
dir_avg = []

for si in 1:length(series)
    global correct
    global wrong
    months = series[si]
    for mi in 1:length(months)
        month = months[mi]
        push!(real_class, month.shoreline_change_class)
        push!(real_change, norm_to_range(month.shoreline_change, -1.0, 1.0, -15.0, 15.0))
        push!(tp_min, minimum(month.tp))
        push!(tp_max, maximum(month.tp))
        push!(tp_avg, mean(month.tp))
        push!(hs_min, minimum(month.hs))
        push!(hs_max, maximum(month.hs))
        push!(hs_avg, mean(month.hs))
        push!(dir_min, minimum(month.direction))
        push!(dir_max, maximum(month.direction))
        push!(dir_avg, mean(month.direction))
    end
end

println("Max change: ", maximum(real_change))
println("Min change: ", minimum(real_change))

# p_tp_min = scatter(tp_min, preds, xlabel="min(Tp)", ylabel="Pred (1=True, 0=False)")
# p_tp_max = scatter(tp_max, preds, xlabel="max(Tp)", ylabel="Pred (1=True, 0=False)")
# p_tp_avg = scatter(tp_avg, preds, xlabel="mean(Tp)", ylabel="Pred (1=True, 0=False)")
# tp_plot = plot(p_tp_min, p_tp_max, p_tp_avg, layout=(1, 3))

# p_hs_min = scatter(hs_min, preds, xlabel="min(Hs)", ylabel="Pred (1=True, 0=False)")
# p_hs_max = scatter(hs_max, preds, xlabel="max(Hs)", ylabel="Pred (1=True, 0=False)")
# p_hs_avg = scatter(hs_avg, preds, xlabel="mean(Hs)", ylabel="Pred (1=True, 0=False)")
# hs_plot = plot(p_hs_min, p_hs_max, p_hs_avg, layout=(1, 3))

# p_dir_min = scatter(dir_min, preds, xlabel="min(Dir)", ylabel="Pred (1=True, 0=False)")
# p_dir_max = scatter(dir_max, preds, xlabel="max(Dir)", ylabel="Pred (1=True, 0=False)")
# p_dir_avg = scatter(dir_avg, preds, xlabel="mean(Dir)", ylabel="Pred (1=True, 0=False)")
# dir_plot = plot(p_dir_min, p_dir_max, p_dir_avg, layout=(1, 3))

# display(plot(tp_plot, hs_plot, dir_plot, layout=(3, 1), size = (1600, 800)))

# display(plot(real_change, legend=false, ylabel="Shoreline change [m]", xlabel="Months", size=(1280, 300)))
# readline()
# exit()


################## dx correlations ##################
p_tp_min = scatter(tp_min, real_change, xlabel="min(Tp)", ylabel="Real change", title=string("Correlation: ", round(cor(tp_min, real_change), digits=2)))
p_tp_max = scatter(tp_max, real_change, xlabel="max(Tp)", ylabel="Real change", title=string("Correlation: ", round(cor(tp_max, real_change), digits=2)))
p_tp_avg = scatter(tp_avg, real_change, xlabel="mean(Tp)", ylabel="Real change", title=string("Correlation: ", round(cor(tp_avg, real_change), digits=2)))
tp_plot = plot(p_tp_min, p_tp_max, p_tp_avg, layout=(1, 3), bottom_margin=20mm)

p_hs_min = scatter(hs_min, real_change, xlabel="min(Hs)", ylabel="Real change", title=string("Correlation: ", round(cor(hs_min, real_change), digits=2)))
p_hs_max = scatter(hs_max, real_change, xlabel="max(Hs)", ylabel="Real change", title=string("Correlation: ", round(cor(hs_max, real_change), digits=2)))
p_hs_avg = scatter(hs_avg, real_change, xlabel="mean(Hs)", ylabel="Real change", title=string("Correlation: ", round(cor(hs_avg, real_change), digits=2)))
hs_plot = plot(p_hs_min, p_hs_max, p_hs_avg, layout=(1, 3), bottom_margin=20mm)

p_dir_min = scatter(dir_min, real_change, xlabel="min(dir)", ylabel="Real change", title=string("Correlation: ", round(cor(dir_min, real_change), digits=2)))
p_dir_max = scatter(dir_max, real_change, xlabel="max(dir)", ylabel="Real change", title=string("Correlation: ", round(cor(dir_max, real_change), digits=2)))
p_dir_avg = scatter(dir_avg, real_change, xlabel="mean(dir)", ylabel="Real change", title=string("Correlation: ", round(cor(dir_avg, real_change), digits=2)))
dir_plot = plot(p_dir_min, p_dir_max, p_dir_avg, layout=(1, 3))

dx_correlations_plot = plot(tp_plot, hs_plot, dir_plot, layout=(3, 1), size = (1600, 800), margin=10mm)

################## tp correlations ##################
p_hs_min = scatter(hs_min, tp_min, xlabel="min(Hs)", ylabel="min(Tp)", title=string("Correlation: ", round(cor(hs_min, tp_min), digits=2)))
p_hs_max = scatter(hs_max, tp_max, xlabel="max(Hs)", ylabel="max(Tp)", title=string("Correlation: ", round(cor(hs_max, tp_max), digits=2)))
p_hs_avg = scatter(hs_avg, tp_avg, xlabel="mean(Hs)", ylabel="mean(Tp)", title=string("Correlation: ", round(cor(hs_avg, tp_avg), digits=2)))
hs_plot = plot(p_hs_min, p_hs_max, p_hs_avg, layout=(1, 3), bottom_margin=20mm)

p_dir_min = scatter(dir_min, tp_min, xlabel="min(dir)", ylabel="min(Tp)", title=string("Correlation: ", round(cor(dir_min, tp_min), digits=2)))
p_dir_max = scatter(dir_max, tp_max, xlabel="max(dir)", ylabel="max(Tp)", title=string("Correlation: ", round(cor(dir_max, tp_max), digits=2)))
p_dir_avg = scatter(dir_avg, tp_avg, xlabel="mean(dir)", ylabel="mean(Tp)", title=string("Correlation: ", round(cor(dir_avg, tp_avg), digits=2)))
dir_plot = plot(p_dir_min, p_dir_max, p_dir_avg, layout=(1, 3))

tp_correlations_plot = plot(hs_plot, dir_plot, layout=(2, 1), size = (1600, 600), margin=10mm)
display(plot(tp_correlations_plot))
readline()
