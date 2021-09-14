using MAT
using Dates
using GLPK
using JuMP
using Statistics
using CartesianGeneticProgramming
using Cambrian
import Cambrian.mutate
using Plots
using Plots.PlotMeasures
using StatsPlots
using DataFrames
using GLM


function get_wave_t_in_range(param, min, max)
    indices_more = findall(>=(min), param)
    indices_less = findall(<(max), param)

    indices_to_keep = intersect(indices_more, indices_less)
    param = param[indices_to_keep]

    param, indices_to_keep
end

function fit_poly_inf_norm(x,y,n)
	m = Model(GLPK.Optimizer)

    @variable(m, alpha[1:n+1])
	@variable(m, t >= 0.)
	@constraint(m, [i=1:length(y)], sum(alpha[j]*x[i]^(j-1) for j in 1:n+1) - y[i] <= t);
	@constraint(m, [i=1:length(y)], -t <= sum(alpha[j]*x[i]^(j-1) for j in 1:n+1) - y[i]);
	@objective(m, Min, t)

	optimize!(m) # added new 2020-04-08

    return JuMP.value.(alpha)
end

f2(x, args) = (args[2] * x) .+ args[1]

const MATLAB_EPOCH = Dates.DateTime(-1,12,31)

date2num(d::Dates.DateTime) = Dates.value(d-MATLAB_EPOCH)/(1000*60*60*24)
num2date(n::Number) =  MATLAB_EPOCH + Dates.Millisecond(round(Int64, n*1000*60*60*24))

function plot_bathy(xs, profile; show_points=false, x_last_positive=NaN, x_first_negative=NaN, value_last_positive=NaN, value_first_negative=NaN, x_at_zero=NaN, title=NaN)
    p1 = plot(xs, profile, label="profile_t0")
    if show_points
        scatter!([x_last_positive, x_first_negative], [value_last_positive, value_first_negative], label="A, B")
        scatter!([x_at_zero], [0], label="Shoreline")
    end
    display(plot(p1, title=title))
end

function find_x_at_zero(profile, xs, i)
    t0_index_last_positive = findlast(profile.>0)
    if t0_index_last_positive == nothing t0_index_last_positive = 1 end

    t0_index_first_negative = findfirst(profile.<0)
    if t0_index_first_negative == nothing t0_index_first_negative = 1 end

    t0_x_last_positive = xs[t0_index_last_positive]
    t0_x_first_negative = xs[t0_index_first_negative]
    t0_value_last_positive = profile[t0_index_last_positive]
    t0_value_first_negative = profile[t0_index_first_negative]

    t0_coeffs = fit_poly_inf_norm([t0_value_last_positive, t0_value_first_negative], [t0_x_last_positive, t0_x_first_negative], 1)
    t0_x_at_zero = f2([0], t0_coeffs)[1]
    t0_x_at_zero
end

dir_path = "/media/mn/WD4TB/topo/grand_popo/"
kalman_bathy_vars = matread(string(dir_path, "Kalman_bathy.mat"))
bathy_xb = kalman_bathy_vars["xb"][1, :]# 1×213 Matrix{Float64}
# bathy_tn = kalman_bathy_vars["tn"]      # 1×1214 Matrix{Float64}
bathy_hbK = kalman_bathy_vars["hbK"]    # 213×1214 Matrix{Float64}

min_target_matdate = minimum(filter(!isnan, kalman_bathy_vars["tn"]))
max_target_matdate = maximum(filter(!isnan, kalman_bathy_vars["tn"]))

wave_vars = matread(string(dir_path, "Wave_GPP_ERA.mat"))
wave_t = wave_vars["t"]
wave_Tm = wave_vars["Tm"]
wave_hs = wave_vars["hs"]
wave_dir = wave_vars["dir"]

wave_t, indices = get_wave_t_in_range(wave_t, min_target_matdate, max_target_matdate)
wave_Tm = wave_Tm[indices]
wave_hs = wave_hs[indices]
wave_dir = wave_dir[indices]

bathy_dates = map(num2date, collect(min_target_matdate:max_target_matdate))  # collect() instead of using bathy_tn because bathy_tn includes NaN's, even though bathy_hbK has no NaN
wave_dates = map(num2date, wave_t)

days = []
n_offshore = 0
n_landward = 0
n_none = 0

mean_tp = []
min_tp = []
max_tp = []
std_tp = []

mean_hs = []
min_hs = []
max_hs = []
std_hs = []

mean_dir = []
min_dir = []
max_dir = []
std_dir = []

changes = []
for i in 2:length(bathy_dates)
    _d0 = bathy_dates[i-1]
    _d1 = bathy_dates[i]

    _t0 = date2num(_d0)
    _t1 = date2num(_d1)

    _bathy_t0 = bathy_hbK[:, i-1]
    _bathy_t1 = bathy_hbK[:, i]

    _wave_t = deepcopy(wave_t)
    _wave_Tm = deepcopy(wave_Tm)
    _wave_dir = deepcopy(wave_dir)
    _wave_hs = deepcopy(wave_hs)

    _wave_t, _indices = get_wave_t_in_range(_wave_t, _t0, _t1)
    _wave_Tm = _wave_Tm[_indices]
    _wave_hs = _wave_hs[_indices]
    _wave_dir = _wave_dir[_indices]

    t0_x_at_zero = find_x_at_zero(_bathy_t0, bathy_xb, i)
    t1_x_at_zero = find_x_at_zero(_bathy_t1, bathy_xb, i)
    change = (t1_x_at_zero - t0_x_at_zero)# / 100

    # if change > 1 || change < -1
    #     println(i, "    ", change)
    #     exit()
    # end

    if change > 0
        class = (1, 0, 0)  # (offshore, landward, none)
    elseif change < 0
        class = (0, 1, 0)
    else
        class = (0, 0, 1)
    end

    # _wave_Tm = _wave_Tm ./ 100
    # _wave_hs = _wave_hs ./ 10
    # _wave_dir = _wave_dir ./ 100

    push!(mean_tp, mean(_wave_Tm))
    push!(min_tp, minimum(_wave_Tm))
    push!(max_tp, maximum(_wave_Tm))
    push!(std_tp, std(_wave_Tm))
        
    push!(mean_hs, mean(_wave_hs))
    push!(min_hs, minimum(_wave_hs))
    push!(max_hs, maximum(_wave_hs))
    push!(std_hs, std(_wave_hs))

    push!(mean_dir, mean(_wave_dir))
    push!(min_dir, minimum(_wave_dir))
    push!(max_dir, maximum(_wave_dir))
    push!(std_dir, std(_wave_dir))

    push!(changes, change)
end

train_df = DataFrame(
    Tp_mean=convert(Vector{Float64}, mean_tp[1:800]), 
    Tp_min=convert(Vector{Float64}, min_tp[1:800]), 
    Tp_max=convert(Vector{Float64}, max_tp[1:800]), 
    Tp_std=convert(Vector{Float64}, std_tp[1:800]), 
    
    Hs_mean=convert(Vector{Float64}, mean_hs[1:800]), 
    Hs_min=convert(Vector{Float64}, min_hs[1:800]), 
    Hs_max=convert(Vector{Float64}, max_hs[1:800]), 
    Hs_std=convert(Vector{Float64}, std_hs[1:800]), 
    
    Dir_mean=convert(Vector{Float64}, mean_dir[1:800]), 
    Dir_min=convert(Vector{Float64}, min_dir[1:800]), 
    Dir_max=convert(Vector{Float64}, max_dir[1:800]), 
    Dir_std=convert(Vector{Float64}, std_dir[1:800]),

    Change=convert(Vector{Float64}, changes[1:800])
)

test_df = DataFrame(
    Tp_mean=convert(Vector{Float64}, mean_tp[800:end]), 
    Tp_min=convert(Vector{Float64}, min_tp[800:end]), 
    Tp_max=convert(Vector{Float64}, max_tp[800:end]), 
    Tp_std=convert(Vector{Float64}, std_tp[800:end]), 
    
    Hs_mean=convert(Vector{Float64}, mean_hs[800:end]), 
    Hs_min=convert(Vector{Float64}, min_hs[800:end]), 
    Hs_max=convert(Vector{Float64}, max_hs[800:end]), 
    Hs_std=convert(Vector{Float64}, std_hs[800:end]), 
    
    Dir_mean=convert(Vector{Float64}, mean_dir[800:end]), 
    Dir_min=convert(Vector{Float64}, min_dir[800:end]), 
    Dir_max=convert(Vector{Float64}, max_dir[800:end]), 
    Dir_std=convert(Vector{Float64}, std_dir[800:end]),

    Change=convert(Vector{Float64}, changes[800:end])
)

println(describe(train_df))

model = lm(@formula(Change ~ Tp_mean + Tp_min + Tp_max + Tp_std + Hs_mean + Hs_min + Hs_max + Hs_std + Dir_mean + Dir_min + Dir_max + Dir_std), train_df)
println(model)
println("R2 score: ", r2(model))

ypredicted_test = predict(model, test_df)
ypredicted_train = predict(model, train_df)

rmse(target, preds) = mean((target .- preds).^2)

println("Train RMSE: ", rmse(changes[1:800], ypredicted_train))
println("Test RMSE: ", rmse(changes[800:end], ypredicted_test))

corr_train = round(cor(changes[1:800], ypredicted_train), digits=2)
corr_test = round(cor(changes[800:end], ypredicted_test), digits=2)
println("Train corr: ", corr_train)
println("Test corr: ", corr_test)

change_xs = collect(1:length(changes))
train_xs = collect(1:800)
test_xs = collect(800:length(changes))

p1 = plot(change_xs, changes, label="Raw data")
plot!(train_xs, ypredicted_train, label="Train set preds")
plot!(test_xs, ypredicted_test, label="Test set preds")
# display(plot(p1))
# readline()

plot!(size=(1920, 1080))
savefig("/home/mn/JuliaProjects/old_cgp-morphcast/plots/grand_popo/linear_regression.pdf")
