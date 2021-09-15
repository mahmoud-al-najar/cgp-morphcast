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
include("../utils.jl")


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

struct Day
    bathy0
    bathy1

    tp
    hs
    direction

    shoreline_change
    shoreline_change_class
    
    input_row
end

dir_path = get_dataset_path(;dataset_name="grand_popo")
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

all_tm = []
all_hs = []
all_dir = []
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
    change = (t1_x_at_zero - t0_x_at_zero) / 100

    if change > 1 || change < -1
        println(i, "    ", change)
        exit()
    end

    if change > 0
        class = (1, 0, 0)  # (offshore, landward, none)
        global n_offshore
        n_offshore += 1
    elseif change < 0
        class = (0, 1, 0)
        global n_landward
        n_landward += 1
    else
        class = (0, 0, 1)
        global n_none
        n_none += 1
    end

    _wave_Tm = _wave_Tm ./ 100
    _wave_hs = _wave_hs ./ 10
    _wave_dir = _wave_dir ./ 100

    input_row = vcat(
        mean(_wave_Tm),
        minimum(_wave_Tm),
        maximum(_wave_Tm),
        std(_wave_Tm),
        mean(_wave_hs),
        minimum(_wave_hs),
        maximum(_wave_hs),
        std(_wave_hs),
        mean(_wave_dir),
        minimum(_wave_dir),
        maximum(_wave_dir),
        std(_wave_dir)
    )

    push!(all_tm, vcat(mean(_wave_Tm), minimum(_wave_Tm), maximum(_wave_Tm), std(_wave_Tm))...)
    push!(all_hs, vcat(mean(_wave_hs), minimum(_wave_hs), maximum(_wave_hs), std(_wave_hs))...)
    push!(all_dir, vcat(mean(_wave_dir), minimum(_wave_dir), maximum(_wave_dir), std(_wave_dir))...)

    day = Day(_bathy_t0, _bathy_t1, _wave_Tm, _wave_hs, _wave_dir, change, class, input_row)
    push!(days, day)    
end

println("normalized(all_tm)     ", minimum(all_tm), "  ", maximum(all_tm))      # RAW: 0.007926686233555585     14.506896347331828
println("normalized(all_hs)     ", minimum(all_hs), "  ", maximum(all_hs))      # RAW: 0.0007081265658494296    2.6867510628524585
println("normalized(all_dir)     ", minimum(all_dir), "  ", maximum(all_dir))   # RAW: -13.212138751720701      17.384923209385946

println("n_offshore: ", n_offshore, ", n_landward: ", n_landward, ", n_none: ", n_none)  # n_offshore: 274, n_landward: 243, n_none: 696

function evaluate(ind::CGPInd, days::Vector{Any}, cfg::NamedTuple)
    errors = Vector{Float64}(undef, length(days))
    reset!(ind)
    for di in 1:length(days)
        day = days[di]
        
        new_change = process(ind, day.input_row)[1]
        error = (new_change*100 - day.shoreline_change*100)^2

        errors[di] = error
    end
    [mean(errors) * -1]
end

function make_configs(template)
    lst_n_cols = collect(20:10:50)
    lst_n_pop = collect(40:20:200)
    lst_f_recur = collect(0:0.1:1)
    println(length(lst_n_cols))
    println(length(lst_n_pop))
    println(length(lst_f_recur))
    configs = []
    for n_pop in lst_n_pop
        n_elite = 0.1 * n_pop
        for n_cols in lst_n_cols
            for f_recur in lst_f_recur
                config = deepcopy(template)
                two_arity = falses(length(config["functions"]))
                functions = Array{Function}(undef, length(config["functions"]))
                for i in eachindex(config["functions"])
                    fname = config["functions"][i]
                    if CGPFunctions.arity[fname] == 2
                        two_arity[i] = true
                    end
                    functions[i] = eval(Meta.parse(string("CGPFunctions.", fname)))
                end
                config["two_arity"] = two_arity
                config["functions"] = functions
            
                config["n_population"] = n_pop
                config["n_elite"] = n_elite
                config["columns"] = n_cols
                config["recur"] = f_recur
                config["id"] = string("auto_", "npop", n_pop, "ncols", n_cols, "frecur", f_recur)
                push!(configs, deepcopy(config))
            end
        end
    end
    configs
end

include("cfg_template.jl")
configs = make_configs(template)
println(length(configs))

cfg = configs[1]
mutate(i::CGPInd) = goldman_mutate(cfg, i)

for config in configs
    global cfg
    cfg = (; (Symbol(k)=>v for (k, v) in config)...)
    println(cfg)
    fit(i::CGPInd) = evaluate(i, days, cfg)
    e = CGPEvolution(cfg, fit)
    run!(e)
end
