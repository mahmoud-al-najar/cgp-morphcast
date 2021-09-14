using MAT
using Dates
using GLPK
using JuMP
using Plots
using Plots.PlotMeasures


function get_wave_t_in_range(param, min, max)
    indices_to_keep = findall(>=(min), param)
    param = param[indices_to_keep]
    
    indices_to_keep = findall(<(max), param)
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
    if t0_index_last_positive == nothing
        println(i, ": t0_index_last_positive")
        t0_index_last_positive = 1
    end
    # if t0_index_last_positive == nothing t0_index_last_positive = 1 end

    t0_index_first_negative = findfirst(profile.<0)
    if t0_index_first_negative == nothing
        println(i, ": t0_index_first_negative") 
        t0_index_first_negative = 1 
    end

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
    push!(changes, t1_x_at_zero - t0_x_at_zero)

end

# p1 = plot(bathy_dates[2:end], changes, ylabel="Shoreline change", legend=false)
# p2 = plot(wave_dates, wave_Tm, ylabel="Tm", legend=false)
# p3 = plot(wave_dates, wave_hs, ylabel="Hs", legend=false)
# p4 = plot(wave_dates, wave_dir, ylabel="Dir", legend=false)

# full_plot = plot(p1, p2, p3, p4, layout=(4, 1), legend=false)#, title="Grand Popo shoreline change")
# plot!(size=(1920, 1080), left_margin=10mm) # 1920, 1080  # 1280, 720
# savefig("/home/mn/JuliaProjects/old_cgp-morphcast/plots/grand_popo/wavesandshorelinechange.pdf")

