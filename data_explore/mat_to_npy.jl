using MAT
using NPZ
using Dates
using GLPK
using JuMP
using Plots
using Plots.PlotMeasures
include("../utils.jl")


function get_wave_t_in_range(param, min, max)
    indices_to_keep = findall(>=(min), param)
    param = param[indices_to_keep]
    
    indices_to_keep = findall(<(max), param)
    param = param[indices_to_keep]
    param, indices_to_keep
end

dir_path = get_dataset_path(;dataset_name="grand_popo")
kalman_bathy_vars = matread(string(dir_path, "Kalman_bathy.mat"))
bathy_xb = kalman_bathy_vars["xb"][1, :]
bathy_tn = kalman_bathy_vars["tn"][1, :]
bathy_hbK = kalman_bathy_vars["hbK"]

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

outdir = "/Users/mn/JuliaProjects/cgp-morphcast/output/grand_popo/data/"

npzwrite(string(outdir, "bathy_xb.npy"), bathy_xb)
npzwrite(string(outdir, "bathy_tn.npy"), bathy_tn)
npzwrite(string(outdir, "bathy_hbK.npy"), bathy_hbK)

npzwrite(string(outdir, "wave_t.npy"), wave_t)
npzwrite(string(outdir, "wave_Tm.npy"), wave_Tm)
npzwrite(string(outdir, "wave_hs.npy"), wave_hs)
npzwrite(string(outdir, "wave_dir.npy"), wave_dir)
