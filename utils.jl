using JSON

get_dataset_path(;dataset_name="dataset_path") = JSON.parsefile("paths.json")[dataset_name]

norm(arr, min, max) = (arr .- min) ./ (max - min)
denorm(arr, min, max) = arr .* (max - min) + min
norm_to_range(arr, min_x, max_x, range_a, range_b) = (range_b - range_a) .* ((arr .- min_x) ./ (max_x - min_x)) .+ range_a

has_inf(arr) = sum(isinf.(arr)) > 0
has_nan(arr) = sum(isnan.(arr)) > 0
