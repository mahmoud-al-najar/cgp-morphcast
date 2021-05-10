using JSON

function norm(arr, min, max)
    (arr .- min) ./ (max - min)
end

function denorm(arr, min, max)
    arr .* (max - min) + min
end

function norm_to_range(arr, min_x, max_x, range_a, range_b)
    (range_b - range_a) .* ((arr .- min_x) ./ (max_x - min_x)) .+ range_a
end

function get_dataset_path()
    JSON.parsefile("paths.json")["dataset_path"]
end
