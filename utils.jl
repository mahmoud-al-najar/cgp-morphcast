using JSON

function min_max_norm(arr, min, max)
    b = arr
    for i in 1:length(arr)
        b[i] = (arr[i] - min) / (max - min)
    end
    b
end

function get_dataset_path()
    JSON.parsefile("paths.json")["dataset_path"]
end