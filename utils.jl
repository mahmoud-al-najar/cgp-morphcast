using JSON

function norm(arr, min, max)
    b = arr
    for i in 1:length(arr)
        b[i] = (arr[i] - min) / (max - min)
    end
    b
end

function denorm(arr, min, max)
    b = arr
    for i in 1:length(arr)
        b[i] = arr[i] * (max - min) + min
    end
    b
end

function get_dataset_path()
    JSON.parsefile("paths.json")["dataset_path"]
end
