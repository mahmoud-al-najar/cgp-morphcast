using CSV
using Plots
using DataFrames

headers = ["timestamp", "cambrian", "info", "gen", "max", "avg", "std"]
df = CSV.read("path/to/csv", DataFrame, header=headers)

plot(df.gen, df.avg, xlabel="Generation", ylabel="Negative fitness", label="label")
