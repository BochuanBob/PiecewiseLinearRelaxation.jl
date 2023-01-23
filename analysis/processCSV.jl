using CSV
using DataFrames
using StatsBase


methodList = ["pwl", "pwlMerge", "pwr"]
pwlMethods = Dict("pwl"=>["Incremental", "CC", "MC", "Logarithmic", "ZigZag", "ZigZagInteger"],
                  "pwlMerge"=>["Incremental", "Logarithmic", "SOS2", "ZigZag", "ZigZagInteger"],
                  "pwr"=>["BRGC", "balancedCode", "biclique"])

df = CSV.read(ARGS[1], DataFrame)

results = unique(df[!, "Experiments"])
m, n = size(df)
d = Dict()
for i in 1:m
    d[(df[i, "Methods"], df[i , "Approaches"])] = []
end
names = Vector{String}([])
for res in results
    sub_df = df[df.Experiments .== res, :]
    sub_m, sub_n = size(sub_df)
    append!(names, [res, "", "", ""])
    for i in 1:sub_m
        append!(d[(sub_df[i, "Methods"], sub_df[i, "Approaches"])], sub_df[i, "Means"])
        append!(d[(sub_df[i, "Methods"], sub_df[i, "Approaches"])], sub_df[i, "Stds"])
        append!(d[(sub_df[i, "Methods"], sub_df[i, "Approaches"])], sub_df[i, "Wins"])
        append!(d[(sub_df[i, "Methods"], sub_df[i, "Approaches"])], sub_df[i, "Fails"])
    end
end

df_res = DataFrame(names=names)

for method in methodList
    for form in pwlMethods[method]
        if (haskey(d,(method, form)))
            df_new = DataFrame(a=d[method, form])
            rename!(df_new, :a=> "$(method)_$(form)")
            global df_res = hcat(df_res, df_new)
        end
    end
end

CSV.write("$(ARGS[1])_out.csv", df_res)
