using CSV
using DataFrames
using StatsBase


const timeLimit = 600

function main()
    expName = "shareOfChoiceNLP_results_2"

    fd = "../results/$(expName)/result/"

    csvFiles = readdir("$(fd)")
    nums = length(csvFiles)
    res = zeros(0,0)
    print(csvFiles)
    for csvName in csvFiles
        df = CSV.read(string(fd, csvName), DataFrame)

        m, n = size(df)
        arr = zeros(1, 1+m * 3)
        for i in 1:m
            arr[1, 1] = df[i, "Indices"]
            arr[1, ((i-1)*3+2):(3*i+1)] = [df[i, "Objectives"],
                                        df[i, "Bounds"],
                                        df[i, "SolverTimes"]]
        end
        if (res == zeros(0,0))
            res = arr
        else
            res = vcat(res, arr)
        end
    end

    df_res = DataFrame(res, :auto)
    sort!(df_res, [:x1])
    CSV.write("$(expName)_out.csv", df_res)
end

main()
