using CSV
using DataFrames
using StatsBase


const timeLimit = 600
fd = "../results/robot_results_1/"
folders = readdir("$(fd)")

for folder in folders
    folderName = string(fd, folder, "/")

    csvFiles = readdir("$(folderName)")
    methodList = [:pwl, :pwlMerge, :pwr]
    pwlMethods = Dict(:pwl=>[:Incremental, :CC, :MC, :Logarithmic, :ZigZag, :ZigZagInteger],
                          :pwlMerge=>[:Incremental, :Logarithmic, :SOS2, :ZigZag, :ZigZagInteger],
                          :pwr=>[:BRGC, :balancedCode, :biclique])

    result = Dict()
    for method in methodList
        for approach in pwlMethods[method]
            result[(method, approach)] = Dict("arr"=>[], "win"=>0, "fail"=>0)
        end
    end

    for csvName in csvFiles
        df = CSV.read(string(folderName, csvName), DataFrame)
        m, n = size(df)
        winner = (:None, :None)
        usedTime = timeLimit * 2
        for i in 1:m
            method = Symbol(df[i, "Methods"])
            approach = Symbol(df[i, "Approaches"])
            time = df[i, "GurobiTimes"]
            append!(result[(method, approach)]["arr"], time)
            if (time < usedTime)
                winner = (method, approach)
                usedTime = time
            end
            if (time >= timeLimit)
                result[(method,approach)]["fail"] += 1
            end
        end
        result[winner]["win"] += 1
    end

    for method in methodList
        for approach in pwlMethods[method]
            println("==============================")
            println(folderName)
            println("Method ", method, " Approach ", approach)
            d = result[(method, approach)]
            println("Win ", d["win"], " Fail ", d["fail"])
            println("Mean ", mean(d["arr"]), " Std ", StatsBase.std(d["arr"]))
            println("==============================")
        end
    end
end
