using CSV
using DataFrames
using StatsBase


const timeLimit = 1800
expName = "shareOfChoice_results_2"
fd = "../results/$(expName)/"
folders = readdir("$(fd)")

expOut, methodOut, approachOut, winOut,
        failOut, meanOut, stdOut = [], [], [], [], [], [], []

for folder in folders
    folderName = string(fd, folder, "/")

    csvFiles = readdir("$(folderName)")
    methodList = [:pwlMerge, :pwr]
    pwlMethods = Dict(:pwlMerge=>[:Incremental, :DLog, :Logarithmic, :SOS2, :ZigZag, :ZigZagInteger],
                      :pwr=>[:DLog, :BRGC, :balancedCode, :biclique])

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
            println(folder)
            println("Method ", method, " Approach ", approach)
            d = result[(method, approach)]
            println("Win ", d["win"], " Fail ", d["fail"])
            println("Mean ", mean(d["arr"]), " Std ", StatsBase.std(d["arr"]))
            println("==============================")

            push!(expOut, folder)
            push!(methodOut, method)
            push!(approachOut, approach)
            append!(winOut, d["win"])
            append!(failOut, d["fail"])
            append!(meanOut, mean(d["arr"]))
            append!(stdOut, StatsBase.std(d["arr"]))
        end
    end
end

df = DataFrame(Experiments=expOut, Methods=methodOut, Approaches=approachOut,
                Wins=winOut, Fails=failOut, Means=meanOut, Stds=stdOut)
CSV.write(string(expName, ".csv"), df)
