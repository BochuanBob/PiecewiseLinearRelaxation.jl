function nonlinearRelaxation!(m::JuMP.Model, x::VarOrAff, z::VarOrAff,
            lower::Float64, upper::Float64, act::Function, dAct::Function,
            points::Int64, breaks::Vector{Float64}; method=:pwl, times::Int64=1,
            pwl_method=:Logarithmic, threshold=1e-6)
    if (upper - lower <= threshold)
        @constraint(m, z == act((lower+upper)/2))
        return
    end
    xList = [lower + (i-1) * (upper - lower) / (points-1) for i = 1:points]
    for b in breaks
        if !(b in xList) && (b >= lower) && (b <= upper)
            append!(xList, b)
        end
    end
    sort!(xList)
    xList = removeClosePoints(xList)
    # println(xList)
    zList = [act(x) for x in xList]
    if (method == :pwl)
        xLow, zLow, xUpp, zUpp = getLowUppPWL(xList, act, dAct, times=times)
        zVarLow = piecewiselinear(m, x, xLow, zLow, method=pwl_method)
        zVarUpp = piecewiselinear(m, x, xUpp, zUpp, method=pwl_method)
        @constraint(m, z >= zVarLow)
        @constraint(m, z <= zVarUpp)
    elseif (method == :pwr)
        xList, zList, SS = getPWR(xList, act, dAct, times=times)
        if (pwl_method == :DLog) || (pwl_method == :Incremental)
            getDisaggregatedPWR!(xList, zList, SS)
            if (pwl_method == :DLog)
                piecewiseLinearRelaxation!(m, x, z, xList, zList, SS, method=:BRGC)
            else # pwl_method == :Incremental
                lambda = @variable(m, [1:length(xList)], lower_bound=0.0)

                @constraint(m, sum(lambda) == 1)
                @constraint(m, sum(lambda[i] * xList[i] for i in 1:length(xList)) == x)
                @constraint(m, sum(lambda[i] * zList[i] for i in 1:length(zList)) == z)

                sLen = length(SS) - 1
                u = JuMP.@variable(m, [1:(sLen-1)], lower_bound=0, upper_bound=1)
                y = @variable(m, [1:sLen], binary=true)
                for i in 1:sLen
                    @constraint(m, sum(lambda[j] for j in SS[i]) == y[i])
                end
                @constraint(m, y[1] == 1 - u[1])
                for i in 1:(sLen-2)
                    @constraint(m, u[i] >= u[i+1])
                    @constraint(m, y[i+1] == u[i] - u[i+1])
                end
                @constraint(m, y[sLen] == u[sLen-1])
            end
        else
            piecewiseLinearRelaxation!(m, x, z, xList, zList, SS, method=pwl_method)
        end
    elseif (method == :pwlMerge)
        xLow, zLow, xUpp, zUpp = getLowUppPWL(xList, act, dAct, times=times)
        xUnion, zLowU, zUppU = unionLowUppPWL(xLow, zLow, xUpp, zUpp)
        zVarLow, zVarUpp = piecewiselinear2(m, x, xUnion, zLowU, zUppU, method=pwl_method)
        @constraint(m, z >= zVarLow)
        @constraint(m, z <= zVarUpp)
    else
        error("Unrecognized method $method")
    end
end

function activation(x::Float64; name::String="sigmoid")
    if (name == "sigmoid")
        return 1 / (1 + exp(-x))
    elseif (name == "tanh")
        return (exp(x) - exp(-x)) / (exp(x) + exp(-x))
    # elseif (name == "silu")
    #     return x / (1 + exp(-x))
    elseif (name == "sin")
        return sin(x)
    elseif (name == "cos")
        return cos(x)
    end
end

sigmoid = x -> 1 / (1 + exp(-x))
dsigmoid = x -> exp(-x) / (1 + exp(-x))^2
dsin = x -> cos(x)
dcos = x -> -sin(x)
dtanh = x -> 4 / (exp(x) + exp(-x))^2

function dActivation(x::Float64; name::String="sigmoid")
    if (name == "sigmoid")
        return exp(-x) / (1 + exp(-x))^2
    elseif (name == "tanh")
        return 4 / (exp(x) + exp(-x))^2
    # elseif (name == "silu")
    #     return (1 + exp(-x) + exp(-x) * x)/ (1 + exp(-x))^2
    elseif (name == "sin")
        return cos(x)
    elseif (name == "cos")
        return -sin(x)
    end
end


function getIntersection(x1::Float64, z1::Float64, slope1::Float64,
                        x2::Float64, z2::Float64, slope2::Float64)
    # println([x1, z1, slope1, x2, z2, slope2])
    return inv([-slope1 1; -slope2 1]) * [z1 - slope1 * x1, z2 - slope2 * x2]
end

function getMultiIntersection(x1::Float64, x2::Float64, act::Function,
                        dAct::Function; times::Int64=1, check::Bool=false)
    xList = Vector{Float64}([])
    zList = Vector{Float64}([])
    intersectionHelper!(x1, x2, xList, zList, act, dAct, times, x1, x2)
    if (check)
        @assert length(xList) == length(zList)
        if (length(xList) > 0)
            zList2 = [act(x) for x in xList]
            @assert (sum(zList .<= zList2) == length(zList)) ||
                    (sum(zList .>= zList2) == length(zList))
        end
    end
    return xList, zList
end

function intersectionHelper!(left::Float64, right::Float64,
            xList::Vector{Float64}, zList::Vector{Float64},
            act::Function, dAct::Function, times::Int64,
            x1::Float64, x2::Float64; threshold::Float64=1e-8)
    if abs(dAct(left) - dAct(right)) <= threshold
        if (left != x1) && (left != x2) && !(left in xList)
            append!(xList, left)
            append!(zList, act(left))
        end

        if (right != x1) && (right != x2) && !(right in xList)
            append!(xList, right)
            append!(zList, act(right))
        end
        return
    end
    xx, zz = getIntersection(left, act(left), dAct(left),
                            right, act(right), dAct(right))
    if (times == 1)
        append!(xList, xx)
        append!(zList, zz)
    else
        intersectionHelper!(left, xx, xList, zList, act, dAct, times-1, x1, x2)
        # append!(xList, xx)
        # append!(zList, act(xx))
        intersectionHelper!(xx, right, xList, zList, act, dAct, times-1, x1, x2)
    end
    return
end

function removeLinear(xArr::Vector{Float64}, zArr::Vector{Float64}; threshold::Float64=1e-8)
    if (length(xArr) == 0)
        return Vector{Float64}([]), Vector{Float64}([])
    end
    xRes = Vector{Float64}([xArr[1]])
    zRes = Vector{Float64}([zArr[1]])
    if (length(xArr) == 0)
        return xRes, zRes
    end

    for i in 2:(length(xArr)-1)
        x1, x2, x3 = xArr[i-1], xArr[i], xArr[i+1]
        z1, z2, z3 = zArr[i-1], zArr[i], zArr[i+1]
        if abs((x1 - x2) * (z2 - z3) - (x2 - x3) * (z1 - z2)) > threshold * (x2 - x1) * (x3 - x2)
            append!(xRes, x2)
            append!(zRes, z2)
        end
    end
    append!(xRes, xArr[length(xArr)])
    append!(zRes, zArr[length(zArr)])
    return xRes, zRes
end

function getLowUppPWL(xList::Vector{Float64}, act::Function, dAct::Function;
                    times::Int64=1, threshold::Float64=1e-8)
    zList = [act(x) for x in xList]
    xLow, zLow = [x for x in xList], [z for z in zList]
    xUpp, zUpp = [x for x in xList], [z for z in zList]
    xLen = length(xList)
    for i in 1:(xLen-1)
        xx, zz = getMultiIntersection(xList[i], xList[i+1], act, dAct, times=times)
        if (dAct(xList[i]) > dAct(xList[i+1]))
            append!(xUpp, xx)
            append!(zUpp, zz)
        else
            append!(xLow, xx)
            append!(zLow, zz)
        end
    end
    p = sortperm(xLow);
    xLow = xLow[p]; zLow = zLow[p];
    xLow, zLow = removeLinear(xLow, zLow)

    p = sortperm(xUpp);
    xUpp = xUpp[p]; zUpp = zUpp[p];
    xUpp, zUpp = removeLinear(xUpp, zUpp)

    return xLow, zLow, xUpp, zUpp
end

function unionLowUppPWL(xLow::Vector{Float64}, zLow::Vector{Float64},
                        xUpp::Vector{Float64}, zUpp::Vector{Float64})
    xUnion = union(xLow, xUpp)
    sort!(xUnion)
    zLowUnion = Vector{Float64}(undef, length(xUnion))
    i, j = 1, 1
    while (i <= length(xUnion))
        if xLow[j] == xUnion[i]
            zLowUnion[i] = zLow[j]
            i += 1
            j += 1
        elseif xLow[j] < xUnion[i]
            i += 1
        else
            zLowUnion[i] = (xUnion[i] - xLow[j-1]) * (zLow[j] - zLow[j-1]) / (xLow[j] - xLow[j-1]) + zLow[j-1]
            i += 1
        end
    end
    zUppUnion = Vector{Float64}(undef, length(xUnion))
    i, j = 1, 1
    while (i <= length(xUnion))
        if xUpp[j] == xUnion[i]
            zUppUnion[i] = zUpp[j]
            i += 1
            j += 1
        elseif xUpp[j] < xUnion[i]
            i += 1
        else
            zUppUnion[i] = (xUnion[i] - xUpp[j-1]) * (zUpp[j] - zUpp[j-1]) / (xUpp[j] - xUpp[j-1]) + zUpp[j-1]
            i += 1
        end
    end
    return xUnion, zLowUnion, zUppUnion
end

function getDisaggregatedPWR!(xList::Vector{Float64}, zList::Vector{Float64},
                            SS::Vector{Set{Int64}})
    sLen = length(SS)
    count = length(xList) + 1
    for i in 1:(sLen-1)
        S = intersect(SS[i], SS[i+1])
        setdiff!(SS[i+1], S)
        union!(SS[i+1], Set(count:(count+length(S)-1)))
        count += length(S)
        for ele in S
            append!(xList, xList[ele])
            append!(zList, zList[ele])
        end
    end
    return
end

function getPWR(xListOrg::Vector{Float64}, act::Function, dAct::Function;
                times::Int64=1, threshold::Float64=1e-8)
    xList = [x for x in xListOrg]
    zList = [act(x) for x in xList]
    xLen = length(xList)
    count = xLen + 1
    SS = Vector{Set{Int64}}([])
    for i in 1:(xLen - 1)
        xx, zz = getMultiIntersection(xList[i], xList[i+1], act, dAct, times=times)
        append!(xList, xx)
        append!(zList, zz)
        S = [i, i+1]
        append!(S, count:(count+length(xx)-1))
        append!(SS, [Set(S)])
        count += length(xx)
        # if abs(dAct(xList[i]) - dAct(xList[i+1])) > threshold
        #     xx, zz = getIntersection(xList[i], zList[i], dAct(xList[i]),
        #                 xList[i+1], zList[i+1], dAct(xList[i+1]))
        #
        #     append!(SS, [Set([i, i+1, count])])
        #     count += 1
        # else
        #     append!(SS, [Set([i, i+1])])
        # end
    end
    return xList, zList, SS
end

function piecewiseLinearRelaxation!(m::JuMP.Model, x::VarOrAff, z::VarOrAff,
                            xList, zList, SS::Vector{Set{Int64}}; method=:BRGC)
    J = Set{Int64}()
    for S in SS
        J = union(J, S)
    end
    bc = getBicliqueCoverLinearTree(SS::Vector{Set{Int64}}, method=method)
    @assert(length(xList) == length(zList))
    lambda = @variable(m, [1:length(xList)], lower_bound=0.0)
    y = @variable(m, [1:length(bc)], binary=true)
    @constraint(m, sum(lambda) == 1)
    @constraint(m, sum(lambda[i] * xList[i] for i in 1:length(xList)) == x)
    @constraint(m, sum(lambda[i] * zList[i] for i in 1:length(zList)) == z)

    for i in 1:length(bc)
        @constraint(m, sum(lambda[u] for u in bc[i][1]) <= y[i])
        @constraint(m, sum(lambda[u] for u in bc[i][2]) <= 1 - y[i])

    end
    return
end

function getSeparation(SS::Vector{Set{Int64}})
    ssLen = length(SS)
    if (ssLen <= 1)
        return Vector{Tuple{Set{Int64}, Set{Int64}}}([])
    end
    cutIndex = floor(Int64, ssLen / 2)
    SS1 = SS[1:cutIndex]
    SS2 = SS[(cutIndex+1):ssLen]
    mid = intersect(SS[cutIndex], SS[cutIndex+1])
    L, R = Set{Int64}(), Set{Int64}()
    for S in SS1
        union!(L, S)
    end
    for S in SS2
        union!(R, S)
    end
    setdiff!(L, mid)
    setdiff!(R, mid)
    b = (L, R)
    bc1 = getSeparation(SS1)
    bc2 = getSeparation(SS2)

    insert!(bc1, 1, b)
    append!(bc1, bc2)

    return bc1
end

function isBiclique(g, bMerge)
    for u in bMerge[1]
        for v in bMerge[2]
            if (u == v) || !has_edge(g, u, v)
                return false
            end
        end
    end
    return true
end

function checkAndMerge(g, b1, b2)
    b2diff = (setdiff(b2[1], b1[1]), setdiff(b2[2], b1[2]))
    union!(b1[1], b2diff[1])
    union!(b1[2], b2diff[2])
    if isBiclique(g, b1)
        return true
    else
        setdiff!(b1[1], b2diff[1])
        setdiff!(b1[2], b2diff[2])
    end

    b2diff = (setdiff(b2[1], b1[2]), setdiff(b2[2], b1[1]))
    union!(b1[1], b2diff[2])
    union!(b1[2], b2diff[1])
    if isBiclique(g, b1)
        return true
    else
        setdiff!(b1[1], b2diff[2])
        setdiff!(b1[2], b2diff[1])
    end
    return false
end

function buildConflictGraph(SS::Vector{Set{Int64}})
    g = SimpleGraph()
    J = Set{Int64}()
    for S in SS
        J = union(J, S)
    end
    add_vertices!(g, length(J))
    for S in SS
        sLen = length(S)
        for u in S
            for v in S
                if (u < v)
                    add_edge!(g, u, v)
                end
            end
        end
    end
    return complement(g)
end

function isBicliqueCover(g,bc)
    g2 = SimpleGraph()
    add_vertices!(g2, maximum(vertices(g)))
    for b in bc
        for u in b[1]
            for v in b[2]
                if (u == v) || !has_edge(g, u, v) ||
                        !has_vertex(g, u) || !has_vertex(g, v)
                    return false
                end
                add_edge!(g2, u, v)
            end
        end
    end
    return collect(edges(g)) == collect(edges(g2))
end

function extendMaximalBiclique(g, b)
    flag = true
    vToAdd = setdiff(setdiff(Set(vertices(g)), b[1]), b[2])
    for v in vToAdd
        union!(b[1], Set([v]))
        if isBiclique(g, b)
            continue
        else
            setdiff!(b[1], Set([v]))
        end
        union!(b[2], Set([v]))
        if isBiclique(g, b)
            continue
        else
            setdiff!(b[2], Set([v]))
        end
    end
    return b
end

# Method takes :general, :BRGC, :balancedCode, :biclique
function getBicliqueCoverLinearTree(SS::Vector{Set{Int64}}; method=:biclique, maximal=false, check=false)
    res = Vector{Tuple{Set{Int64}, Set{Int64}}}([])
    if (length(SS) <= 1)
        return res
    end
    if (method == :general)
        g = buildConflictGraph(SS)
        bc = getSeparation(SS)
        bcLen = length(bc)
        append!(res, [bc[1]])
        println("After separation!")
        for i in 2:bcLen
            flag = false
            for j in 1:length(res)
                isMerge = checkAndMerge(g, res[j], bc[i])
                if (isMerge)
                    flag = true
                    break
                end
            end
            if (!flag)
                append!(res, [bc[i]])
            end
        end
    elseif (method == :BRGC) || (method == :balancedCode) || (method == :biclique)
        r = ceil(Int64, log2(length(SS)))
        if (method == :BRGC)
            code = reflectedBinaryCode(r)
        elseif (method == :balancedCode)
            code = bicliqueLinearCode(length(SS))
        elseif (method == :biclique)
            code = bicliqueAlgCode(length(SS))
        else
            error("Unrecognized method $method")
        end

        for i in 1:r
            L, R = Set{Int64}(), Set{Int64}()
            for j in 1:length(SS)
                if code[j, i] == 0
                    union!(L, SS[j])
                elseif code[j, i] == 1
                    union!(R, SS[j])
                end
            end
            mid = intersect(L, R)
            setdiff!(L, mid)
            setdiff!(R, mid)
            append!(res, [(L, R)])
        end
    else
        error("Unrecognized method $method")
    end
    if (maximal)
        for i in length(res)
            res[i] = extendMaximalBiclique(g, res[i])
        end
    end

    if (check) # Check whether it is biclique cover the conflict graph.
        g = buildConflictGraph(SS)
        @assert isBicliqueCover(g, res)
        println("Is it biclique cover ", isBicliqueCover(g, res))
    end
    return res
end

function linearCodeHelper!(change::Vector{Vector{Int64}}, left::Int64,
                            right::Int64, level::Int64)
    if (left == right)
        return
    end
    mid = ceil(Int64, (left+right) / 2)
    append!(change[level], mid)
    linearCodeHelper!(change, left, mid-1, level+1)
    linearCodeHelper!(change, mid, right, level+1)
end

function bicliqueLinearCode(d::Int64)
    r = ceil(Int64, log2(d))
    code = zeros(Int64, d, r)
    change = Vector{Vector{Int64}}([[] for i in 1:r])
    linearCodeHelper!(change, 1, d, 1)
    for i in 1:r
        val = 0
        for j in 1:d
            if (j in change[i])
                val = 1 - val
            end
            code[j, i] = val
        end
    end
    return code
end

function algCodeHelper!(change::Vector{Vector{Int64}},
                    leave::Vector{Vector{Int64}},
                    left::Int64,right::Int64, level::Int64)
    if (left == right)
        append!(leave[level], left)
        return
    end
    mid = ceil(Int64, (left+right) / 2)
    append!(change[level], mid)
    algCodeHelper!(change, leave, left, mid-1, level+1)
    algCodeHelper!(change, leave, mid, right, level+1)
end

function bicliqueAlgCode(d::Int64)
    r = ceil(Int64, log2(d))
    code = zeros(Int64, d, r)
    change = Vector{Vector{Int64}}([[] for i in 1:r])
    leave = Vector{Vector{Int64}}([[] for i in 1:(r+1)])
    algCodeHelper!(change, leave, 1, d, 1)
    for i in 1:r
        val = 0
        for j in 1:d
            if (j in change[i])
                val = 1 - val
            end
            code[j, i] = val
        end
        for j in leave[i]
            code[j, i] = -1
        end
    end
    return code
end

function reflectedBinaryCode(r::Int64)
    @assert (r > 0)
    if (r == 1)
        return [0 ; 1;;]
    else
        code = reflectedBinaryCode(r-1)
        codeRev = reverse(code, dims=1)
        vec = vcat([0 for i in 1:2^(r-1)], [1 for i in 1:2^(r-1)])
        return hcat(vec, vcat(code, codeRev))
    end
end


# Method takes strings "code" or "biclique".
function getSOSkBicliqueCover(b::Int64, k::Int64)
    SS = [Set([j for j in i:(i+k-1)]) for i in 1:(2^b)]
    # println(SS)
    d = Dict{Tuple{Int64, Int64}, Tuple{Set{Int64}, Set{Int64}}}()
    if (length(SS) <= 1)
        return res
    end
    g = buildConflictGraph(SS)

    for i in 0:(b-1)
        for j in 0:(2^i-1)
            # println("i ", i, "j ", j)
            d[(i,j)] = (Set([l for l in (1 + j * 2^(b-i)):((2 * j + 1)*2^(b-i-1))]),
                        Set([l for l in ((2*j+1) * 2^(b-i-1) +k):((j+1) * 2^(b-i) + k-1)]))
            # println(d[(i,j)])
        end
    end
    res = [d[l] for l in keys(d)]
    # println(res)
    println("Before Merge Is it biclique cover ", isBicliqueCover(g, res))
    @assert isBicliqueCover(g, res)


    res = Vector{Tuple{Set{Int64}, Set{Int64}}}([])
    for i in 0:(b-1)
        alphai = ceil(Int64, (k-1 + 2^(b-i-1)) / (2^(b-i)) )
        # println(alphai)

        for p in 0:(min(alphai-1, 2^i - 1))
            L, R = Set{Int64}([]), Set{Int64}([])
            for q in 0:(floor(Int64, (2^i - 1 - p) / (2 * alphai) ))
                union!(L, d[(i, 2*q*alphai+p)][1])
                union!(R, d[(i, 2*q*alphai+p)][2])
            end

            for q in 0:(floor(Int64, (2^i - 1 - p) / (2 * alphai) - (1/2)))
                union!(L, d[(i, (2*q+1)*alphai+p)][2])
                union!(R, d[(i, (2*q+1)*alphai+p)][1])
            end
            append!(res, [(L, R)])
        end
    end
    for ele in res
        # println(isBiclique(g, ele))
        # println(ele)
        @assert isBiclique(g, ele)
    end
    # println(res)
    println("After Merge Is it biclique cover ", isBicliqueCover(g, res))
    @assert isBicliqueCover(g, res)

    println("The value of b: ", b, " k: ", k)
    println("The biclique cover size: ", length(res))
    @assert (length(res) <= b + k - 2)

    return res
end
