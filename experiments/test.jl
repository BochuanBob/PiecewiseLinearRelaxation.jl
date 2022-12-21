# cli.jl

println(map(x->string(x, x), ARGS))

for arg in ARGS
    println(arg, typeof(arg))
end
