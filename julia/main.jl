macro run()
    :(include("main.jl") ; main())
end

using DataFrames
using CSV
using Plots
using Printf
using DataInterpolations
import TOV

function main()
    eoss = Vector{TOV.EoS}()
    neos = 2000
    for i in 1:neos
        push!(eoss, TOV.EoS("../eos/out/eos"*string(i)*".csv", ["P", "ϵ"], :linear_interpolation))
    end

    P0i = 5.0
    P0f = 800.0
    nstars = 100
    u = range(log(P0i), log(P0f), length=nstars)
    P0 = exp.(u).*TOV.MeVfm3

    @time Threads.@threads for i in 1:neos
        eos = eoss[i]
        mrdiagram = TOV.solvemrdiagram(P0, P->eos(P), 1.0TOV.SI_TO_LENGTH_UNIT)
        writedat("out/mrdiagram"*string(i)*".dat", mrdiagram[:,1], mrdiagram[:,2], mrdiagram[:,3])
    end
end

function formatandjoin(v::AbstractVector)::String
    vstrvec = []
    for vi in v
        push!(vstrvec, @sprintf("%.16e", vi))
    end

    return join(vstrvec, ' ')
end

# I can implement a variable space version
function writedat(file::String, columns::AbstractVector...)
    io = open(file, "w")
    try
        for line = zip(columns...)
            linevec = collect(line)
            linestring = formatandjoin(linevec)
            write(io, linestring*'\n')
        end
    finally
        close(io)
    end
end
