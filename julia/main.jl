macro run()
    :(include("main.jl") ; main())
end

using DataFrames
using CSV
using Plots
using Printf
import TOV

function main()
    eos = TOV.EoS("../eos_2.csv")
    P0 = [range(0.1TOV.MeVfm3, 10TOV.MeVfm3, length=50) ; range(10TOV.MeVfm3, 600TOV.MeVfm3, length=200)]
    @time "Solving MR diagram" mrdiagram = TOV.solvemrdiagram(P0, P->eos(P), 1.0TOV.SI_TO_LENGTH_UNIT)
    writedat("out/mrdiagram.dat", mrdiagram[:,1], mrdiagram[:,2], mrdiagram[:,3])
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
