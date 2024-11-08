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
    # eoss = Vector{TOV.EoS}()
    # neos = 100
    # for i in 1:neos
    #     push!(eoss, TOV.EoS("../eos/out/eos"*string(i)*".csv", ["P", "ϵ"], :linear_interpolation))
    # end

    # P0i = 5.0
    # P0f = 800.0
    # nstars = 100
    # u = range(log(P0i), log(P0f), length=nstars)
    # P0 = exp.(u).*TOV.MeVfm3

    # @time Threads.@threads for i in 1:neos
    #     eos = eoss[i]
    #     mrdiagram = TOV.solvemrdiagram(P0, P->eos(P), 1.0TOV.SI_TO_LENGTH_UNIT)
    #     writedat("out/mrdiagram"*string(i)*".dat", mrdiagram[:,1], mrdiagram[:,2], mrdiagram[:,3])
    # end
    
    eos1 = TOV.EoS("../EOSTRS1.3.csv", ["P", "ϵ"], :linear_interpolation)
    P0i = 5.0
    P0f = 480.0
    nstars = 100
    u = range(log(P0i), log(P0f), length=nstars)
    P0 = exp.(u).*TOV.MeVfm3

    push!(eos1.itp.a, last(eos1.itp.a))
    push!(eos1.itp.b, last(eos1.itp.b))

    mrdiagram = TOV.solvemrdiagram(P0, P->eos1(P), 10.0TOV.SI_TO_LENGTH_UNIT)
    writedat("out/0mrdiagram.dat", mrdiagram[:,1], mrdiagram[:,2], mrdiagram[:,3])
    
    eos2 = TOV.EoS("../EOSMSS1.3.csv", ["P", "ϵ"], :linear_interpolation)

    mrdiagram = TOV.solvemrdiagram(P0, P->eos2(P), 1.0TOV.SI_TO_LENGTH_UNIT)
    writedat("out/1mrdiagram.dat", mrdiagram[:,1], mrdiagram[:,2], mrdiagram[:,3])

    plot(P0, eos1.(P0), label="EOS TRS1.3", xaxis="Pressure [MeV/fm³]", yaxis="Energy density [MeV/fm³]", title="EOS TRS1.3")
    plot!(P0, eos2.(P0), label="EOS MSS1.3")
    gui()
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
