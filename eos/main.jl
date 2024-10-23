macro run()
    :(include("main.jl") ; main())
end

using DataFrames
using CSV
using Plots
using Printf

function getK(ρ, P, Γ)
    return P / ρ^Γ
end

function P(ρ, K, Λ)
    return K * ρ^Λ
end

function ϵ(ρint, Pint, ϵint, ρ, Γ, K)
    if Γ == 1
        return ϵint/ρint * ρ + K*log(1/ρint)*ρ - K*log(1/ρ)*ρ
    end
    
    return (ϵint/ρint - Pint/(ρint*(Γ - 1)))*ρ + K/(Γ - 1) * ρ^Γ
end

struct Polytrope
    Γvals::Vector{Float64}
    Kvals::Vector{Float64}
    ρintvals::Vector{Float64}
    ϵintvals::Vector{Float64}
    Pintvals::Vector{Float64}
    maxρ::Float64
end

function (pol::Polytrope)(ρ)
    if ρ in pol.ρintvals
        ind = findfirst(x -> x == ρ, pol.ρintvals)
        return pol.Pintvals[ind], pol.ϵintvals[ind]
    end

    for i in 1:length(pol.ρintvals)
        if i != length(pol.ρintvals) && pol.ρintvals[i] < ρ < pol.ρintvals[i+1]
            return P(ρ, pol.Kvals[i], pol.Γvals[i]), ϵ(pol.ρintvals[i], pol.Pintvals[i], pol.ϵintvals[i], ρ, pol.Γvals[i], pol.Kvals[i])
        end
    end

    return P(ρ, last(pol.Kvals), last(pol.Γvals)), ϵ(pol.ρintvals[end], pol.Pintvals[end], pol.ϵintvals[end], ρ, last(pol.Γvals), last(pol.Kvals))
end

function main()
    ρ0 = 0.16 # fm^-3
    ρvals = [1.0, 1.4, 2.2, 3.3, 4.9].*ρ0
    maxρ = 7.4*ρ0
   
    crustdf = CSV.File("crust_sly.csv") |> DataFrame
    crustρ = crustdf[:,1]
    crustP = crustdf[:,2]
    crustϵ = crustdf[:,3]
    crustcs2 = crustdf[:,4]

    ρc = last(crustρ)
    cs2c = last(crustcs2)
    ϵc = last(crustϵ)
    Pc = last(crustP)
    Γc = cs2c * (ϵc + Pc) / Pc
    Kc = Pc / ρc^Γc

    transform(col, val) = val
    transform(col, val::Float64) = @sprintf("%.16e", val)

    neos = 100

    for i in range(1, neos)
        Γvals = 1 .+ 4 .*rand(Float64, length(ρvals))
        Kvals = zeros(length(Γvals))
        Kvals[1] = getK(ρc, Pc, Γvals[1])
        for i in 2:length(Γvals)
            ρint = (ρvals[i])
            Pint = P(ρint, Kvals[i-1], Γvals[i-1])
            Kvals[i] = getK(ρint, Pint, Γvals[i])
        end
        ϵintvals = zeros(length(Γvals))
        Pintvals = zeros(length(Γvals))

        ϵintvals[1] = ϵc
        Pintvals[1] = Pc
        
        for i in 2:length(ϵintvals)
            ρint = (ρvals[i])
            Pint = P(ρint, Kvals[i-1], Γvals[i-1])
            ϵint = ϵ(ρvals[i-1], Pintvals[i-1], ϵintvals[i-1], ρint, Γvals[i-1], Kvals[i-1])
            ϵintvals[i] = ϵint
            Pintvals[i] = Pint
        end

        poly = Polytrope(Γvals, Kvals, ρvals, ϵintvals, Pintvals, maxρ)
        ρvals = range(1.0ρ0, maxρ, length=300)
        Pvals = zeros(length(ρvals))
        ϵvals = zeros(length(ρvals))
        for i in 1:length(ρvals)
            Pval, ϵval = poly(ρvals[i])
            Pvals[i] = Pval
            ϵvals[i] = ϵval
        end

        df = DataFrame(P=[crustP;Pvals], ϵ=[crustϵ;ϵvals])
        CSV.write("out/eos$i.csv", df, writeheader=false, transform=transform)
    end

    plot()
    for i in range(1, neos)
        df = CSV.File("out/eos$i.csv", header=["P", "ϵ"]) |> DataFrame
        plot!(df.ϵ, df.P, label=false)
    end
    gui()
end