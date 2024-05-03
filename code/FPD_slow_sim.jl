using Distributions
using Random
using Plots
using DataFrames
using CSV
using Dates

function FPD_bernoulli_trial_slow(distrY0, distrY1, τ_H0, τ_H1, α, αc,
                                   β, p, shift, zlist, maxn)
    τ_d = τ_H1 - τ_H0
    z_α = quantile(Normal(), 1 - α)
    z_β = quantile(Normal(), 1 - β)
    nbar = (z_α + z_β)^2 / (τ_d^2 * (p * (1-p)))
    a = 1 / sqrt((1-p) * var(distrY0) + p * var(distrY1))
    ρ_cons = bin_search(100, .0000001, .00000001, nbar, αc, true)
    ρ_av = bin_search(100, .0000001, .00000001, nbar, α, true)

    N = Array{Int64}(undef, 4)
    ATE_EST = Array{Float64}(undef, 4)
    REJECT = Array{Bool}(undef, 4)
    STOP = [false, false, false, false]
    
    W, Y, n0, n1, n = initialize_data(distrY0, distrY1, p, shift, a)

    while (minimum(STOP) == false) & (n <= maxn)
        κ_n = n / (n0 * n1)
        μ1hat_n = mean(Y[W .== 1])
        μ0hat_n = mean(Y[W .== 0])
        σ2hat_1n = mean((Y[W .== 1] .- μ1hat_n).^2)
        σ2hat_0n = mean((Y[W .== 0] .- μ0hat_n).^2)
        σ2hat_pn = n0/n * σ2hat_1n + n1/n * σ2hat_0n 
        σ2_τn = σ2hat_pn * κ_n
        τ_Hat_n = μ1hat_n - μ0hat_n
        phat_n = n1 / n

        if STOP[1] == false
            if σ2_τn ≤ τ_d^2 / (z_α + z_β)^2
                STOP[1] = true
                N[1] = n
                ATE_EST[1] = τ_Hat_n
                REJECT[1] = τ_Hat_n - τ_d * z_α / (z_α + z_β) > τ_H0
            end
        end

        if STOP[2] == false
            Z = W .* (Y .- μ1hat_n).^2 .+ (1 .- W) .* (Y .- μ0hat_n).^2
            vZ = mean((Z .- σ2hat_pn).^2)
            ϕ_n = get_phi(ρ_cons, αc, n)
            σ2_p_ub_n = σ2hat_pn + sqrt(vZ) * ϕ_n
            if σ2_p_ub_n * κ_n ≤ τ_d^2 / (z_α + z_β)^2
                STOP[2] = true
                N[2] = n
                ATE_EST[2] = τ_Hat_n
                REJECT[2] = τ_Hat_n - τ_d * z_α / (z_α + z_β) > τ_H0
            end
        end

        if STOP[3] == false
            v_n = sqrt(1 / phat_n * (σ2hat_1n + μ1hat_n^2) + 
                       1 / (1-phat_n) * (σ2hat_0n + μ0hat_n^2) -
                       (μ1hat_n - μ0hat_n)^2)
            ϕ_n = get_phi(ρ_av, α, n)
            if -(τ_Hat_n -τ_H0) / v_n ≤ -ϕ_n
                STOP[3] = true
                N[3] = n
                ATE_EST[3] = τ_Hat_n
                REJECT[3] = true
            end
        end

        if STOP[4] == false
            if -τ_Hat_n / sqrt(σ2_τn) ≤ -zlist[n]
                STOP[4] = true
                N[4] = n
                ATE_EST[4] = τ_Hat_n
                REJECT[4] = true
            end
        end

        if minimum(STOP) == false
            if n < maxn
                w = rand(Bernoulli(p))
                if w == true
                    y = rand(distrY1) * a + shift
                    n1 += 1
                elseif w == false
                    y = rand(distrY0) * a
                    n0 += 1
                end
                n = n0 + n1
                push!(W, w)
                push!(Y, y)
            elseif n == maxn
                n += 1
            end
        end
    end
    for i = 1:4
        if STOP[i] == false
            μ1hat_n = mean(Y[W .== 1])
            μ0hat_n = mean(Y[W .== 0])
            τ_Hat_n = μ1hat_n - μ0hat_n
            N[i] = n - 1
            ATE_EST[i] = τ_Hat_n
            REJECT[i] = false
        end
    end
    return hcat(N, ATE_EST, REJECT)
end

include("auxiliary_functions.jl")

α = .05
αc = .1
p = .5
β = .2
z_α = quantile(Normal(), 1 - α)
z_β = quantile(Normal(), 1 - β)
τ_H0 = 0

distrlist = simulation_distributions()

OUTFILE = "simulations\\out_FPD_slow.csv"
ZFILE = "simulations\\gst_z.csv"

REPS = parse(Int64, ARGS[1])
τ_H1MIN = parse(Float64, ARGS[2])
τ_H1MAX = parse(Float64, ARGS[3])
τ_H1STEP = parse(Float64, ARGS[4])

dfz = CSV.read(ZFILE, DataFrame)
Random.seed!(12345)
for τ_H1 = τ_H1MAX:-τ_H1STEP:τ_H1MIN
    for mult_nbar = [1, 1.5, 2]
        nbar = (z_α + z_β)^2 / ((τ_H1 - τ_H0)^2 * (p * (1-p)))
        maxn = Int(round(nbar * mult_nbar))
        zlist = dfz[(dfz[!, :tau_H1] .== τ_H1) .&
                    (dfz[!, :mult_maxn] .== mult_nbar), :zval]
        OUT2 = Array{Float64}(undef, length(distrlist), 12)
        for j = 1:length(distrlist)
            OUT = Array{Float64}(undef, REPS, 12)
            START_TIME = time_ns()
            for i = 1:REPS
                q = FPD_bernoulli_trial_slow(distrlist[j][1], distrlist[j][2], τ_H0,
                    τ_H1, α, αc, β, p, distrlist[j][3], zlist, maxn)
                OUT[i, :] = q
            end
            df = DataFrame(mean(OUT, dims=1), :auto)
            rename!(df, [:n1, :n2, :n3, :n4, :est1, :est2, :est3, :est4,
                        :power1, :power2, :power3, :power4])
            df[!, :distrYid] .= j
            df[!, :maxn] .= maxn
            df[!, :tau_H1] .= τ_H1
            if isfile(OUTFILE)
                CSV.write(OUTFILE, df, append=true)
            else
                CSV.write(OUTFILE, df)
            end
            println(τ_H1, " ", mult_nbar, " ", j)
        end
    end
end
