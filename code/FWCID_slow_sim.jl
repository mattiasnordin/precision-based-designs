using Distributions
using Random
using Plots
using DataFrames
using CSV
using Dates

function FWCID_bernoulli_trial_slow(distrY0, distrY1, dlist, α, αc, p, shift)
    nbarlist = (4 * quantile(Normal(), 1 - α / 2)^2) ./ dlist.^2
    a = 1 / sqrt((1-p) * var(distrY0) + p * var(distrY1))
    nr_d = length(dlist)
    z_α = quantile(Normal(), 1-α/2)
    ρ_clist = bin_search.(100, .0000001, .00000001, nbarlist, αc, true)
    ρ_avlist = bin_search.(100, .0000001, .00000001, nbarlist, α, false)
    τ = (mean(distrY1) - mean(distrY0)) * a + shift

    N = Array{Int64}(undef, nr_d, 3)
    ATE_EST = Array{Float64}(undef, nr_d, 3)
    COVERAGE = Array{Bool}(undef, nr_d, 3)
    NEXT = [1, 1, 1]

    W, Y, n0, n1, n = initialize_data(distrY0, distrY1, p, shift, a)

    while minimum(NEXT) <= nr_d
        n = n0 + n1
        κ_n = n / (n0 * n1)
        μhat_1n = mean(Y[W .== 1])
        μhat_0n = mean(Y[W .== 0])
        σ2hat_1n = mean((Y[W .== 1] .- μhat_1n).^2)
        σ2hat_0n = mean((Y[W .== 0] .- μhat_0n).^2)
        σ2hat_pn = n0/n * σ2hat_1n + n1/n * σ2hat_0n 
        σ2_τn = σ2hat_pn * κ_n
        τhat_n = μhat_1n - μhat_0n
        phat_n = n1 / n


        for i = NEXT[1]:nr_d
            if σ2_τn ≤ dlist[i]^2 / z_α^2
                NEXT[1] += 1
                N[i, 1] = n
                ATE_EST[i, 1] = τhat_n
                COVERAGE[i, 1] = (τhat_n - dlist[i] .<= τ .<=
                                  τhat_n + dlist[i])
            end
        end

        for i = NEXT[2]:nr_d
            Z = W .* (Y .- μhat_1n).^2 .+ (1 .- W) .* (Y .- μhat_0n).^2
            vZ = mean((Z .- σ2hat_pn).^2)
            ρ_c = ρ_clist[i]
            ϕ_n = get_phi(ρ_c, αc, n)
            σ2_p_ub_n = σ2hat_pn + sqrt(vZ) * ϕ_n
            if σ2_p_ub_n * κ_n ≤ dlist[i]^2 / z_α^2
                NEXT[2] += 1
                N[i, 2] = n
                ATE_EST[i, 2] = τhat_n
                COVERAGE[i, 2] = (τhat_n - dlist[i] .<= τ .<=
                                  τhat_n + dlist[i])
            end
        end

        for i = NEXT[3]:nr_d
            v2_n = 1 / phat_n * (σ2hat_1n + μhat_1n^2) +
                   1 / (1 - phat_n) * (σ2hat_0n + μhat_0n^2) - 
                   (μhat_1n - μhat_0n)^2
            ρ_av = ρ_avlist[i]
            ψ_n = get_psi(ρ_av, α, n)
            if v2_n ≤ dlist[i]^2 / ψ_n^2
                NEXT[3] += 1
                N[i, 3] = n
                ATE_EST[i, 3] = τhat_n
                COVERAGE[i, 3] = (τhat_n - dlist[i] .<= τ .<=
                                  τhat_n + dlist[i])
            end
        end

        if minimum(NEXT) <= nr_d
            w = rand(Bernoulli(p))
            if w == true
                y = rand(distrY1) * a + shift
                n1 += 1
            elseif w == false
                y = rand(distrY0) * a
                n0 += 1
            end
            push!(W, w)
            push!(Y, y)
        end
    end
    return hcat(N, ATE_EST, COVERAGE, dlist)
end

include("auxiliary_functions.jl")

α = .1
αc = .1
p = .5

distrlist = simulation_distributions()

OUTFILE = "simulations\\out_FWCID_slow.csv"
OUTFILE_AGG = "simulations\\out_FWCID_slow_agg.csv"

REPS = parse(Int64, ARGS[1])
dMIN = parse(Float64, ARGS[2])
dMAX = parse(Float64, ARGS[3])
dSTEP = parse(Float64, ARGS[4])

dlist = dMAX:-dSTEP:dMIN
Random.seed!(12345)

START_TIME = time_ns()
for i = 1:REPS
    for j = 1:length(distrlist)
        DY0 = distrlist[j][1]
        DY1 = distrlist[j][2]
        sh = distrlist[j][3]
        v = FWCID_bernoulli_trial_slow(DY0, DY1, dlist, α, αc, p, sh)
        df = DataFrame(v, :auto)
        rename!(df, [:n_naive, :n_cons, :n_av,
                     :te_naive, :te_cons, :te_av,
                     :cov_naive, :cov_cons, :cov_av, :d])
        df[!, :distrYid] .= j
        if isfile(OUTFILE)
            CSV.write(OUTFILE, df, append=true)
        else
            CSV.write(OUTFILE, df)
        end
    end
    print_time(i, REPS, START_TIME)
end

df = CSV.read(OUTFILE, DataFrame)
group_vars = [:d, :distrYid]
gdf = groupby(df, group_vars)
df_agg = combine(gdf, valuecols(gdf) .=> mean)
df_agg[:, :nbar] = (4 * quantile(Normal(), 1 - α / 2)^2) ./ df_agg[:, :d].^2

CSV.write(OUTFILE_AGG, df_agg)
