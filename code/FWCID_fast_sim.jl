using Distributions
using Random
using Plots
using DataFrames
using CSV
using Dates

function FWCID_bernoulli_trial_fast(distrY0, distrY1, dlist, α, αc, p, shift=0)
    expnlist = (4 * quantile(Normal(), 1 - α / 2)^2) ./ dlist.^2
    nr_d = length(dlist)
    a = 1 / sqrt((1-p) * var(distrY0) + p * var(distrY1))
    n0, n1, n = 0, 0, 0
    sY11, sY12, sY13, sY14 = 0, 0, 0, 0
    sY01, sY02, sY03, sY04 = 0, 0, 0, 0
    rholist2 = bin_search.(100, .000001, .0000001, expnlist, αc, true)
    rholist3 = bin_search.(100, .000001, .0000001, expnlist, α, false)
    τ = (mean(distrY1) - mean(distrY0)) * a .+ shift
    goal_tau_varlist = dlist.^2 ./ (quantile(Normal(), 1-α/2)^2)

    STARTj1 = 0
    STARTj2 = 0
    STARTj3 = 0
    STOPn1 = zeros(nr_d)
    STOPn2 = zeros(nr_d)
    STOPn3 = zeros(nr_d)
    TAU_EST1 = Array{Float64}(undef, nr_d)
    TAU_EST2 = Array{Float64}(undef, nr_d)
    TAU_EST3 = Array{Float64}(undef, nr_d)
    COVERAGE1 = Array{Bool}(undef, nr_d)
    COVERAGE2 = Array{Bool}(undef, nr_d)
    COVERAGE3 = Array{Bool}(undef, nr_d)
    stop1, stop2, stop3 = false, false, false
                   
    mY0, mY1 = 0, 0
    while (n0 < 2) | (n1 < 2)
        w = rand(Bernoulli(p))
        if w == true
            y = rand(distrY1) * a + shift
            n1 += 1
            n += 1
            sY11 += y
            sY12 += y^2
            sY13 += y^3
            sY14 += y^4
            mY1 = sY11 / n1
        elseif w == false
            y = rand(distrY0) * a
            n0 += 1
            n += 1
            sY01 += y
            sY02 += y^2
            sY03 += y^3
            sY04 += y^4
            mY0 = sY01 / n0
        end
    end
    D = Bernoulli(p)
    while (stop1  == false) | (stop2 == false) | (stop3 == false)
        phat = n1 / n
        κ = n / (n0 * n1)
        sigma21 = sY12 / n1 - mY1^2
        sigma20 = sY02 / n0 - mY0^2
        sigma2p = (1 - phat) * sigma21 + phat * sigma20
        sigma2tau = sigma2p * κ
        j = STARTj1
        contd = true
        if stop1 == false
            while contd & (j < nr_d)
                j += 1
                if STOPn1[j] == 0
                    if sigma2tau <= goal_tau_varlist[j]
                        STOPn1[j] = n
                        tau_est = mY1 - mY0
                        TAU_EST1[j] = tau_est
                        COVERAGE1[j] = (tau_est - dlist[j] .<= τ .<=
                                        tau_est + dlist[j])
                        STARTj1 += 1
                        if j == nr_d
                            stop1 = true
                        end
                    else
                        contd = false
                    end
                end
            end
        end
        j = STARTj2
        contd = true
        if stop2 == false
            while contd & (j < nr_d)
                j += 1
                rho = rholist2[j]
                if STOPn2[j] == 0
                    a1 = sY14 + 6 * mY1^2 * sY12 - 4 * mY1 * sY13 - 3 * n1 * mY1^4 + n1 * sigma2p^2 - 2*sigma2p * (sY12 - n1 * mY1^2)
                    a2 = sY04 + 6 * mY0^2 * sY02 - 4 * mY0 * sY03 - 3 * n0 * mY0^4 + n0 * sigma2p^2 - 2*sigma2p * (sY02 - n0 * mY0^2)
                    var_σ2p_hat = (a1 + a2) / n
                    if var_σ2p_hat < 0
                        println("Error ", "n")
                        var_σ2p_hat = 100000
                    end
                    av_se = sqrt(var_σ2p_hat)*sqrt((2*(n*rho^2+1))/(n^2*rho^2) *
                                    log(1+(sqrt(n*rho^2+1))/(2*αc)))
                    if (sigma2p + av_se) * κ <= goal_tau_varlist[j]
                        STOPn2[j] = n
                        tau_est = mY1 - mY0
                        TAU_EST2[j] = tau_est
                        COVERAGE2[j] = (tau_est - dlist[j] .<= τ .<=
                                        tau_est + dlist[j])
                        STARTj2 += 1
                        if j == nr_d
                            stop2 = true
                        end
                    else
                        contd = false
                    end
                end
            end
        end
        j = STARTj3
        contd = true
        if stop3 == false
            while contd & (j < nr_d)
                j += 1
                rho = rholist3[j]
                if STOPn3[j] == 0
                    vR = (1 / (n1 * phat)) * sY12 - (2 - phat) * mY1^2 + phat * mY0^2 + (1 - phat) * 2 * mY1 * mY0 +
                         1 / (n * (1 - phat)^2) * sY02 + (1-phat) * mY1^2 - (1 + phat) * mY0^2 + 2 * phat * mY0 * mY1
                    stop_crit1 = sqrt(vR)*sqrt((2*(n*rho^2+1))/(n^2*rho^2) *
                                 log((sqrt(n*rho^2+1))/α))
                    if stop_crit1 <= dlist[j]
                        STOPn3[j] = n
                        tau_est = mY1 - mY0
                        TAU_EST3[j] = tau_est
                        COVERAGE3[j] = (tau_est - dlist[j] .<= τ .<=
                                        tau_est + dlist[j])
                        STARTj3 += 1
                        if j == nr_d
                            stop3 = true
                        end
                    else
                        contd = false
                    end
                end
            end
        end
        if (stop1 == false) | (stop2 == false) | (stop3 == false)
            w = rand(D)
            if w == true
                y = rand(distrY1) * a + shift
                n1 += 1
                n += 1
                sY11 += y
                sY12 += y^2
                sY13 += y^3
                sY14 += y^4
                mY1 = sY11 / n1
            elseif w == false
                y = rand(distrY0) * a
                n0 += 1
                n += 1
                sY01 += y
                sY02 += y^2
                sY03 += y^3
                sY04 += y^4
                mY0 = sY01 / n0
            end
        end
    end
    return hcat(STOPn1, STOPn2, STOPn3, TAU_EST1, TAU_EST2, TAU_EST3, COVERAGE1, COVERAGE2, COVERAGE3, dlist)
end          

include("auxiliary_functions.jl")

α = .1
αc = .1
p = .5

distrlist = simulation_distributions()

OUTFILE = "simulations\\out_FWCID_fast.csv"
OUTFILE_AGG = "simulations\\out_FWCID_fast_agg.csv"

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
        v = FWCID_bernoulli_trial_fast(DY0, DY1, dlist, α, αc, p, distrlist[j][3])
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
    if i % 10 == 0
        print_time(i, REPS, START_TIME)
    end
end

df = CSV.read(OUTFILE, DataFrame)
group_vars = [:d, :distrYid]
gdf = groupby(df, group_vars)
df_agg = combine(gdf, valuecols(gdf) .=> mean)
df_agg[:, :nbar] = (4 * quantile(Normal(), 1 - α / 2)^2) ./ df_agg[:, :d].^2

CSV.write(OUTFILE_AGG, df_agg)
q = CSV.write(OUTFILE_AGG, df_agg)
