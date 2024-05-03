using Distributions
using Random
using Plots
using DataFrames
using CSV
using Dates

function FPD_bernoulli_trial_fast(distrY0, distrY1, τ_H0, τ_d, α, αc,
    z_β, p, goal_tau_var, expn, maxn, shift, zlist)
    a = 1 / sqrt((1-p) * var(distrY0) + p * var(distrY1))
    n0, n1, n = 0, 0, 0
    sY11, sY12, sY13, sY14 = 0, 0, 0, 0
    sY01, sY02, sY03, sY04 = 0, 0, 0, 0
    rho2 = bin_search(100, .0000001, .00000001, expn, αc, true)
    rho3 = bin_search(100, .0000001, .00000001, expn, α, true)
    STOP = [false, false, false, false]
    STOPn = Array{Int64}(undef, 4)
    TAU_EST = Array{Float64}(undef, 4)
    REJECT = Array{Bool}(undef, 4)
                   
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
    while (minimum(STOP) == false) & (n <= maxn)
        phat = n1 / n
        κ = n / (n0 * n1)
        sigma21 = sY12 / n1 - mY1^2
        sigma20 = sY02 / n0 - mY0^2
        sigma2p = (1 - phat) * sigma21 + phat * sigma20
        sigma2tau = sigma2p * κ
        tau_est = mY1 - mY0
        if STOP[1] == false
            if sigma2tau <= goal_tau_var
                STOP[1] = true
                STOPn[1] = n
                TAU_EST[1] = tau_est
                REJECT[1] = tau_est - τ_d * z_α / (z_α + z_β) > τ_H0
            end
        end
        if STOP[2] == false
            rho = rho2
            a1 = sY14 + 6 * mY1^2 * sY12 - 4 * mY1 * sY13 - 3 * n1 * mY1^4 + n1 * sigma2p^2 - 2*sigma2p * (sY12 - n1 * mY1^2)
            a2 = sY04 + 6 * mY0^2 * sY02 - 4 * mY0 * sY03 - 3 * n0 * mY0^4 + n0 * sigma2p^2 - 2*sigma2p * (sY02 - n0 * mY0^2)
            var_σ2p_hat = (a1 + a2) / n
            if var_σ2p_hat < 0
                println("Error ", "n")
                var_σ2p_hat = 100000
            end
            av_se = sqrt(var_σ2p_hat)*sqrt((2*(n*rho^2+1))/(n^2*rho^2) *
                                    log(1+(sqrt(n*rho^2+1))/(2*αc)))
            if (sigma2p + av_se) * κ <= goal_tau_var
                STOP[2] = true
                STOPn[2] = n
                TAU_EST[2] = tau_est
                REJECT[2] = tau_est - τ_d * z_α / (z_α + z_β) > τ_H0
            end
        end
        if STOP[3] == false
            rho = rho3
            vR = (1 / (n1 * phat)) * sY12 - (2 - phat) * mY1^2 + phat * mY0^2 + (1 - phat) * 2 * mY1 * mY0 +
                         1 / (n * (1 - phat)^2) * sY02 + (1-phat) * mY1^2 - (1 + phat) * mY0^2 + 2 * phat * mY0 * mY1
            stop_crit1 = sqrt(vR)*sqrt((2*(n*rho^2+1))/(n^2*rho^2) *
                                 log(1 + (sqrt(n*rho^2+1))/(2 * α)))
            if tau_est - stop_crit1 >= τ_H0
                STOP[3] = true
                STOPn[3] = n
                TAU_EST[3] = tau_est
                REJECT[3] = true
            end
        end
        if STOP[4] == false
            zstat = tau_est / sqrt(sigma2tau)
                    stop_crit1 = 0
            if zstat >= zlist[n]
                STOP[4] = true
                STOPn[4] = n
                TAU_EST[4] = tau_est
                REJECT[4] = true
            end
        end
        if minimum(STOP) == false
            if n < maxn
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
            elseif n == maxn
                n += 1
            end
        end
    end
    tau_est = mY1 - mY0
    for i = 1:4
        if STOP[i] == false
            STOPn[i] = n - 1
            TAU_EST[i] = tau_est
            REJECT[i] = false
        end
    end
    OUT = vcat(STOPn, TAU_EST, REJECT)
    return OUT 
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

OUTFILE = "simulations\\out_FPD_fast.csv"
ZFILE = "simulations\\gst_z.csv"

REPS = parse(Int64, ARGS[1])
τ_H1MIN = parse(Float64, ARGS[2])
τ_H1MAX = parse(Float64, ARGS[3])
τ_H1STEP = parse(Float64, ARGS[4])

dfz = CSV.read(ZFILE, DataFrame)
Random.seed!(12345)
for τ_H1 = τ_H1MAX:-τ_H1STEP:τ_H1MIN
    τ_d = τ_H1 - τ_H0
    goal_tau_var = τ_d^2 / (z_α + z_β)^2
    for mult_maxn = [1, 1.5, 2]
        expn = (z_α + z_β)^2 / (τ_d^2 * (p * (1-p)))
        maxn = Int(round(expn * mult_maxn))
        zlist = dfz[(dfz[!, :tau_H1] .== τ_H1) .&
                    (dfz[!, :mult_maxn] .== mult_maxn), :zval]
        OUT2 = Array{Float64}(undef, length(distrlist), 12)
        for j = 1:length(distrlist)
            OUT = Array{Float64}(undef, REPS, 12)
            for i = 1:REPS
                q = FPD_bernoulli_trial_fast(distrlist[j][1], distrlist[j][2], τ_H0,
                             τ_d, α, αc, z_β, p, goal_tau_var,
                             expn, maxn, distrlist[j][3], zlist)
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
            println(τ_H1, " ", mult_maxn, " ", j)
        end
    end
end
