function eval_fun(n, ρ, α, one_sided)
    # derivative of the boundary function
    if one_sided
        return n^2*ρ^2*((2*n*ρ^2 + 2)*log(1 + (n*ρ^2 + 1)^0.5/(2*α))/(n^2*ρ^2))^0.5*(
               2.0*log(1 + (n*ρ^2 + 1)^0.5/(2*α))/(n*ρ) - 1.0*(2*n*ρ^2 + 2)*log(1 + (n*ρ^2 + 1)^0.5/(2*α))/(n^2*ρ^3) + 
               0.25*(2*n*ρ^2 + 2)/(α*n*ρ*(1 + (n*ρ^2 + 1)^0.5/(2*α))*(n*ρ^2 + 1)^0.5))/(
               (2*n*ρ^2 + 2)*log(1 + (n*ρ^2 + 1)^0.5/(2*α)))
    else
        return n^2*ρ^2*((2*n*ρ^2 + 2)*log((n*ρ^2 + 1)^0.5/α)/
               (n^2*ρ^2))^0.5*(0.5*(2*n*ρ^2 + 2)/(n*ρ*(n*ρ^2 + 1)^1.0) +
               2.0*log((n*ρ^2 + 1)^0.5/α)/(n*ρ) -
               1.0*(2*n*ρ^2 + 2)*log((n*ρ^2 + 1)^0.5/α)/(n^2*ρ^3))/
               ((2*n*ρ^2 + 2)*log((n*ρ^2 + 1)^0.5/α))
    end
end

function bin_search(uv, lv, tol, n, α, one_sided)
    #Find optimal ρ for a given n
    @assert eval_fun(n, lv, α, one_sided) + tol < 0
    @assert eval_fun(n, uv, α, one_sided) - tol > 0
    while true
        nt = (lv + uv) / 2
        v = eval_fun(n, nt, α, one_sided)
        if abs(v) < tol
            return nt
        elseif v < 0
            lv = nt
        elseif v > 0
            uv = nt
        end
    end
end

function print_time(rep_nr, tot_reps, start_time)
    """Output information on how many simulations that have been performed
    and when the simulation is expected to finish."""
    elapsed_time = time_ns() - start_time
    perc_compl = (rep_nr / tot_reps) * 100
    exp_tot_time = elapsed_time / (perc_compl / 100)
    eft = Dates.now() + Dates.Nanosecond(floor(exp_tot_time - elapsed_time))
    perc_compl = round(perc_compl, sigdigits=3)
    println(perc_compl, "% completed, expect to finish ", eft)
end

function initialize_data(distrY0, distrY1, p, shift, a)
    n0 = 0
    n1 = 0
    W = []
    Y = []
    while (n0 < 2) | (n1 < 2)
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
    n = n0 + n1
    return W, Y, n0, n1, n
end

function get_phi(ρ, α, n)
    ϕ = sqrt(
        2 * (n * ρ^2 + 1) / (n^2 * ρ^2) * 
        log(1 + sqrt(n * ρ^2 + 1) / (2 * α))
    )
    return ϕ
end

function get_psi(ρ, α, n)
    ψ = sqrt(
        2 * (n * ρ^2 + 1) / (n^2 * ρ^2) * 
        log(sqrt(n * ρ^2 + 1) / α)
    )
    return ψ
end

function simulation_distributions()
    return [[Normal(0,1), Normal(0,1), 0],
            [LogNormal(0,1), LogNormal(0,1), 0],
            [Normal(0,1), Normal(0.2,1) , 0],
            [LogNormal(0,1), LogNormal(0,1), .2],
            [LogNormal(1,1), LogNormal(1.2660058352236003, 1), 0],
            [Gamma(1,1), LogNormal(-0.09498792407321321, 3/4), 0],
            [Gamma(1,1), Gamma(1,2/(2-.2^2) + sqrt(4 / (2-.2^2)^2 - 1)), 0],
            [Gamma(1,1), Gamma(1.2102498439450078, 1), 0]
           ]
end