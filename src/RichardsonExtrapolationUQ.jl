module RichardsonExtrapolationUQ

using Statistics
using FinEtools
using FinEtools.AlgoBaseModule: richextrapol

"""
    richextrapol_uq(solutions, elementsizes; W = 2/3)

Quantify uncertainty of Richardson extrapolation.

This function it takes exactly four solutions and their corresponding element sizes, 
and computes the Richardson extrapolated solution and convergence rate, along 
with confidence intervals for both quantities. The confidence intervals are computed 
using a weighted combination of the best extrapolated solution and the median 
of all four extrapolated solutions, with the weight `W` controlling the balance 
between these two estimates.

Reference: VVUQ-22-1017 - Confidence Intervals for 
Richardson Extrapolation in Solid Mechanics.
"""
function richextrapol_uq(solutions, elementsizes; W = 2/3)
    @assert length(solutions) == 4
    @assert length(elementsizes) == 4
    @assert length(elementsizes)  == length(solutions)
    # Four possible combinations of results
    c = [[1, 2, 3], [1, 3, 4], [1, 2, 4], [2, 3, 4]]

    edat = []
    for j in eachindex(c)
        e = (NaN, NaN, NaN, NaN)
        err = ""
        try
            e = richextrapol(solutions[c[j]], elementsizes[c[j]])
            if e[4] > minimum(abs.(solutions)) / 1.0e-6
                e = (NaN, NaN, NaN, NaN)
                err = "Richardson extrapolation failed: large residual ($e[4])"
            end
        catch
            e = (NaN, NaN, NaN, NaN)
            err = "Richardson extrapolation failed"
        end
        push!(edat, (solnestim=e[1], beta=e[2], c=e[3], maxresidual=e[4], err=err, 
                     data=(solutions=solutions[c[j]], elementsizes=elementsizes[c[j]])))
    end
    anyerr = any([e.err != "" for e in edat])
    if anyerr
        return (success=false, q_star=NaN, q_star_ci=NaN, beta_star=NaN, beta_star_ci=NaN, edat=edat)
    else
        extrsols = [e.solnestim for e in edat]
        q_m = median(extrsols)
        q_mad = 1.4826 * median(abs.(extrsols .- q_m))
        q_star = W * edat[end].solnestim + (1 - W) * q_m
        beta_m = median([e.beta for e in edat])
        beta_mad = median(abs.([e.beta for e in edat] .- beta_m))
        beta_star = W * edat[end].beta + (1 - W) * beta_m
        return (success=true, q_star=q_star, q_star_ci=2*q_mad, beta_star=beta_star, beta_star_ci=2*beta_mad, edat=edat)
    end
end 

end # module RichardsonExtrapolationUQ
