module RichardsonExtrapolationUQ

using Statistics
using FinEtools
using FinEtools.AlgoBaseModule: richextrapol

"""
    richextrapol_uq(solutions, elementsizes; W = 2/3)

Quantify uncertainty of Richardson extrapolation.

Reference: VVUQ-22-1017 - Confidence Intervals for 
Richardson Extrapolation in Solid Mechanics.
"""
function richextrapol_uq(solutions, elementsizes; W = 2/3)
    @assert length(solutions) > 3
    @assert length(elementsizes) > 3
    @assert length(elementsizes)  == length(solutions)
    # Four possible combinations of results
    c = [[1, 2, 3], [1, 3, 4], [1, 2, 4], [2, 3, 4]]

    # If we have more results than four, we can use a sliding window 
    # of four results to compute the extrapolation repeatedly.
    results = []
    for i in 1:length(solutions)-3
        ess = elementsizes[i:i+3]
        qs = solutions[i:i+3]

        edat = []
        for j in eachindex(c)
            e = (NaN, NaN, NaN, NaN, "")
            try
                e = richextrapol(qs[c[j]], ess[c[j]])
                if e[4] > minimum(abs.(solutions)) / 1.0e-6
                    e = (NaN, NaN, NaN, NaN, "Richardson extrapolation failed: large residual ($e[4])")
                end 
            catch
            end
            # println("extrapolation $(qs[c[j]]) $(e[1])")
            # println("convergence rate $(e[2])")
            push!(edat, (solnestim = e[1], beta = e[2], c = e[3], maxresidual = e[4], err = e[5], data = (qs[c[j]], ess[c[j]])))
        end
        extrsols = [e.solnestim for e in edat]
        q_m = median(extrsols)
        q_mad = 1.4826 * median(abs.(extrsols .- q_m))
        q_star = W * edat[end].solnestim + (1 - W) * q_m
        beta_m = median([e.beta for e in edat])
        beta_mad = median(abs.([e.beta for e in edat] .- beta_m))
        beta_star = W * edat[end].beta + (1 - W) * beta_m
        push!(results, (q_star = q_star, q_star_ci = 2*q_mad, beta_star = beta_star, beta_star_ci = 2*beta_mad, elementsize = ess[end-1], edat = edat)) 
    end 
    return results
end 

end # module RichardsonExtrapolationUQ
