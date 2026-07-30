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

    results = []
    for j in 1:length(solutions)-3
        ess = elementsizes[j:j+3]
        qs = solutions[j:j+3]

        extrapolations = []
        for j in eachindex(c)
            e = (Inf * sign(qs[end]), 0.0, 0.0, Inf)
            try
                e = richextrapol(qs[c[j]], ess[c[j]])
                if e[4] > minimum(abs.(solutions)) / 1.0e-6
                    error("Richardson extrapolation failed: large residual ($e[4])")
                end 
            catch
            end
            # println("extrapolation $(qs[c[j]]) $(e[1])")
            # println("convergence rate $(e[2])")
            push!(extrapolations, (solnestim = e[1], beta = e[2], c = e[3], maxresidual = e[4]))
        end
        extrsols = [e.solnestim for e in extrapolations]
        q_m = median(extrsols)
        q_m_ad = 1.4826 * median(abs.(extrsols .- q_m))
        q_star = W * extrapolations[end].solnestim + (1 - W) * q_m
        beta_m = median([e.beta for e in extrapolations])
        beta_m_ad = median(abs.([e.beta for e in extrapolations] .- beta_m))
        beta_star = W * extrapolations[end].beta + (1 - W) * beta_m
        push!(results, (estim = q_star, estim_ad_x_2 = 2*q_m_ad, beta = beta_star, beta_ad_x_2 = 2*beta_m_ad, elementsize = ess[end-1], extrapolations = extrapolations)) 
    end 
    return results
end 

end # module RichardsonExtrapolationUQ
