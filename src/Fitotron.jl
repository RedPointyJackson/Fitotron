__precompile__()

module Fitotron

export CustomModel, LinearModel
export fitmodel, plotfit, plotcost

using Optim
using Calculus
using DataFrames

# Cool plots
using Gadfly
using Compose
using Colors
# Will be used in plots
nicered = colorant"#7E273E"
niceblue = colorant"#006291"

import Base.show

"""
    Container of fit results.
"""
immutable FitResult
    data              ::Matrix{Float64} # Columns with x,y,yerr.
    param_results     ::Vector{Float64} # Fit results.
    param_deviations  ::Vector{Float64} # Deviations found.
    fit_func          ::Function        # Function used to fit.
    fit_dev           ::Function        # 1σ deviation at each point.
    cost              ::Function        # Cost function.
    covariance_matrix ::Matrix{Float64} # Covariance (can be empty).
    dof               ::Int64           # Degrees of freedom.
    rescaling         ::Bool            # Was rescaling applied?
end

function show(io::IO, r::FitResult)
    println(io, "Fit results:")
    println(io, "───────────────────────────────────────────────────────────────")
    N = length(r.param_results)
    for i in 1:N
        param_mean = r.param_results[i]
        param_stdev = r.param_deviations[i]
        param_relerror = 100*param_stdev/param_mean |> abs
        @printf(io,"Param. %d:\t\t\t%+.2e ± %.2e (%.1f%%)\n"
                ,i , param_mean, param_stdev, param_relerror)
    end
    redchisqr = r.cost(r.param_results)/r.dof
    @printf(io, "Reduced χ²\t\t\t%.4lf\n", redchisqr)
    if r.rescaling
        @printf(io, "Errors were rescaled so χ²=d.o.f.\n")
    end
end

include("fitutils.jl")
include("models.jl")
include("fittingfunctions.jl")

end # module
