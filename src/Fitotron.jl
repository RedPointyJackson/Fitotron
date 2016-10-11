__precompile__()

module Fitotron

export fitmodel, plotfit

using Optim
using Calculus
using DataFrames
using Gadfly
import Base.show

# Constants
MAX_PARAM_DEV = 9

# Avoid big function signatures
MaybeFMatrix   = Union{Matrix{Float64}        , Void}
MaybeFVector   = Union{Vector{Float64}        , Void}
MaybeNVector   = Union{Vector{Number}         , Void}
MaybeOptimizer = Union{Vector{Optim.Optimizer}, Void}
MaybeBool      = Union{Bool                   , Void}
MaybeSymbol    = Union{Symbol                 , Void}

# Container of fit results
immutable FitResult
    data::Matrix{Float64}                 # Columns with x,y,yerr[,xerr]
    param_results::Vector{Float64}        # Fit results
    param_deviations::Vector{Float64}     # Deviations found
    fit_func::Function                    # Function used to fit
    optresults::Optim.OptimizationResults # Optim result
    cost::Function                        # Cost function
    nparams::Int64                        # n of fit parameters
    covariance_matrix::MaybeFMatrix       # Covariance (can be empty)
    uncertainty_method::String            # Uncertainty estimation method
    dof::Int64                            # Degrees of freedom
end

function show(io::IO, r::FitResult)
    println(io, "Fit results:")
    println(io, "─────────────────────────────────────────────────────────")
    N = length(r.param_results)
    for i in 1:N
        param_mean = r.param_results[i]
        param_stdev = r.param_deviations[i]
        @printf(io,"Param. %d:\t\t\t%.2e ± %.2e \n"
                ,i , param_mean, param_stdev)
    end
    redchisqr = Optim.minimum(r.optresults)/r.dof
    println(io, "Reduced χ²:\t\t\t", redchisqr)
    println(io, "Parameter estimation method:\t", Optim.method(r.optresults))
    println(io, "Uncertainty estimation method:\t", r.uncertainty_method)
end

include("fitfunctions.jl")
include("fitutils.jl")

end # module
