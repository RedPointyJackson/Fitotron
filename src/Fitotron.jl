__precompile__()

module Fitotron

    export fit_model_multivariate, fit_model_univariate

    using Optim
    using Calculus
    using DataFrames
    import Base.show

    MaybeMatrix = Union{Matrix{Float64},Void}

    # Container of fit results
    immutable FitResult
        param_results::Vector{Float64}        # Fit results
        param_deviations::Vector{Float64}     # Deviations found
        fit_func::Function                    # Function used to fit
        optresults::Optim.OptimizationResults # Optim result
        cost::Function                        # Cost function
        nparams::Int64                        # n of fit parameters
        covariance_matrix::MaybeMatrix        # Estimation of the
                                              # covariance (can be empty)
        uncertainty_method::UTF8String        # Uncertainty estimation method
        dof::Int64                            # Degrees of freedom
    end

    function show(io::IO, r::FitResult)
        println(io, "Fit results:")
        println(io, "────────────────────────────────────────────────")
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

end # module
