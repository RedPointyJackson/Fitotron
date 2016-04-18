module Fitotron

    export fit_model

    using Optim
    using Calculus
    using DataFrames
    import Base.show

    # Container of fit results
    immutable FitResult
        # expand later. We'll need dof and things like that.
        param_results::Vector{Float64} # Fit results
        param_stdevs::Vector{Float64} # Standard deviations at 1 σ
        covariance::Matrix{Float64} # Covariance of the parameters
        resid::Float64 # final sum of residuals
        fit_value::Function # gives value of fit function in x
        fit_stdev::Function # gives stderr of fit function in x
    end


    function show(io::IO, r::FitResult)
        println(io, "Fit results:")
        println(io, "----------------------")
        N = length(r.param_results)
        for i in 1:N
            @printf(io,"Param. %d: %.2e ± %.2e\n",i,r.param_results[i],r.param_stdevs[i])
        end
        println(io, "S² at the minimum: ",r.resid)
    end

    include("functions.jl")

end # module
