# Utility function
function propagate_errors(f,params,param_stdev,x)
    # Given f(x,params) and the stdev in the params, compute the stdev
    # of f in x by error propagation.
    #
    # Take the derivatives of f(x,params) in the point
    g(params) = f(x,params)
    dds = Calculus.gradient(g)(params)
    # The stdev squared is (∂f/∂p₁ ⋅ δp₁)² + ⋯ 
    stdev = sqrt(sumabs2( dds .* param_stdev))
end


# Core of the library
function fit_model(fit_func,df,initial_params)
    N = size(df)[1]
    haserrors = size(df)[2] == 3 # Identify if the fit has uncertainties

    model(xx,params) = [fit_func(x,params)::Float64 for x in xx]

    xx = df[1]
    yy = df[2]
    if haserrors
        ss = df[3] # These are the stdevs in y
    else
        ss = ones(N) 
    end

    function cost(params)
        # Cost function to minimize. Defined as the sum of squares of (y-yᵢ)/σᵢ
        csts = ( model(xx,params)-yy )  ./ ss
        return sum(csts.^2)
    end

    # Find the best parameters
    opt = optimize(cost,initial_arams)

    best_params = opt.minimum

    # Compute the jacobian of the function that returns a vector with the
    # fit of the model in x = xx (the user x vector), as a function of the parameters.
    ypoint(params) = [model(xx,params)[i]::Float64 for i in 1:N]

    J = Calculus.jacobian(ypoint)(opt.minimum)

    # Weight matrix
    W = diagm(1./(ss.^2))

    # Covariance matrix
    C = (J'*W*J)^-1

    # Parameter stdevs
    param_stdev = sqrt(diag(C))

    # Make functions to access the fit with ease
    fit_in_point(x) = fit_func(x,best_params)
    stdev_in_point(x) = propagate_errors(fit_func,best_params,param_stdev,x)

    # Return the fit in its nice container
    return FitResult( best_params, # Fit results
                      param_stdev, # Standard deviations at 1 σ
                      C, # Covariance of the parameters
                      opt.f_minimum, # final sum of residuals
                      fit_in_point, # gives value of fit function in x
                      stdev_in_point) # gives stderr of fit function in x
end
