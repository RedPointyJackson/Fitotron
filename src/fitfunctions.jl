# Utility function
function propagate_errors(f,params,param_stdevs,x)
    # Given f(x,params) and the stdev in the params, compute the stdev
    # of f in x by error propagation.
    #
    # Take the derivatives of f(x,params) in the point
    g(params) = f(x,params)
    if length(params) > 1
        dds = Calculus.gradient(g)(params)
    else
        dds = g'(params)
    end
    # The stdev squared is (∂f/∂p₁ ⋅ δp₁)² + ⋯ 
    stdev = sqrt(sumabs2( dds .* param_stdevs))
end



# function fit_model_multivariate( fit_func::Function
#                                 , df::DataFrame
#                                 , initial_params::Vector{Float64}
#                    ; rescaling::Bool        = false
#                    , method::Optim.Optimizer = Optim.NelderMead()
#                    , xcolumn                 = 1
#                    , ycolumn                 = 2
#                    , yerrcolumn              = Void
#                    , xerrcolumn              = Void
#                    )

#     N, columns = size(df)

#     # First, infer x errors, y errors, x data, etc.
#     @assert columns >= xcolumn "xcolumn ($xcolumn) > columns ($columns) "
#     @assert columns >= ycolumn "ycolumn ($ycolumn) > columns ($columns) "
#     if xerrcolumn != Void
#         @assert columns >= xerrcolumn "xerrcolumn ($xerrcolumn) > columns ($columns) "
#     end
#     if yerrcolumn != Void
#         @assert columns >= yerrcolumn "yerrcolumn ($yerrcolumn) > columns ($columns) "
#     end

#     xx   = df[xcolumn]
#     yy   = df[ycolumn]
#     xerr = xerrcolumn == Void ? ones(xx) : df[xerrcolumn]
#     yerr = yerrcolumn == Void ? ones(yy) : df[yerrcolumn]


#     # Auxiliar functions
#     model(xvalues,params) = [fit_func(x,params)::Float64 for x in xvalues]

#     function cost(params)
#         # Cost function to minimize. Defined as the sum of squares of (y-yᵢ)/σᵢ
#         csts = ( model(xx,params)-yy )  ./ yerr
#         return sum(csts.^2)
#     end

#     # Find the best parameters
#     opt = Optim.optimize(cost, initial_params, method)

#     best_params = Optim.minimizer(opt)

#     # Compute the jacobian of the function which returns a vector with the
#     # fit of the model in x = xx (the user x vector), as a function of the parameters.
#     ypoint(params) = model(xx,params)

#     Jac = Calculus.jacobian(ypoint) 
#     J = Jac(best_params) # Eval at minimum

#     # Weight matrix
#     W = diagm(1./(yerr.^2))

#     # Covariance matrix
#     C = (J'*W*J)^-1

#     # Parameter stdevs
#     if rescaling
#         # Residuals should follow a χ² distribution at the minimum.
#         # The value should be in that case the d.o.f.
#         dof = N - length(initial_params)
#         dof < 1 && warn("dof = $dof < 1")
#         param_stdevs = sqrt(diag(C)*Optim.minimum(opt)/dof)
#     else
#         param_stdevs = sqrt(diag(C))
#     end

#     # Make functions to access the fit with ease
#     fit_in_point(x) = fit_func(x,best_params)
#     stdev_in_point(x) = propagate_errors(fit_func,best_params,param_stdevs,x)

#     # Return the fit in its nice container
#     return FitResult( best_params, # Fit results
#                       param_stdevs, # Standard deviations at 1 σ
#                       C, # Covariance of the parameters
#                       opt.f_minimum, # final sum of residuals
#                       fit_in_point, # gives value of fit function in x
#                       stdev_in_point) # gives stderr of fit function in x
# end


# function fit_model_univariate(fit_func::Function
#                               , df::DataFrame
#                               , lowerbound::Float64
#                               , upperbound::Float64
#                    ; rescaling::Bool         = false
#                    , method::Optim.Optimizer = Optim.Brent()
#                    , xcolumn                 = 1
#                    , ycolumn                 = 2
#                    , yerrcolumn              = Void
#                    , xerrcolumn              = Void
#                    )

#     N, columns = size(df)

#     # First, infer x errors, y errors, x data, etc.
#     @assert columns >= xcolumn "xcolumn ($xcolumn) > columns ($columns) "
#     @assert columns >= ycolumn "ycolumn ($ycolumn) > columns ($columns) "
#     if xerrcolumn != Void
#         @assert columns >= xerrcolumn "xerrcolumn ($xerrcolumn) > columns ($columns) "
#     end
#     if yerrcolumn != Void
#         @assert columns >= yerrcolumn "yerrcolumn ($yerrcolumn) > columns ($columns) "
#     end

#     xx   = df[xcolumn]
#     yy   = df[ycolumn]
#     xerr = xerrcolumn == Void ? ones(xx) : df[xerrcolumn]
#     yerr = yerrcolumn == Void ? ones(yy) : df[yerrcolumn]


#     # Auxiliar functions
#     model(xvalues,params) = [fit_func(x,params)::Float64 for x in xvalues]

#     function cost(params)
#         # Cost function to minimize. Defined as the sum of squares of (y-yᵢ)/σᵢ
#         csts = ( model(xx,params)-yy )  ./ yerr
#         return sum(csts.^2)
#     end

#     # Find the best parameters
#     opt = Optim.optimize(cost,lowerbound,upperbound,method)

#     best_param = Optim.minimizer(opt)

#     # Compute the jacobian of the function which returns a vector with the
#     # fit of the model in x = xx (the user x vector), as a function of the parameters.
#     # In this case, it's just the row vector 
#     # [∂/∂x f(x_1), ⋯ , ∂/∂x  f(x_n) ]ᵀ

#     Jac(p) = [Calculus.derivative(x->fit_func(x,p),ζ)::Float64 for ζ in xx]
#     J = Jac(best_param)

#     # Weight matrix
#     W = diagm(1./(yerr.^2))

#     # Covariance matrix is just an scalar
#     C = (J'*W*J)[1] |> inv

#     # Parameter stdevs
#     if rescaling
#         # Residuals should follow a χ² distribution at the minimum.
#         # The value should be in that case the d.o.f.
#         dof = N - 1 # Just one parameter!
#         dof < 1 && warn("dof = $dof < 1")
#         param_stdevs = sqrt(C*Optim.minimum(opt)/dof)
#     else
#         param_stdevs = sqrt(C)
#     end

#     # Make functions to access the fit with ease
#     fit_in_point(x) = fit_func(x,best_param)
#     stdev_in_point(x) = propagate_errors(fit_func,best_param,param_stdevs,x)

#     # Return the fit in its nice container
#     return FitResult( [best_param], # Fit results
#                       [param_stdevs], # Standard deviations at 1 σ
#                       C*ones(1,1), # Covariance of the parameters
#                       opt.f_minimum, # final sum of residuals
#                       fit_in_point, # gives value of fit function in x
#                       stdev_in_point) # gives stderr of fit function in x
# end


function fit_model_univariate(fit_func::Function
                              , x::Vector{Float64}
                              , y::Vector{Float64}
                              , lowerbound::Float64
                              , upperbound::Float64
                              ; xerr                       = Void
                              , yerr                       = ones(y)
                              , fitmethod::Optim.Optimizer = Optim.Brent()
                              , uncertaintymethod::Symbol  = :chisweep
                              , rescaling                  = false 
                              )

    N = length(x)
    dof = N - 1 # Just one parameter!
    dof < 1 && warn("dof = $dof < 1")

    # To save the method description
    unc_method_string = ""

    # Auxiliar functions
    model(xvalues,param) = [fit_func(x,param)::Float64 for x in xvalues]

    function cost(param)
        # Cost function to minimize. Defined as the sum of squares of (y-yᵢ)/σᵢ
        csts = ( model(x,param)-y )  ./ yerr
        return sum(csts.^2)
    end

    # Find the best parameters
    opt = Optim.optimize(cost,lowerbound,upperbound,fitmethod)

    best_param = Optim.minimizer(opt)

    if xerr != Void
        error("Orear's effective variances not implemented yet")
    end

    # Time to stimate the errors.
    if uncertaintymethod == :jacobian
        unc_method_string = "Jacobian"
        # Compute the jacobian of the function which returns a vector
        # with the fit of the model in x = xx (the user x vector), as
        # a function of the parameters.
        # In the 1D case, it's just the row vector 
        # [∂/∂x f(x_1), ⋯ , ∂/∂x  f(x_n) ]ᵀ

        Jac(p) = [Calculus.derivative(x->fit_func(x,p),ζ)::Float64 for ζ in x]
        J = Jac(best_param)

        # Weight matrix
        W = diagm(1./(yerr.^2))

        # Covariance matrix is just an scalar
        C = (J'*W*J)[1] |> inv
        C = fill(C,1,1)

        # Stdevs follow from C
        param_stdevs = [sqrt(C[1])]
    elseif uncertaintymethod == :chisweep
        unc_method_string = "χ² sweeping"
        # We won't have covariance
        C = nothing
        # The variance of χ² is just twice it's mean, or we hope so.
        # Find that point in the cost function, wich is ∼χ²:
        onestddeviation = √2*Optim.minimum(opt)
        chifit = optimize(x->abs(cost(x)-onestddeviation)
                          ,best_param[1]*(1-10)
                          ,best_param[1]*(1+10)
                          )
        parameterstdev = Optim.minimizer(chifit)
        # We have truncated the search at a maximum deviation of
        # 1000%. Warn the user if that is the case.
        if abs((parameterstdev - best_param)/best_param) > 9
            warn(
"""
The deviation of the parameter is superior to 900%.
It will be truncated at 1000%.
"""
                 )
        end
        # Return stdev
        param_stdevs = [abs(best_param - parameterstdev)]
    else
        error("Unsuported method $uncertaintymethod for the uncertainty estimation")
    end

    if rescaling
        # We need χ² to be the dof at the minimum
        # Using that χ² in the minimum is ∝ 1/yerr², it sufices to
        # multiply yerr by sqrt(χmin/dof). And that implies that the
        # σ's get multiplied by that also.
        factor = sqrt(Optim.minimum(opt)/dof)
        if C != nothing 
            C *= factor
        end
        param_stdevs *= factor
        unc_method_string = "$unc_method_string (rescaled)"
    else
        unc_method_string = "$unc_method_string (unrescaled)"
    end


    # Parameter stdevs
    return FitResult(
                       [best_param]        # Fit results
                     , param_stdevs        # Deviations found
                     , fit_func            # Function used to fit
                     , opt                 # Optim result
                     , cost                # Cost function
                     , 1                   # n of fit parameters
                     , C                   # Estimation of the covariance
                     , unc_method_string   # Uncertainty estimation method
                     , dof                 # degrees of freedom
                     )
end



# function fit_model_multivariate(fit_func::Function
#                                 , x::Vector{Float64}
#                                 , y::Vector{Float64}
#                                 , initialparams::Vector{Float64}
#                                 ; xerr                       = Void
#                                 , yerr                       = ones(y)
#                                 , fitmethod::Optim.Optimizer = Optim.NelderMead()
#                                 , rescaling                  = false 
#                                 )

#     # Auxiliar functions
#     model(xvalues,param) = [fit_func(x,param)::Float64 for x in xvalues]

#     function cost(param)
#         # Cost function to minimize. Defined as the sum of squares of (y-yᵢ)/σᵢ
#         csts = ( model(x,param)-y )  ./ yerr
#         return sum(csts.^2)
#     end

#     # Find the best parameters
#     opt = Optim.optimize(cost,initialparams,fitmethod)

#     best_param = Optim.minimizer(opt)

#     if xerr != Void
#         error("Orear's effective variances not implemented yet")
#     end

#     # Compute the jacobian of the function which returns a vector with the
#     # fit of the model in x (the user x vector), as a function of the parameters.

#     Jac = Calculus.jacobian(param->model(x,param)) 
#     J = Jac(best_param) # Eval at minimum

#     # Weight matrix
#     W = diagm(1./(yerr.^2))

#     # Covariance matrix
#     C = (J'*W*J)^-1

#     # Parameter stdevs
#     N = length(x)
#     dof = N - length(initialparams)
#     dof < 1 && warn("dof = $dof < 1")
#     if rescaling
#         unc_method_string = "Jacobian with rescaling"
#         # Residuals should follow a χ² distribution at the minimum.
#         # The value should be in that case the d.o.f.
#         param_stdevs = sqrt(diag(C)*Optim.minimum(opt)/dof)
#     else
#         unc_method_string = "Jacobian without rescaling"
#         param_stdevs = sqrt(diag(C))
#     end

#     return FitResult(
#                        best_param          # Fit results
#                      , param_stdevs        # Deviations found
#                      , fit_func            # Function used to fit
#                      , opt                 # Optim result
#                      , cost                # Cost function
#                      , 1                   # n of fit parameters
#                      , C                   # Estimation of the covariance
#                      , unc_method_string   # Uncertainty estimation method
#                      , dof                 # Degrees of freedom
#                      )

# end

