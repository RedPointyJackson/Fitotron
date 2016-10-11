# Utility function
"""

    propagate_errors(f,params,param_stdevs,x)


Given `f(x,params)` and the stdev in the params, compute the stdev
of `f` in `x` by error propagation.
"""
function propagate_errors(f,params,param_stdevs,x)
    # Take the derivatives of f(x,params) in the point
    g(params) = f(x,params)
    if length(params) > 1
        dds = Calculus.gradient(g)(params)
    else
        dds = g'(params)
    end
    # The stdev squared is ∑(∂f/∂pᵢ ⋅ δpᵢ)² + ⋯
    stdev = sqrt(sumabs2( dds .* param_stdevs))
end


"""

    fit_line(x,y,yerr)

Fit `x`,`y` arrays to a straight line, with errors
given by `yerr`. Returns m,n and the covariance
matrix, where the line of fit is mx+n.
"""
function fit_line(x,y,yerr)
    N = length(x)

    # Compute the parameters
    function br(f)
        # By definition, [x] = 1/N Σ x/σ
        if f != :one
            vals = [f[i]/yerr[i]^2 for i in 1:N]/N
            return sum(vals)
        else
            myf = [1.0 for i in 1:N]
            vals = [myf[i]/yerr[i]^2 for i in 1:N]/N
            return sum(vals)
        end
    end

    # Linear regresion parameters
    m_mean = ( br(:one)*br(x.*y) - br(x)*br(y) )/(  br(:one)*br(x.*x) - br(x)*br(x)  )
    n_mean = mean(y)-m_mean*mean(x)

    # Covariance matrix
    covmatrix = zeros(2,2)

    D = br(x.*x)*br(:one) - br(x)*br(x)

    covmatrix[2,2] = 1/(N*D) * br(x.*x)
    covmatrix[1,2] = -1/(N*D) * br(x)
    covmatrix[2,1] = -1/(N*D) * br(x)
    covmatrix[1,1] = 1/(N*D) * br(:one)

    return m_mean, n_mean, covmatrix
end

"""

    fitmodel_raw(...)


General invocation. Usually, wrappers are used.
Arguments needed (in order, compulsory):

fit_func            Function to fit.
x                   x values.
y                   y values.
initialparams       Initial params or 1D range.
xerr                x errors or `nothing`.
yerr                Same for y errors.
fitmethod           Optimizer to use.
uncmethod           Or :chisweep or :Jacobian
rescaling           Boolean.
fit_kind            :multivariate or :univariate

Should not be used except for developing.
"""
function fitmodel_raw(fit_func::Function
                       , x::Vector{Float64}
                       , y::Vector{Float64}
                       , initialparams::Vector{Float64}
                       , xerr::MaybeFVector
                       , yerr::Vector{Float64}
                       , fitmethod
                       , uncmethod::Symbol
                       , rescaling::Bool
                       , fit_kind::Symbol
                       )
    const N = length(x)

    if fit_kind == :univariate
        const nparameters = 1
    else
        const nparameters = length(initialparams)
    end

    dof = N - nparameters # Just one parameter!

    # TODO check how many do you need. More than
    # one for sure...
    dof < 1 && warn("dof = $dof < 1")

    # To save the method description
    unc_method_string = ""

    # Auxiliar functions
    # The [1] fixes a bug with functions of the style p*x instead of
    # p[1]*x.
    model(xvalues,param) =
    [fit_func(x,param)[1]::Float64 for x in xvalues]

    function cost(param)
        # Cost function to minimize. Defined as
        # the sum of squares of (y-yᵢ)/σᵢ
        csts = ( model(x,param)-y )  ./ yerr
        return sum(csts.^2)
    end

    # Find the best parameters
    if fit_kind == :univariate
        lowerbound, upperbound = initialparams
        opt = Optim.optimize(cost,lowerbound,upperbound,fitmethod)
    elseif fit_kind == :multivariate
        opt = Optim.optimize(cost,initialparams,fitmethod)
    end

    best_param = Optim.minimizer(opt)
    # For consistency, ensure it's a vector
    if typeof(best_param) != Vector{Float64}
        best_param = [best_param]
    end

    # If needed, Orear things
    if xerr != nothing
        error("Orear's effective variances not implemented yet")
    end

    # Time to estimate the errors.
    if uncmethod == :jacobian
        unc_method_string = "Jacobian"
        # Compute the jacobian of the function which returns a vector
        # with the fit of the model in x = xx (the user x vector), as
        # a function of the parameters.
        if fit_kind == :univariate
            # In the 1D case, it's just the row vector
            # [∂/∂x f(x_1), ⋯ , ∂/∂x  f(x_n) ]ᵀ
            Jac(p) = [Calculus.derivative(x->fit_func(x,p),ζ)::Float64 for ζ in x]
            J = Jac(best_param[1]) # Eval at minimum
        elseif fit_kind == :multivariate
            Jac = Calculus.jacobian(param->model(x,param))
            J = Jac(best_param)
        end
        # Weight matrix
        W = diagm(1./(yerr.^2))
        # Covariance
        if fit_kind == :univariate
            C = 1/(J'*W*J)[1]
            param_stdevs = [sqrt(C)]
            # Make the covariance a matrix
            C = fill(C,1,1)
        elseif fit_kind == :multivariate
            C = (J'*W*J) |> inv
            param_stdevs = [sqrt(c) for c in diag(C)]
        end
    elseif uncmethod == :chisweep
        unc_method_string = "χ² sweeping"
        # We won't have covariance
        C = nothing
        # The variance of a χ² distributed thing
        # is just twice it's mean, or we hope so.
        param_stdevs = Vector{Float64}(nparameters)
        for i in 1:nparameters
            # Find that point in the cost function,
            # wich is ≃ χ²:
            onestddeviation = √2*Optim.minimum(opt)
            bounds = [best_param[i]*(1-10)
                      best_param[i]*(1+10)]
            sort!(bounds)
            # Build a function constant in all
            # params except the explored one
            function explorer(p)
                v = copy(best_param)
                v[i] = p
                return abs(cost(v)-onestddeviation)
            end
            chifit = optimize(explorer
                              , bounds[1]
                              , bounds[2]
                              )
            # Now we now the variance, which is
            # twice the value found minus the
            # "mean" point:
            paramdev = Optim.minimizer(chifit)
            parameterstdev =
            sqrt(2*abs(paramdev - best_param[i]))
            # We have truncated the search at a maximum deviation of
            # 1000%. Warn the user if that is the case.
            if abs((paramdev - best_param[i])/best_param[i]) > MAX_PARAM_DEV
                warn(
"""
The deviation of the $(i)th parameter is superior to 900%.
It will be truncated at 1000%.
"""
                     )
            end
            # Return stdev
            param_stdevs[i] = abs(best_param[i]-parameterstdev)
        end
    else
        error(
"Unsuported method $uncmethod for the
uncertainty estimation"
              )
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

    if xerr == nothing
        data = [x y yerr]
    else
        data = [x y yerr xerr]
    end

    # Parameter stdevs
    return FitResult(
                       data                # Columns with x,y[,yerr,xerr]
                     , best_param          # Fit results
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


"""

    fitmodel(f, x, y, from, to)

Fit a 1D model. Optional parameters are:
  xerr               errors in x
  yerr               errors in y
  fitmethod          Custom `Optim.Optimizer()` to
                     use. Defaults to Brent().
  uncmethod          Or `:chisweep` (default) or `:jacobian`
  rescaling          Force red.χ² = 1 in the minimum. Enabled
                     by default if there are no x or y errors.
"""
function fitmodel(fit_func::Function
                              , x
                              , y
                              , from
                              , to
                              ; xerr = nothing
                              , yerr = nothing
                              , fitmethod = nothing
                              , uncmethod::MaybeSymbol = nothing
                              , rescaling::MaybeBool = nothing
                              )
    fit_kind = :univariate

    # Enable rescaling if needed
    if rescaling == nothing
        if xerr == nothing && yerr == nothing
            rescaling = true
        else
            rescaling = false
        end
    end
    # We are in a 1D fit
    initialparams = [Float64(from), Float64(to)]

    # Everything to Float64
    x = [Float64(ζ) for ζ in x]
    y = [Float64(ζ) for ζ in y]

    if xerr != nothing
        xerr = [Float64(ζ) for ζ in xerr]
    end

    if yerr != nothing
        yerr = [Float64(ζ) for ζ in yerr]
    else
        yerr = ones(x)
    end

    # Get an optimizer and an uncertainty method
    if fitmethod == nothing
        fitmethod = Optim.Brent()
    end
    if uncmethod == nothing
        uncmethod = :chisweep
    end


    # Invoke the main function
    fitmodel_raw(fit_func::Function
                  , x::Vector{Float64}
                  , y::Vector{Float64}
                  , initialparams::Vector{Float64}
                  , xerr::MaybeFVector
                  , yerr::Vector{Float64}
                  , fitmethod
                  , uncmethod::Symbol
                  , rescaling::Bool
                  , fit_kind::Symbol
                  )
end


"""

    fitmodel(f, x, y, [initalparams])

Fit a ND model, with N>1. Optional parameters are:
  xerr               errors in x
  yerr               errors in y
  fitmethod          Custom `Optim.Optimizer()` to
                     use. Defaults to NelderMead()
  uncmethod          Or `:chisweep` (default) or `:jacobian`
  rescaling          Force red.χ² = 1 in the minimum. Enabled
                     by default if there are no x or y errors.
"""
function fitmodel(fit_func::Function
                       , x
                       , y
                       , initialparams
                       ; xerr = nothing
                       , yerr = nothing
                       , fitmethod = nothing
                       , uncmethod::MaybeSymbol = nothing
                       , rescaling::MaybeBool = nothing
                       )
    fit_kind = :multivariate

    # Enable rescaling if needed
    if rescaling == nothing
        if xerr == nothing && yerr == nothing
            rescaling = true
        else
            rescaling = false
        end
    end


    # Everything to Float64
    x = [Float64(ζ) for ζ in x]
    y = [Float64(ζ) for ζ in y]

    if xerr != nothing
        xerr = [Float64(ζ) for ζ in xerr]
    end

    if yerr != nothing
        yerr = [Float64(ζ) for ζ in yerr]
    else
        yerr = ones(x)
    end

    initialparams = [Float64(ζ) for ζ in initialparams]

    # Get an optimizer and an uncertainty method
    if fitmethod == nothing
        fitmethod = Optim.NelderMead()
    end
    if uncmethod == nothing
        uncmethod = :chisweep
    end


    # Invoke the main function
    fitmodel_raw(fit_func::Function
                  , x::Vector{Float64}
                  , y::Vector{Float64}
                  , initialparams::Vector{Float64}
                  , xerr::MaybeFVector
                  , yerr::Vector{Float64}
                  , fitmethod
                  , uncmethod::Symbol
                  , rescaling::Bool
                  , fit_kind::Symbol
                  )

end
