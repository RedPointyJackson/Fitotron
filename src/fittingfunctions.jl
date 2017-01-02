"""

    fit(LinearModel)

Fit to a straight line using analytical results.
"""
function fitmodel(m::LinearModel)
    # Extract data
    x    = m.x
    y    = m.y
    yerr = m.yerr

    L = length(x)

    # Precompute things
    brone = br(ones(L), yerr)
    brx   = br(x      , yerr)
    brxx  = br(x.*x   , yerr)
    bry   = br(y      , yerr)
    brxy  = br(x.*y   , yerr)

    # Means
    m_mean = ( brone*brxy - brx*bry )/
             ( brone*brxx - brx*brx )
    n_mean = mean(y) - m_mean * mean(x)

    # Covariance matrix
    covmatrix = zeros(2,2)

    D = brxx*brone - brx*brx

    covmatrix[2,2] = +1/(L*D) * brxx
    covmatrix[1,2] = -1/(L*D) * brx
    covmatrix[2,1] = -1/(L*D) * brx
    covmatrix[1,1] = +1/(L*D) * brone

    # If needed, rescale
    fit_func(x)     = m_mean*x + n_mean
    fit_func_p(x,p) = p[1]*x + p[2]
    cost(p)         = sumabs2((y[i]-fit_func_p(x[i],p))/yerr[i] for i in 1:L)
    dof             = L-2
    redχ²           = cost([m_mean;n_mean])/dof
    if m.rescale
        covmatrix *= redχ²/dof
    end
    m_dev, n_dev = covmatrix |> diag |> sqrt

    fit_dev(x)      = sqrt( (m_dev*x)^2 + n_dev^2 )

    FitResult(
               [x y yerr]       # Columns with x,y,yerr.
              ,[m_mean; n_mean] # Fit results.
              ,[m_dev; n_dev]   # Deviations found.
              ,fit_func         # Function used to fit.
              ,fit_dev          # 1σ deviation at each point.
              ,cost             # Cost function.
              ,covmatrix        # Covariance (can be empty).
              ,dof              # Degrees of freedom.
              ,m.rescale        # Was rescaling applied?
              )
end


"""

    fit(CustomModel [;from=:beginning,to=:end])

Fit to a custom function, optimizing the cost function and estimating
the deviations with a Jacobian linearization. Alternatively, you can
specify the x range to fit with `from` and `to`. `:beginning` and
`:end` will be substituted by the minimum and maximum x.
"""
function fitmodel(m::CustomModel; from=:beginning, to=:end)

    if isa(from,Real)
        firstx = convert(Float64,from)
    elseif from == :beginning
        firstx = minimum(m.x)
    else
        error("Please, use `:beginning` or a real number for `from`")
    end


    if isa(to,Real)
        lastx = convert(Float64,to)
    elseif to == :end
        lastx = maximum(m.x)
    else
        error("Please, use `:end` or a real number for `to`")
    end

    firstidx = findfirst(x->x==firstx ,m.x)
    lastidx  = findfirst(x->x==lastx  ,m.x)

    x           = m.x[firstidx:lastidx]
    y           = m.y[firstidx:lastidx]
    yerr        = m.yerr
    f           = m.func
    bounds      = m.bounds
    rescale     = m.rescale
    dims        = m.dims
    optimizator = m.optimizator

    N = length(x)

    dof = N - dims
    dof < 1 && warn("dof = $dof < 1")

    # Auxiliar functions
    # The [1] fixes a bug with functions of the style p*x instead of
    # p[1]*x.
    vmodel(xvalues,param) =
        [f(x,param)[1]::Float64 for x in xvalues]

    # Cost function to minimize. Defined as
    # the sum of squares of (y-yᵢ)/σᵢ
    cost(param) = sumabs2( ( vmodel(x,param)-y ) ./ yerr)

    # Find the best parameters
    if dims == 1
        lowerbound, upperbound = bounds
        opt = Optim.optimize(cost,lowerbound,upperbound,m.optimizator)
    else
        opt = Optim.optimize(cost,bounds,optimizator)
    end

    param_results = Optim.minimizer(opt)
    # For consistency, ensure it's a vector
    if typeof(param_results) != Vector{Float64}
        param_results = [param_results]
    end

    # Now we have the mean values! To obtain deviations, use a
    # Jacobian linearization.

    # Compute the jacobian of the function which returns a vector
    # with the fit of the model in x = xx (the user x vector), as
    # a function of the parameters.
    if dims == 1
        # In the 1D case, it's just the row vector
        # [∂/∂x f(x_1), ⋯ , ∂/∂x  f(x_n) ]ᵀ
        Jac(p) = [Calculus.derivative(x->fit_func(x,p),ζ)::Float64 for ζ in x]
        J = Jac(best_param[1]) # Eval at minimum
    else
        Jac = Calculus.jacobian(param->vmodel(x,param))
        J = Jac(param_results)
    end
    # Weight matrix
    W = diagm(1./(yerr.^2))
    # Covariance
    if dims == 1
        C = 1/(J'*W*J)[1]
        param_deviations = [sqrt(C)]
        # Make the covariance a matrix
        C = fill(C,1,1)
    else
        C = (J'*W*J) |> inv
        param_deviations = [sqrt(c) for c in diag(C)]
    end

    if rescale
        # We need χ² to be the dof at the minimum
        # Using that χ² in the minimum is ∝ 1/yerr², it sufices to
        # multiply yerr by sqrt(χmin/dof). And that implies that the
        # σ's get multiplied by that also.
        α = sqrt(Optim.minimum(opt)/dof)
        C *= α
        param_deviations *= α
    end

    # Get some functions and export everything.
    fit_func(x) = f(x,param_results)
    fit_dev(x) = propagate_errors(f,param_results,param_deviations,x)
    FitResult(
               [x y yerr]       #Columns with x,y,yerr.
              ,param_results    #Fit results.
              ,param_deviations #Deviations found.
              ,fit_func         #Function used to fit.
              ,fit_dev          #1σ deviation at each point.
              ,cost             #Cost function.
              ,C                #Covariance (can be empty).
              ,dof              #Degrees of freedom.
              ,rescale          #Was rescaling applied?
              )
end
