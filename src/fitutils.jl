"""

    plotfit(res)

Returns a Gadfly plot of `res`.
"""
function plotfit(fit::FitResult)
    # Customization of the plot style
    Gadfly.push_theme(:dark)
    # Number of points when discretizing the functions.
    const fun_resolution = 1000
    # Extract the data
    x      = fit.data[:,1]
    y      = fit.data[:,2]
    yerr   = fit.data[:,3]
    xmin   = minimum(x)
    xmax   = maximum(x)
    xrange = linspace(xmin,xmax,fun_resolution)
    fitfun = fit.fit_func
    params = fit.param_results
    devs   = fit.param_deviations
    # Fit function layer
    fitlayer = layer(x->fitfun(x,params)
                     ,xmin, xmax)
    # Fit deviation layer
    errorfun(x) = Fitotron.propagate_errors(fitfun,params,devs,x)
    error_vals = map(errorfun, xrange)
    fit_vals = map(x->fitfun(x,params), xrange)
    error_df = DataFrame(  x    = xrange
                         , ymax = fit_vals + error_vals
                         , ymin = fit_vals - error_vals
                           )
    devlayer = layer(error_df,x=:x,ymax=:ymax,ymin=:ymin
                     , Geom.ribbon)
    # Data layer
    if !all(x->x==1,yerr)
        # Errorbars
        ymax = y+yerr/2
        ymin = y-yerr/2
        datalayer = layer( x    = x
                          ,y    = y
                          ,ymax = ymax
                          ,ymin = ymin
                          ,Geom.point
                          ,Geom.errorbar
                          )
    else
        # Raw dots
        datalayer = layer(x=x, y=y, Geom.point)
    end
    # Combine them and plot.
    p = plot(datalayer, fitlayer, devlayer)
    Gadfly.pop_theme()
    return p
end


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
br(x,σ) computes [x] = 1/N Σ xᵢ/σ.
"""
br(x,σ) = sum(  x./(σ.^2) / length(x) )
