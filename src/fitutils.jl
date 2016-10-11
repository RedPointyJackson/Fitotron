using Gadfly
using Fitotron


"""

    plotfit(res)

Returns a Gadfly plot of `res`.
"""
function plotfit(fit::FitResult)
    # Customization of the plot style
    Gadfly.push_theme(:dark)
    # Number of points when discretizing the functions.
    const fun_resolution = 1000
    # Ignore Orear things.
    if size(fit.data,2) == 4
        error("Can't plot (yet) a fit with Orear method.")
    end
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
