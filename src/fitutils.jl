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
    # Fit function layer
    fitlayer = layer(fit.fit_func
                     ,xmin, xmax)
    # Fit deviation layer
    error_vals = map(fit.fit_dev, xrange)
    fit_vals = map(fit.fit_func, xrange)
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

    plotcost(res)

Returns a Gadfly contour plot of the cost function of `res`. It must
have two parameters maximum. The 1σ, 2σ and 3σ contours will be drawn,
with the obtained parameters.

TODO: when sigma contours are outside the picture, it fails.
"""
function plotcost(fit::FitResult)
    Gadfly.push_theme(:dark)
    if length(fit.param_results) == 1
        error("TODO for 1D fits")
    elseif length(fit.param_results) == 2
        μ1, μ2 = fit.param_results
        σ1, σ2 = fit.param_deviations
        sigmarange = 6
        xrange = linspace(μ1-sigmarange*σ1, μ1+sigmarange*σ1, 100)
        yrange = linspace(μ2-sigmarange*σ2, μ2+sigmarange*σ2, 100)
        fitcost(x,y) = fit.cost([x; y])
        mincost = fitcost(μ1,μ2)
        colors = linspace(nicered,niceblue,4)
        p =
        plot(z=fitcost, x=xrange, y=yrange
             , Geom.contour(
                            # Variance of red. χ² is 2N → σ = √(2N)
                            levels = mincost + [0.5;1;2;3]*sqrt(2*mincost)/2
                            )
             , Scale.color_discrete_manual(colors...)
             , Guide.manual_color_key("Coberture"
                                      ,["½σ";"1σ";"2σ";"3σ"]
                                      ,collect(colors)
                                      )
             ,Guide.annotation(
                              compose(context()
                                      , circle(μ1,μ2,1mm)
                                      , line([(μ1,μ2),(μ1+σ1,μ2)])
                                      , line([(μ1,μ2),(μ1-σ1,μ2)])
                                      , line([(μ1,μ2),(μ1,μ2+σ2)])
                                      , line([(μ1,μ2),(μ1,μ2-σ2)])
                                      , fill(colorant"orange"),
                                      stroke(colorant"orange")
                                      )
                               )
             , Guide.xlabel("Parameter 1")
             , Guide.ylabel("Parameter 2")
             , Guide.title("Cost function")
             )
    else
        error("Fit has too many parameters for a 2D screen.")
    end
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

"""
    thing2array(anything)

Gets something and tries its best to obtain a vector of Float64
"""
function thing2array(thing)::Vector{Float64}
    if length(thing) == 0
        error("Can't use an empty array")
    else
        return vec([Float64(x) for x in thing])
    end
end
