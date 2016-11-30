"""

    fit(LinearModel)

Fit to a straight line.
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
