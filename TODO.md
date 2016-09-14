- The jacobian method is nice because it gives a covariance, but fails
  giving uncertainties for constant functions. MWE:
  
    x = Float64[1,3,4,7]
    y = [2.3, 2.1, 2.0, 1.8]
    yerr = 0.2*ones(4)
    f(x,p) = p
    # Create the fit
    fit = fit_model_univariate(f, x, y, 0.0, 10.0
				, yerr = yerr)
  
  The uncertainty is Inf.
 


- Implement Orear's method

- Be consistent on notation (foobar vs foo_bar)

- Add test for the uncertainty estimation
