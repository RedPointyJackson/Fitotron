using Base.Test
using Fitotron


info("Error propagation tests")
# Error propagation, 1D, handmade example
f(x,p) = x*p
@test Fitotron.propagate_errors(f,1.0,π,0) ≈ 0.0
f(x,p) = x^p
@test Fitotron.propagate_errors(f,1.0,2.0,3.0) ≈ log(3.0)*exp(1.0*log(3.0))*2.0 

# Error propagation, 2D, handmade example
f(x,p) = x*p[1]+p[2]
@test Fitotron.propagate_errors(f,ones(2),ones(2),0) ≈ 1.0
@test Fitotron.propagate_errors(f,ones(2),ones(2),1) ≈ √2

# Error propagation, 3D, handmade example
f(x,p) = x*p[1]+p[2]*p[3]
result(x,p,δp) = sqrt( (x*δp[1])^2 + (p[3]*δp[2])^2 + (p[2]*δp[3])^2 )
@test Fitotron.propagate_errors(f,ones(3),ones(3),0) ≈ result(0,ones(3),ones(3))
@test Fitotron.propagate_errors(f,Float64[1,2,3],0.5*ones(3),0.2) ≈ result(0.2,Float64[1,2,3],0.5*ones(3))


info("Fitting tests")
# Fit a 1D function, handmade example
x = 0.001:0.1:1 |> collect 
# We don't use 0 because finite differences will use -0.001 and f
# isn't defined there
y = x.^4.2
f(x,p) = x^p

# Simple fit
res = fit_model_univariate(f, x, y, -10.0, 100.0).param_results
@test norm(res-4.2) < 0.000001 # beware, despite being L=1 is still a vector

# Testing interface
res = fit_model_univariate(f, x, y, -10.0, 100.0
                           , fitmethod = Optim.Brent()
                           ).param_results
@test norm(res-4.2) < 1e-5 

res = fit_model_univariate(f, x, y, -10.0, 100.0
                           , fitmethod = Optim.GoldenSection()
                           ).param_results
@test norm(res-4.2) <  1e-5

# Fit a 2D function, handmade example
x = 0:0.1:1 |> collect
y = x.^2
f(x,p) = x^p[1] + p[2]
res = fit_model_multivariate(f, x, y, ones(2)).param_results
@test norm(res-[2.0,0.0]) < 1e-5

# Fit a 3D function, handmade example
x = 0:0.1:1 |> collect
y = x.^1.7 + π
f(x,p) = p[1]*x^p[2] + p[3]
res = fit_model_multivariate(f, x, y, ones(3)).param_results
@test norm(res-[1.0,1.7,π]) < 1e-4


info("All test correct ✓")
