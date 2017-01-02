# Fitotron

[![Build Status](https://travis-ci.org/RedPointyJackson/Fitotron.jl.svg?branch=master)](https://travis-ci.org/RedPointyJackson/Fitotron.jl)

Fitotron provides ordinary least squares for 2D data with a simple
interface. It uses by default a Nelder-Mead method (or a Brent method,
in univariate fits) to find the optimum parameters of the fit, thanks
to the `Optim` library.

It can rescale the parameter uncertainties using the minimum value of
the sum of residuals function; this is the default behaviour when y
errors are not provided, although it can be disabled.

## Usage
In a REPL, type:

```julia
using Fitotron

N = 50
x = linspace(0,2π,N)
y = sin(x) + 0.5randn(N)

fun(x,p) = p[1] + p[2]*x
model = CustomModel(fun,2,x,y)
fit = fitmodel(model)

p = plotfit(fit)
c = plotcost(fit)
```

`p,c` will be two `Gadfly.Plot` like the following:

![fit result](https://github.com/RedPointyJackson/Fitotron.jl/blob/master/fitresult.png)

![fit cost](https://github.com/RedPointyJackson/Fitotron.jl/blob/master/fitcost.png)

And a `fit` a `FitResult` that shows like
```
Fit results:
───────────────────────────────────────────────────────────────
Param. 1:                       +1.14e+00 ± 1.86e-01 (16.3%)
Param. 2:                       -3.49e-01 ± 5.11e-02 (14.6%)
Reduced χ²                      0.4467
Errors were rescaled so χ²=d.o.f.
```
in the REPL. A fit to a sine would be better, indeed:

```julia
sine(x,p) = p[1] + p[2]*sin(p[3]*x+p[4])
model = CustomModel(sine,4,x,y)
fit = fitmodel(model)

p = plotfit(fit)
```

Which gives:
```
Fit results:
───────────────────────────────────────────────────────────────
Param. 1:                       +4.48e-02 ± 7.30e-02 (163.0%)
Param. 2:                       +1.07e+00 ± 1.04e-01 (9.8%)
Param. 3:                       +9.27e-01 ± 5.42e-02 (5.9%)
Param. 4:                       +2.11e-01 ± 1.96e-01 (92.8%)
Reduced χ²                      0.2645
Errors were rescaled so χ²=d.o.f.
```
and

![sine fit result](https://github.com/RedPointyJackson/Fitotron.jl/blob/master/fitresult_sine.png)

# TODO/Bugs:
## Cost function plot
- It has two legends
- Works only for functions of 2 parameters, it
  should let the user choose which dimensions should it choose for
  models of higher dimensionality and be a 1D plot for 1D models.
## Linear models
- Uncertainties are way smaller than the `CustomModel` ones. Is that
  the expected behaivour?
