# Fitotron

Fitotron provides ordinary least squares with a simple interface. It
uses a Nelder-Mead method (or a Brent method, in univariate fits) to
find the optimum parameters of the fit.

It can rescale the parameter uncertainties using the
minimum value of the sum of residuals function; this is the default
behaivour when y errors are not provided, although it can be disabled.

## Usage
In a REPL, type:

```jl
using Fitotron

srand(42)
const N = 50
x = linspace(0,2π,N) |> collect
y = sin(x) + 0.5randn(N)

fun(x,p) = p[1] +p[2]*x
fit = fitmodel(fun,x,y,[1,1])

p = plotfit(fit)
```

`p` will be a `Gadfly.Plot` like the following:

![fit result](https://github.com/RedPointyJackson/Fitotron/blob/master/fitresult.png)


And a `fit` a `FitResult` that shows like
```
Fit results:
─────────────────────────────────────────────────────────
Param. 1:                       7.57e-01 ± 5.31e-01
Param. 2:                       -2.49e-01 ± 1.77e-01
Reduced χ²:                     0.4950205376301822
Parameter estimation method:    Nelder-Mead
Uncertainty estimation method:  χ² sweeping (rescaled)
```
in the REPL. A fit to a sine would be better, indeed.

There are two posible invocations, for univariate and multivariate
fits. Use `?fitmodel` to see the documentation.

The function used to fit should be in the form `f(x,p)` where `p` is the vector containing the parameters.

## Caveats
The uncertainty estimation is very different to `gnuplot`'s one, for
example, and a lot of times very pessimistic.
