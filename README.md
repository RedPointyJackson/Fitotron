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
x = linspace(0,1,N) |> collect
y = x + 0.05randn(N)

fun(x,p) = p[1] + p[2]*x
fit = fitmodel(fun,x,y,[1,1])

p = plotfit(fit)
```

`p` will be a `Gadfly.Plot` like the following:
![fit result](https://github.com/RedPointyJackson/Fitotron/blob/master/fitresult.png)


And a `fit` a `FitResult` that shows like
```
Fit results:
─────────────────────────────────────────────────────────
Param. 1:                       -1.40e-02 ± 1.55e-02
Param. 2:                       1.02e+00 ± 3.76e-02
Reduced χ²:                     0.0030915027388933136
Parameter estimation method:    Nelder-Mead
Uncertainty estimation method:  χ² sweeping (rescaled)
```
in the REPL.

## Caveats
The uncertainty estimation is very different to `gnuplot`'s one, for
example, and a lot of times very pessimistic.
