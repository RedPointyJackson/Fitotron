# Fitotron

Fitotron provides ordinary least squares with a simple interface. It
uses a Nelder-Mead method (or a Brent method, in univariate fits) to
find the optimum parameters of the fit.

It can rescale the parameter uncertainties using the
minimum value of the sum of residuals function; this is the default
behaivour when y errors are not provided, although it can be disabled.

## Usage
In a REPL, type:

```julia
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
───────────────────────────────────────────────────────────────
Param. 1:                       +7.57e-01 ± 5.31e-01 (70.1%)
Param. 2:                       -2.49e-01 ± 1.77e-01 (71.0%)
Reduced χ²:                     0.4950205376301822
Parameter estimation method:    Nelder-Mead
Uncertainty estimation method:  χ² sweeping (rescaled)
```
in the REPL. A fit to a sine would be better, indeed:

```julia
sine(x,p) = p[1] + p[2]*sin(p[3]*x+p[4])
fit = fitmodel(sine,x,y,ones(4))
plotfit(fit)
```

Which gives:
```
Fit results:
───────────────────────────────────────────────────────────────
Param. 1:                       -3.04e-02 ± 1.89e-02 (62.1%)
Param. 2:                       +9.16e-01 ± 5.09e-01 (55.6%)
Param. 3:                       +1.03e+00 ± 6.97e-01 (67.9%)
Param. 4:                       +7.75e-02 ± 4.11e-02 (53.0%)
Reduced χ²:                     0.30947383110883564
Parameter estimation method:    Nelder-Mead
Uncertainty estimation method:  χ² sweeping (rescaled)
```
and

![sine fit result](https://github.com/RedPointyJackson/Fitotron/blob/master/fitresult_sine.png)

The Jacobian method of error estimation (`uncmethod=:jacobian` in
`fitmodel`) gives errors that are a lot
less conservative:
```
Fit results:
───────────────────────────────────────────────────────────────
Param. 1:                       -3.04e-02 ± 7.92e-02 (260.8%)
Param. 2:                       +9.16e-01 ± 1.16e-01 (12.7%)
Param. 3:                       +1.03e+00 ± 6.06e-02 (5.9%)
Param. 4:                       +7.75e-02 ± 2.31e-01 (298.2%)
Reduced χ²:                     0.30947383110883564
Parameter estimation method:    Nelder-Mead
Uncertainty estimation method:  Jacobian (rescaled)
```
and

![sine fit result (jacobian)](https://github.com/RedPointyJackson/Fitotron/blob/master/fitresult_sine_jac.png)

There are two posible invocations, one for univariate fits and other for multivariate
ones. Use `?fitmodel` to see the documentation.

The function used to fit should be in the form `f(x,p)` where `p` is the vector containing the parameters.

## Caveats
The chi sweep uncertainty estimation is very different to `gnuplot`'s one, for
example, and a lot of times very pessimistic.
