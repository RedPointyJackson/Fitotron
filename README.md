Fitotron
===========

Fitotron provides ordinary least squares with a simple interface. It
uses a Nelder-Mead method to find the optimum parameters of the fit.

It can (optionaly) rescale the parameter uncertainties using the minimum value of the sum of residuals function.

Usage
-----------

We load the necesary libs:

```jl
using Fitotron
using DataFrames
using Gadfly # To see the results
```
We read the data from a plain text file and create the fit function:

```jl
# Read the data
df = readtable("data.dat",header=false,separator='\t')

# Create the fit model
fit_fun(x,params) = params[1] + params[2]*cos(params[3]*x-params[4])
```

We call the fitting function:

```jl
# Fit it 
fit_result = fit_model(fit_fun,df,[1.0,1.0,0.1,1.0])
```
An optional argument `rescale`, that can be true or false, makes the minimum of
the sum of the residues squared be the number of degrees of freedom (minus the
number of parameters) at the minimum. In other words; if enabled, the parameter
standard deviations will be multiplied by the root of S/(dof), where dof is the
number of datapoints minus the number of parameters and S is the residue of the
cost function minimization.

The returned `fit_result` contains in its fields all the available data:


* `param_results`  Fit results
* `param_stdevs`  Standard deviations at 1 σ
* `covariance`  Covariance of the parameters
* `resid` final sum of residuals
* `fit_value` gives value of fit function in x
* `fit_stdev` gives stderr of fit function in x

We can create a plot to see the results, with the help of the `fit_value` and `fit_stdev` functions:

```jl
# See the result. Create a dataframe with the fit results:
xx = linspace(minimum(df[1]),maximum(df[1]),100)
yy = [fit_result.fit_value(x)::Float64 for x in xx]
ss = [fit_result.fit_stdev(x)::Float64 for x in xx]


df_fit = DataFrame( x = xx,
		    y = yy,
		    ymin = yy-ss,
		    ymax = yy+ss)

# We also need the data in a dataframe:
df_data = DataFrame( x = df[1],
		    y = df[2],
		    ymin = df[2]-df[3],
		    ymax = df[2]+df[3])

plot(
     layer(df_data,x=:x,y=:y,ymin=:ymin,ymax=:ymax,Geom.errorbar),
     layer(df_fit,x=:x,y=:y,ymin=:ymin,ymax=:ymax,Geom.ribbon))
```

![Fitted!](http://i.imgur.com/mp9XHYw.png)
