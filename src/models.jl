"""
Most general representation of a fitting model.
"""
abstract AbstractModel

"""
    defaultOptim(dims)
Gives the default optimizer to use in fits of `dims` parameters.
"""
defaultOptim(dims) = dims == 1 ? Optim.Brent() : Optim.NelderMead()

"""

    CustomModel(func, dims, x, y, yerr, rescale, optimizator)

    CustomModel(func, dims, x, y [; yerr        = nothing
                                  , rescale     = true unless yerr given
                                  , optimizator = Brent/NelderMead
                                  , bounds      = [-1e8,1e8] / [1,1,…]
                                 ]
               )

General fitting model, with user provided custom function of `dims`
parameters. If `rescale`, force χ² = d.o.f. in its minimum by
rescaling `yerr` accordingly.

# Parameters
* `func`: Function to fit the model, of the form `f(x,p)` where `p` is
  the vector of parameters. `length(p)` must be `dims`.
* `x`,`y`: Vectors of `x` and `y` data.
* `yerr`: Vector of y errors. If omited, will default to a vector of
  ones and rescaling will be enabled ,unless specified otherwise with
  the `rescale` parameter.
* `rescale`: If enabled, χ² will be the d.o.f. in the minimum by
  accordingly rescaling `yerr`. By default is enabled unless the user
  specifies `yerr`.
* `optimizator`: `Optim.Optimizator` to use. Defaults to `Optim.Brent`
  for `dims==1` and `Optim.NelderMead` for `dims>1`.
* `bounds`: For 1 parameter fits, bounds of the parameter (defaults to
  -1e8 to 1e8), for fits with many parameters is a vector with initial
  values (defaults to a vector of ones).
"""
immutable CustomModel <: AbstractModel
    func    ::Function
    dims    ::Int64
    x       ::Vector{Float64}
    y       ::Vector{Float64}
    yerr    ::Vector{Float64}
    rescale ::Bool
    optimizator
    bounds :: Vector{Float64}

    function CustomModel(func,dims
                ,x,y,yerr
                ,rescale,optimizator
                ,bounds
                )
        # User input may (will) be messy
        if dims < 1
            error("Models must have dims>1.")
        end
        if !(length(x) == length(y) == length(yerr))
            error("length of x,y,yerr must be the same.")
        end
        if dims == 1
            if length(bounds) != 2
                error("For dims=1, bounds must be a length 2 array.")
            end
        else
            if length(bounds) != dims
                error("For dims != 1, length(bounds) must be dims.")
            end
        end
        return new(
                   func
                   ,Int64(dims)
                   ,thing2array(x)
                   ,thing2array(y)
                   ,thing2array(yerr)
                   ,rescale
                   ,optimizator
                   ,thing2array(bounds)
                   )
    end
    # Convenience method
    CustomModel(func,dims::Int64,x,y
                ; yerr        = nothing
                , rescale     = yerr==nothing
                , optimizator = defaultOptim(dims)
                , bounds      = dims==1? [-1e8;1e8] : ones(dims)
                ) = CustomModel(
                                  func
                                , dims
                                , x, y
                                , yerr == nothing? ones(x) : yerr
                                , rescale, optimizator
                                , bounds
                                )
end

function show(io::IO, m::CustomModel)
    println(io, "Custom model of $(m.dims) parameter$(m.dims==1?' ':'s')")
    println(io, "   Rescale: $(m.rescale)")
    println(io, "   Function: $(m.func)")
end

"""
    LinearModel(x,y,yerr,rescale)
    LinearModel(y)
    LinearModel(x,y)
    LinearModel(x,y,yerr)

Linear fitting model, of function y = mx+n.
If omitted,

* `x` will be 1,2,…
* `yerr` will be a vector of ones
* `rescale` will be true except if `yerr` is provided.
"""
immutable LinearModel <: AbstractModel
    x       ::Vector{Float64}
    y       ::Vector{Float64}
    yerr    ::Vector{Float64}
    rescale ::Bool

    # User input may (will) be messy
    function LinearModel(x,y,yerr,rescale::Bool)
        if !(length(x) == length(y) == length(yerr))
            error("length of x,y,yerr must be the same.")
        end
        return new( thing2array(x)
                    ,thing2array(y)
                    ,thing2array(yerr)
                    ,rescale)
    end
    # Convenience methods
    LinearModel(y)                 = LinearModel(1:length(y),y,ones(y),true)
    LinearModel(y,rescale::Bool)   = LinearModel(1:length(y),y,ones(y),rescale)

    LinearModel(x,y)               = LinearModel(x,y,ones(x),true)
    LinearModel(x,y,rescale::Bool) = LinearModel(x,y,ones(x),rescale)

    LinearModel(x,y,yerr)          = LinearModel(x,y,yerr,false)
end

function show(io::IO, m::LinearModel)
    if m.rescale
        println(io, "Linear model with rescaling")
    else
        println(io, "Linear model without rescaling")
    end
end
