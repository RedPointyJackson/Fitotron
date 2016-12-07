"""
Most general representation of a fitting model.
"""
abstract AbstractModel



defaultOptim(dims) = dims == 1 ? Optim.Brent() : Optim.NelderMead()
"""

    CustomModel(func, dims, x, y, yerr, rescale, optimizator)

General fitting model, with user provided custom function. If
`rescale`, force χ² = d.o.f. in its minimum.
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

    # User input may (will) be messy
    CustomModel(func,dims
                ,x,y,yerr
                ,rescale,optimizator
                ,bounds
                ) = new(
                        func
                        ,Int64(dims)
                        ,thing2array(x)
                        ,thing2array(y)
                        ,thing2array(yerr)
                        ,rescale
                        ,optimizator
                        ,thing2array(bounds)
                        )
    # Convenience methods
    CustomModel(func,x,y,dims::Int64
                ) = CustomModel(
                                func, dims
                                , x, y, ones(y)
                                , true, defaultOptim(dims)
                                , dims==1? [-1.8;1e8] : ones(dims)
                                )
end

function show(io::IO, m::CustomModel)
    println(io, "Custom model")
    println(io, "   $(m.dims) parameters")
    println(io, "   rescale: $(m.rescale)")
    println(io, "   function: $(m.func)")
    # Print the data
end

"""
    LinearModel(x, y, yerr, rescale)

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
    LinearModel(x,y,yerr,rescale) = new( thing2array(x)
                                        ,thing2array(y)
                                        ,thing2array(yerr)
                                        ,rescale)
    # Convenience methods
    LinearModel(y)                = LinearModel(1:length(y), y ,ones(y), true)
    LinearModel(x,y)              = LinearModel(x, y, ones(x), true)
    LinearModel(x,y,yerr)         = LinearModel(x, y, yerr, false)
end

function show(io::IO, m::LinearModel)
    if m.rescale
        println(io, "Linear model with rescaling")
    else
        println(io, "Linear model without rescaling")
    end
end


"""

    QuadraticModel(x, y, yerr, rescale)

Quadratic fitting model, of function y = ax²+bx+c.
"""
immutable QuadraticModel <: AbstractModel
    x       ::Vector{Float64}
    y       ::Vector{Float64}
    yerr    ::Vector{Float64}
    rescale ::Bool
end
