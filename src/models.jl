"""
Most general representation of a fitting model.
"""
abstract AbstractModel

"""

    CustomModel(func, dims, x, y, yerr, rescale)

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
    LinearModel(x,y,yerr,rescale) = new( convert(Vector{Float64},1.0*x)
                                        ,convert(Vector{Float64},1.0*y)
                                        ,convert(Vector{Float64},1.0*yerr)
                                        ,rescale)
    # Convenience methods
    LinearModel(y)                = new(1:length(y), y ,ones(y), true)
    LinearModel(x,y)              = new(x, y, ones(x), true)
    LinearModel(x,y,yerr)         = new(x, y, yerr, false)
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
