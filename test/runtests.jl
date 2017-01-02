using Base.Test
using Fitotron
using Gadfly

srand(42)

@testset "Error propagation" begin
    # Error propagation, 1D, handmade example
    f(x,p) = x*p
    @test Fitotron.propagate_errors(f,1.0,π,0) ≈ 0.0
    g(x,p) = x^p
    @test Fitotron.propagate_errors(g,1.0,2.0,3.0) ≈ 6log(3.0)

    # Error propagation, 2D, handmade example
    h(x,p) = x*p[1]+p[2]
    @test Fitotron.propagate_errors(h,ones(2),ones(2),0) ≈ 1.0
    @test Fitotron.propagate_errors(h,ones(2),ones(2),1) ≈ √2

    # Error propagation, 3D, handmade example
    i(x,p) = x*p[1]+p[2]*p[3]
    result(x,p,δp) = sqrt( (x*δp[1])^2 + (p[3]*δp[2])^2 + (p[2]*δp[3])^2 )
    @test Fitotron.propagate_errors(i,ones(3),ones(3),0) ≈ result(0,ones(3),ones(3))
    @test Fitotron.propagate_errors(i,Float64[1,2,3],0.5*ones(3),0.2) ≈ result(0.2,Float64[1,2,3],0.5*ones(3))
end

@testset "Custom model creation" begin
    f(x,p) = p[1] + x*p[2]
    x = rand(20)
    y = rand(20)
    yerr = rand(20)

    # Basic model, no options specified
    @test CustomModel(f, 2, x, y).x       == x
    @test CustomModel(f, 2, x, y).y       == y
    @test CustomModel(f, 2, x, y).func    == f
    @test CustomModel(f, 2, x, y).yerr    == ones(x)
    @test CustomModel(f, 2, x, y).rescale == true
    # Test bounds length as a function of dims
    # Note that dims must be the length of p in f(x,p) to do a correct
    # fit.
    @test length(CustomModel(f, 1, x, y).bounds) == 2
    @test length(CustomModel(f, 2, x, y).bounds) == 2
    @test length(CustomModel(f, 3, x, y).bounds) == 3

    # User given yerr
    @test CustomModel(f, 2, x, y, yerr=yerr).yerr    == yerr
    @test CustomModel(f, 2, x, y, yerr=yerr).rescale == false
end

@testset "Linear model creation" begin
    x = rand(20)
    y = rand(20)
    yerr = rand(20)

    @test LinearModel(x,y,yerr,true).x        == x
    @test LinearModel(x,y,yerr,true).y        == y
    @test LinearModel(x,y,yerr,true).yerr     == yerr
    @test LinearModel(x,y,yerr,true).rescale  == true
    @test LinearModel(x,y,yerr,false).rescale == false

    @test LinearModel(y).x             == Fitotron.thing2array(1:length(y))
    @test LinearModel(y).y             == y
    @test LinearModel(y).yerr          == ones(y)
    @test LinearModel(y).rescale       == true
    @test LinearModel(y,true).rescale  == true
    @test LinearModel(y,false).rescale == false

    @test LinearModel(x,y).x             == x
    @test LinearModel(x,y).y             == y
    @test LinearModel(x,y).yerr          == ones(y)
    @test LinearModel(x,y).rescale       == true
    @test LinearModel(x,y,true).rescale  == true
    @test LinearModel(x,y,false).rescale == false

    @test LinearModel(x,y,yerr).x             == x
    @test LinearModel(x,y,yerr).y             == y
    @test LinearModel(x,y,yerr).yerr          == yerr
    @test LinearModel(x,y,yerr).rescale       == false
    @test LinearModel(x,y,yerr,true).rescale  == true
    @test LinearModel(x,y,yerr,false).rescale == false

end

@testset "Analytical line fit" begin
    # Compare to the same algorithm, implemented in another software.
    #    Trivial
    x    = linspace(0,10,20) |> collect
    y    = linspace(0,10,20) |> collect
    yerr = ones(20)
    model = LinearModel(x,y,yerr)
    fit = fitmodel(model)
    @test fit.cost(fit.param_results) ≈ 0
    m ,  n = fit.param_results
    σm, σn = fit.param_deviations
    @test  m ≈ 1
    @test σm ≈ 0.07367883976130075
    @test  n ≈ 0
    @test σn ≈ 0.43094580368566743

    #    Not trivial
    x     = [1,2,3,4,5,6,7,8,9]*1.0
    y     = [1,2,3,8,5,6,7,8,9]*1.0
    yerr  = [0.1,0.2,0.3,0.8,0.5,0.6,0.7,0.8,0.9]
    model = LinearModel(x,y,yerr)
    fit   = fitmodel(model)
    m ,  n = fit.param_results
    σm, σn = fit.param_deviations
    @test  m ≈ 1.03898314436689
    @test σm ≈ 0.05287908290229
    @test  n ≈ 0.24952872260996
    @test σn ≈ 0.12430728667796

    # Try the plotting functions
    plotfit(fit)
    plotcost(fit)
end

@testset "Custom line fit" begin
    x    = linspace(0,2π,20)
    y    = sin(x) + randn(20)
    model = CustomModel(
                        (x,p)->p[1]*sin(p[2]*x), 2
                        ,x,y
                        )
    fit = fitmodel(model)

    # Try the plotting functions
    plotfit(fit)
    plotcost(fit)

    # Check if for a line it gives the same results as an analytical
    # fit.
    x = linspace(0,1,100)
    y = 2x+1+randn(100)
    linmodel = LinearModel(x,y)
    cusmodel = CustomModel((x,p)->p[1]*x+p[2],2,x,y)
    linparams = fitmodel(linmodel).param_results
    cusparams = fitmodel(cusmodel).param_results
    @test norm(linparams-cusparams)/norm(linparams) < 1e4
end
