using Base.Test
using Fitotron

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
end

# @testset "Param correctness" begin
#     ############################
#     #                          #
#     #             1D           #
#     #                          #
#     ############################
#     # Fit a 1D function, handmade example
#     x = linspace(0.01,1,10) |> collect
#     # We don't use 0 because finite differences
#     # will use -0.001 and f isn't defined there
#     y = x.^4.2
#     f(x,p) = x^p[1]

#     # With fitmodel_raw
#     ft = Fitotron.fitmodel_raw(
#                                f
#                                , x
#                                , y
#                                , [1.0, 100.0]
#                                , nothing
#                                , ones(x)
#                                , Optim.Brent()
#                                , :chisweep
#                                , false
#                                , :univariate
#                                )
#     res = ft.param_results
#     @test norm(res-4.2) < 1e-6
#     # Notice the norm, despite being L=1 is still a vector

#     # Do it also with the user interface
#     ft = Fitotron.fitmodel(
#                            f
#                            , x
#                            , y
#                            , -10
#                            , 100
#                            )
#     res = ft.param_results
#     @test norm(res-4.2) < 1e-6

#     ############################
#     #                          #
#     #             2D           #
#     #                          #
#     ############################
#     # Fit a 2D function, handmade example
#     x = linspace(1,100,100) |> collect
#     y = x.^2
#     g(x,p) = x^p[1] + p[2]
#     # With fitmodel_raw
#     ft = Fitotron.fitmodel_raw(
#                                g
#                                , x
#                                , y
#                                , [1.0, 1.0]
#                                , nothing
#                                , ones(x)
#                                , Optim.NelderMead()
#                                , :chisweep
#                                , false
#                                , :multivariate
#                                )
#     res = ft.param_results
#     @test norm(res-[2.0,0.0]) < 1e-5
#     # # Do it also with the user interface
#     ft = Fitotron.fitmodel(
#                            g
#                            , x
#                            , y
#                            , [1, 1]
#                            )
#     res = ft.param_results
#     @test norm(res-[2.0,0.0]) < 1e-5


#     ############################
#     #                          #
#     #             3D           #
#     #                          #
#     ############################
#     # # Fit a 3D function, handmade example
#     x = linspace(0.01,1,1000) |> collect
#     y = x.^1.7 + π
#     h(x,p) = p[1]*x^p[2] + p[3]
#     # With fitmodel_raw
#     ft = Fitotron.fitmodel_raw(
#                                h
#                                , x
#                                , y
#                                , [1.0, 1.0, 1.0]
#                                , nothing
#                                , ones(x)
#                                , Optim.NelderMead()
#                                , :chisweep
#                                , false
#                                , :multivariate
#                                )
#     res = ft.param_results
#     @test norm(res-[1.0, 1.7, π]) < 1e-3
#     # # Do it also with the user interface
#     ft = Fitotron.fitmodel(
#                            h
#                            , x
#                            , y
#                            , [1, 1, 1]
#                            )
#     res = ft.param_results
#     @test norm(res-[1.0, 1.7, π]) < 1e-5

# end
# @testset "Uncertainty" begin
#     ############################
#     #                          #
#     #             1D           #
#     #                          #
#     ############################
#     x = linspace(0.01,1,10) |> collect
#     y = x.^3
#     f(x,p) = x^p[1]
#     ftchi = Fitotron.fitmodel(
#                               f ,x ,y, -10 ,100
#                               , uncmethod = :chisweep
#                               )
#     ftjac = Fitotron.fitmodel(
#                               f ,x ,y, -10 ,100
#                               , uncmethod = :jacobian
#                               )
#     # Results are reescaled (no yerr) so the
#     # #error must be ridiculous.
#     devchi = ftchi.param_deviations
#     devjac = ftjac.param_deviations
#     @test norm(devchi) < 1e-8
#     @test norm(devjac) < 1e-8
#     ############################
#     #                          #
#     #             2D           #
#     #                          #
#     ############################
#     x = linspace(0.01,1,1000) |> collect
#     y = x.^3
#     g(x,p) = p[1]*abs(x)^p[2]
#     # χ sweep
#     ftchi = Fitotron.fitmodel(
#                               g ,x ,y, ones(2)
#                               , uncmethod = :chisweep
#                               )
#     ftjac = Fitotron.fitmodel(
#                               g ,x ,y, ones(2)
#                               , uncmethod = :jacobian
#                               )
#     devchi = ftchi.param_deviations
#     devjac = ftjac.param_deviations
#     @test norm(devchi) < 1e-5
#     @test norm(devjac) < 1e-5
# end
# @testset "Bugs" begin
#     x    = randn(50)
#     y    = 2x + randn(50)
#     yerr = randn(50)
#     # yerr are ignored (last seen in 12b1a9b).
#     f(x,p) = p[1]*x+p[2]
#     g(x,p) = p[1]*x
#     data_multi = Fitotron.fitmodel(f,x,y,[1,1]
#                                    ,yerr=yerr).data
#     data_mono = Fitotron.fitmodel(g,x,y,0,3
#                                   ,yerr=yerr).data
#     @test data_multi[:,3] == yerr
#     @test data_mono[:,3] == yerr
#     # using non-vectorized functions for univariate fit (last seen
#     # in 12b1a9b).
#     h(x,α) = α*x
#     ft = Fitotron.fitmodel(h,x,y,0,3)
#     @test isa(ft,Fitotron.FitResult)
# end

# @testset "Plotting tests" begin
#     const N = 50

#     x1 = linspace(0,1,N) |> collect
#     y1 = x1.^3 + 0.1randn(N)
#     f1(x,p) = p[1]*x^p[2]
#     fit1 = fitmodel(f1,x1,y1,[1,1])
#     @test isa(plotfit(fit1),Gadfly.Plot)

#     x2 = linspace(0,1,N) |> collect
#     y2 = x2 + 0.1randn(N) + 2
#     f2(x,p) = p[1]*x + p[2]
#     fit2 = fitmodel(f2,x2,y2,[1,1])
#     @test isa(plotfit(fit2),Gadfly.Plot)
# end

# @testset "Cost function" begin

#     @testset "Linefit cost" begin
#         const N = 50
#         x = linspace(0,1,N) |> collect
#         y = 2x + 1
#         expected_cost(m,n)=sumabs2(y[i]-(m*x[i]+n) for i ∈ 1:50)
#         line(x,p) = p[1]*x+p[2]
#         ft = Fitotron.fitmodel(line,x,y,[1,1])
#         computed_cost = ft.cost
#         for m in randn(10)
#             for n in randn(10)
#                 @test computed_cost([m,n]) ≈ expected_cost(m,n)
#             end
#         end
#     end

#     @testset "Constantfit cost" begin
#         const N = 50
#         # Fit yᵢ=1 to α. The cost function
#         # should be ∑(1-α)² = N×(1-α)
#         x = linspace(0,1,N)
#         y = ones(N)
#         expected_cost(α) = N*(1-α[1])^2
#         ft = Fitotron.fitmodel((x,p)->p[1],x,y,0,2)
#         computed_cost = ft.cost
#         for α in randn(100)
#             @test computed_cost([α]) ≈ expected_cost(α)
#         end
#     end

#     @testset "Known χ² dev." begin
#         # For y=[1,2], a fit to f(x)=α gives
#         # cost(α) = 2α² - 6α + 5
#         # In the minimum (k), it gives f(1.5)=0.5.
#         # The variance is 2k, so the desvest is √(2k).
#         # That is, σ = √(1) = 1
#         # The desvest is the α with cost(α) = 1, which
#         # is α = {1,2}. So, from 1.5 to 1 or
#         # 2 there are 0.5 units,
#         # which is the thing we are going to check:
#         x = [1.0; 2.0]
#         y = [1.0; 2.0]
#         fit_fun(x,p) = p[1]
#         ft = Fitotron.fitmodel(fit_fun,x,y
#                                ,0,4,rescaling=false)
#         @test ft.param_deviations ≈ [0.5]
#     end

# end
