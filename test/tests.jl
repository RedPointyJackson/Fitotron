using Base.Test
using Fitotron


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
    # trivial
    x    = linspace(0,10,20)
    y    = linspace(0,10,20)
    yerr = 0.02ones(20)
    m,n,C = Fitotron.fit_line(x,y,yerr)
    analyticC = 
        [2.171428571428572e-6 -1.0857142857142861e-5
         -1.0857142857142861e-5 7.428571428571431e-5]
    @test m ≈ 1
    @test n ≈ 0
    @test C ≈ analyticC

    # not trivial
    x    = [1,2,3,4,5,6,7,8,9]
    y    = [1,2,3,8,5,6,7,8,9]
    yerr = [0.1,0.2,0.3,0.8,0.5,.06,0.7,0.8,0.9]
    m,n,C = Fitotron.fit_line(x,y,yerr)
    analyticC = 
        [ 0.0004785877394248638 -0.0021591099755238854
          -0.0021591099755238854 0.012097533052456966 ]
    @test m ≈ 0.998470256138598
    @test n ≈ 0.4520931637514547
    @test C ≈ analyticC
end

@testset "Fitting tests" begin
    @testset "Param correctness" begin
        ############################
        #                          #
        #             1D           #
        #                          #
        ############################
        # Fit a 1D function, handmade example
        x = linspace(0.01,1,10) |> collect 
        # We don't use 0 because finite differences
        # will use -0.001 and f isn't defined there
        y = x.^4.2
        f(x,p) = x^p[1]

        # With fitmodel_raw
        ft = Fitotron.fitmodel_raw(
                                    f
                                    , x
                                    , y
                                    , [1.0, 100.0]
                                    , nothing
                                    , ones(x)
                                    , Optim.Brent()
                                    , :chisweep
                                    , false
                                    , :univariate
                                    )
        res = ft.param_results
        @test norm(res-4.2) < 1e-6
        # Notice the norm, despite being L=1 is still a vector

        # Do it also with the user interface
        ft = Fitotron.fitmodel(
                                f
                                , x
                                , y
                                , -10
                                , 100
                                )
        res = ft.param_results
        @test norm(res-4.2) < 1e-6

        ############################
        #                          #
        #             2D           #
        #                          #
        ############################
        # Fit a 2D function, handmade example
        x = linspace(1,100,100) |> collect
        y = x.^2
        g(x,p) = x^p[1] + p[2]
        # With fitmodel_raw
        ft = Fitotron.fitmodel_raw(
                                    g
                                    , x
                                    , y
                                    , [1.0, 1.0]
                                    , nothing
                                    , ones(x)
                                    , Optim.NelderMead()
                                    , :chisweep
                                    , false
                                    , :multivariate
                                    )
        res = ft.param_results
        @test norm(res-[2.0,0.0]) < 1e-5
        # # Do it also with the user interface
        ft = Fitotron.fitmodel(
                                g
                                , x
                                , y
                                , [1, 1]
                                )
        res = ft.param_results
        @test norm(res-[2.0,0.0]) < 1e-5


        ############################
        #                          #
        #             3D           #
        #                          #
        ############################
        # # Fit a 3D function, handmade example
        x = linspace(0.01,1,1000) |> collect
        y = x.^1.7 + π
        h(x,p) = p[1]*x^p[2] + p[3]
        # With fitmodel_raw
        ft = Fitotron.fitmodel_raw(
                                    h
                                    , x
                                    , y
                                    , [1.0, 1.0, 1.0]
                                    , nothing
                                    , ones(x)
                                    , Optim.NelderMead()
                                    , :chisweep
                                    , false
                                    , :multivariate
                                    )
        res = ft.param_results
        @test norm(res-[1.0, 1.7, π]) < 1e-3
        # # Do it also with the user interface
        ft = Fitotron.fitmodel(
                                h
                                , x
                                , y
                                , [1, 1, 1]
                                )
        res = ft.param_results
        @test norm(res-[1.0, 1.7, π]) < 1e-5

    end

    @testset "Uncertainty" begin
        ############################
        #                          #
        #             1D           #
        #                          #
        ############################
        x = linspace(0.01,1,10) |> collect 
        y = x.^3
        f(x,p) = x^p[1]
        ft = Fitotron.fitmodel( 
                                 f ,x ,y, -10 ,100
                                 , uncmethod = :chisweep
                                )
        # Results are reescaled (no yerr) so the
        # #error must be ridiculous.
        dev = ft.param_deviations
        @test norm(dev) < 1e-8
        ############################
        #                          #
        #             2D           #
        #                          #
        ############################
        x = linspace(0.01,1,1000) |> collect 
        y = x.^3
        g(x,p) = p[1]*abs(x)^p[2]
        # χ sweep
        ft = Fitotron.fitmodel( 
                                 g ,x ,y, ones(2)
                                 , uncmethod = :chisweep
                                )
        dev = ft.param_deviations
        @test norm(dev) < 1e-5
    end
end
