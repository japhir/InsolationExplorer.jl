using InsolationExplorer
using Test

@testset "InsolationExplorer.jl" begin
    @test InsolationExplorer.insolation(0.016705, 0.4090928042223287, 1.7962486166737615;
                                        longitude = pi/2, latitude = deg2rad(65),
                                        S0 = 1360.7, H = nothing) == 509.9881957773474
end
