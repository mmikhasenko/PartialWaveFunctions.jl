using PartialWaveFunctions
using Test

@testset "PartialWaveFunctions.jl" begin
    # Clebsch Gordan coefficient
    @test clebschgordan_doublearg(3,  3, 1,  1, 2*2, 2*2) == 1
    @test clebschgordan_doublearg(3,  3, 1, -1, 2*1, 2*1) ≈ sqrt(3)/2
    @test clebschgordan_doublearg(3, -1, 1,  1, 2*1, 2*0) ≈ -sqrt(2)/2

    # Wigner Function
    @test wignerD(3, 2, 1, π, 0, -π) ≈ -sqrt(10)/8 + 0im
    @test wignerd(1, 1, 0, 0.3) ≈ -sqrt(1-0.3^2)/sqrt(2)
end
