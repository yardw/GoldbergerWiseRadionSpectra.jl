using GoldbergerWiseRadionSpectra
using Test

@testset "GoldbergerWiseRadionSpectra.jl" begin
    @test m2s_analytical(1e-3, 1e3) isa Number
    @test isapprox(search(1e-3, 1e3)[1], m2s_analytical(1e-3, 1e3), rtol=1e-1)
end
