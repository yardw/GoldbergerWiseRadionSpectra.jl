using GoldbergerWiseRadionSpectra
using Test

@testset "GoldbergerWiseRadionSpectra.jl" begin
    l2 = 1e-3
    g2 = 1e3
    @test m2s_analytical(l2, g2) isa Number
    @test isapprox(search(l2, g2)[1], m2s_analytical(l2, g2), rtol=1e-1)
    @test extract_sol([l2, g2, search(l2, g2)...]) isa Function
end
