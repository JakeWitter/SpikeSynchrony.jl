using Test, SpikeSynchrony

@testset "SPIKE distance" begin
    include("SPIKEdistance_tests.jl")
end
@testset "vanRossum" begin
    include("vanRossum_tests.jl")
end
@testset "VictorPurpura" begin
    include("VictorPurpura_tests.jl")
end
