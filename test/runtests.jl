using CamiXon
using Test

@testset "CamiXon.jl" begin
    @test find_indices(collect("ahsgh")) == [[1], [2, 5], [3], [4]]
    @test find_indices(collect("ahsgh"),'h') == [[2, 5]]
    @test find_indices([1,2,3,4,2]) == [[1], [2, 5], [3], [4]]
    @test find_indices([1,2,3,4,2],2) == [[2, 5]]
end
