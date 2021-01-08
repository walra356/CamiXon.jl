using CamiXon
using Test

@testset "CamiXon.jl" begin
    @test indices(collect("ahsgh")) == [[1], [2, 5], [3], [4]]
    @test indices(collect("ahsgh"),'h') == [[2, 5]]
    #@test indices_cnt(collect("ahsgh")) == [1, 2, 1, 1]
    #@test indices_cnt(collect("ahsgh"),'h') == [2]
    #@test permutations_cnt(collect("ahsgh")) == 120
    #@test permutations_cnt(collect("ahsgh");unique=true) == 60
end

