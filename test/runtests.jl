using CamiXon
using Test

@testset "CamiXon.jl" begin
    @test indices(collect("ahsgh")) == [[1], [2, 5], [3], [4]]
    @test indices(collect("ahsgh"),'h') == [[2, 5]]
    @test indices_cnt(collect("ahsgh")) == [1, 2, 1, 1]
    @test indices_cnt(collect("ahsgh"),'h') == [2]
    @test canonical_partitions([6],1; trailer=true)==[[1, 1, 1, 1, 1, 1], [2, 2, 2], [3, 3], [4, 2], [5, 1], [6]]
    @test canonical_partitions([6],1)==[[1, 1, 1, 1, 1, 1], [2, 2, 2], [3, 3], [4, 2], [5, 1]]
    @test canonical_partitions([4,4,4,4],3; trailer=true)==[[4, 4, 1, 1, 1, 1, 1, 1, 1, 1], [4, 4, 2, 2, 2, 2], [4, 4, 3, 3, 2], [4, 4, 4, 4]]
    @test next_partitions([[4,4,4,4]],2)==[[4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [4, 2, 2, 2, 2, 2, 2], [4, 3, 3, 3, 3]]
    @test all_partitions(7)==[[1, 1, 1, 1, 1, 1, 1], [2, 2, 2, 1], [3, 3, 1], [4, 3], [5, 2], [6, 1], [7], [2, 1, 1, 1, 1, 1], [3, 1, 1, 1, 1], [3, 2, 2], [4, 1, 1, 1], [4, 2, 1], [5, 1, 1], [2, 2, 1, 1, 1], [3, 2, 1, 1]]
    @test permutations_cnt(collect("ahsgh")) == 120
    @test permutations_cnt(collect("ahsgh");unique=true) == 60
end
