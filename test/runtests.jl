using CamiXon
using Test

@testset "CamiXon.jl" begin
    @test indices(collect("ahsgh")) == [[1], [2, 5], [3], [4]]
    @test indices(collect("ahsgh"),'h') == [[2, 5]]
#end

