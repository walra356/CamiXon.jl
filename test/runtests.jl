using CamiXon
using Test

@testset "CamiXon.jl" begin
    @test find_all([:📑, :📌,:📢,:📌,:📞]) == [[1], [2, 4], [3], [5]]
    @test find_all([:📑, :📌,:📢,:📌,:📞]; count=true) == [1, 2, 1, 1]
    @test find_all([:📑, :📌,:📢,:📌,:📞], :📌) == [[2, 4]]
    @test find_all([:📑, :📌,:📢,:📌,:📞], :📌; count=true) == [2]
    @test find_all("aβcβd") == [[1], [2, 4], [3], [5]]
    @test find_all("aβcβd"; count=true) == [1, 2, 1, 1]
    @test find_all("aβcβd", 'β') == [[2, 4]]
    @test find_all("aβcβd", 'β'; count=true) == [2]
    @test find_first([:📑,:📌,:📢,:📌,:📞]) == [1,2,3,5]
    @test find_first([:📑,:📌,:📢,:📌,:📞]; dict=true) == [:📑 => 1,:📌 => 2,:📢 => 3,:📞 => 5]
    @test find_first([:📑,:📌,:📢,:📌,:📞], :📌) == [2]
    @test find_first([:📑,:📌,:📢,:📌,:📞], :📌; dict=true) == [:📌 => 2]
    @test find_first([:📑,:📌,:📢,:📌,:📞]) == find_first([1,2,3,2,5]) == find_first("aβcβd")
    @test find_last([:📑,:📌,:📢,:📌,:📞]) == [1,4,3,5]
    @test find_last([:📑,:📌,:📢,:📌,:📞]; dict=true) == [:📑 => 1,:📌 => 4,:📢 => 3,:📞 => 5]
    @test find_last([:📑,:📌,:📢,:📌,:📞], :📌) == [4]
    @test find_last([:📑,:📌,:📢,:📌,:📞], :📌; dict=true) == [:📌 => 4]
    @test find_last([:📑,:📌,:📢,:📌,:📞]) == find_last([1,2,3,2,5]) == find_last("aβcβd")
    @test canonical_partitions(6; header=true) == [[1, 1, 1, 1, 1, 1], [2, 2, 2], [3, 3], [4, 2], [5, 1], [6]]
    @test canonical_partitions(6) == [[1, 1, 1, 1, 1, 1], [2, 2, 2], [3, 3], [4, 2], [5, 1], [6]]
    @test canonical_partitions(6; header=true, reverse=false) == [[6], [5, 1], [4, 2], [3, 3], [2, 2, 2], [1, 1, 1, 1, 1, 1]]
    @test integer_partitions(7) == [[1, 1, 1, 1, 1, 1, 1], [2, 2, 2, 1], [2, 2, 1, 1, 1], [2, 1, 1, 1, 1, 1], [3, 3, 1], [3, 2, 2], [3, 2, 1, 1], [3, 1, 1, 1, 1], [4, 3], [4, 2, 1], [4, 1, 1, 1], [5, 2], [5, 1, 1], [6, 1], [7]]
    @test integer_partitions(9;count=true) == 30
    @test integer_partitions(9) == unique(integer_partitions(9))
    @test integer_partitions(7,4) == [[4, 3], [4, 2, 1], [4, 1, 1, 1]]
    @test integer_partitions(7,4;count=true) == 3
    @test integer_partitions(7; transpose=true) == [[7], [4, 3], [5, 2], [6, 1], [3, 2, 2], [3, 3, 1], [4, 2, 1], [5, 1, 1], [2, 2, 2, 1], [3, 2, 1, 1], [4, 1, 1, 1], [2, 2, 1, 1, 1], [3, 1, 1, 1, 1], [2, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1]]
    @test log10_mantissa.([3,30,300]) == [0.47712125471966244, 0.4771212547196624, 0.4771212547196626]
    @test log10_characteristic_power.([3,30,300]) == [0, 1, 2]
    @test step125.([5,10,21.3,50,100.1]) == [1, 2, 5, 10, 20]
    @test select125(1:100) == [20, 40, 60, 80, 100]
    @test ticks125.([5,10,21.3,50,100.1]) == [-5:1:5, -10:2:10, -20:5:20, -50:10:50, -100:20:100]
    @test centerticks(Base.OneTo(100)) == -100:20:100
    @test centerticks(UnitRange(-200:100)) == -200:50:200
    @test centerticks(-200..100) == -200:50:200
    @test centerticks([1,4,2,5]) == Float32[0.5, 3.0, 6.0, 9.5]
    @test edgeticks(Base.OneTo(100)) == -100:20:100
    @test edgeticks(UnitRange(-200:100)) == -200:50:200
    @test edgeticks(-200..100) == -200:50:200
    @test edgeticks([1,4,2,5]) == [0, 1, 5, 7, 12]
    @test centers(Base.OneTo(100)) == Base.OneTo(100)
    @test centers(UnitRange(-200:100)) == -200:100
    @test centers(-200..100) == -200..100
    @test centers([1,4,2,5]) == [0, 1, 5, 7, 12
    @test edges(Base.OneTo(100)) == 0.5..99.5
    @test edges(UnitRange(-200:100)) == 0.5..99.5
    @test edges(-200..100) == -200.5..99.5
    @test edges([1,4,2,5]) == [0, 1, 5, 7, 12]
    #@test fits_info("T01.fits") == "T01.fits: file was found (for more information set info=true)"
    #@test fits_copy("T01.fits") == "T01.fits was saved as T01 - Copy.fits"
    #@test fits_copy("T01.fits","T01a.fits";protect=false) == "T01.fits was saved as T01a.fits"
    @test fits_create()
    @test fits_extend()
    @test fits_read()
    @test fits_add_key()
    @test fits_edit_key()
    @test fits_delete_key()
    @test fits_rename_key()
end
