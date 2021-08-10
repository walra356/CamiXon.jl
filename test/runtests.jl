using CamiXon
using IntervalSets
using LinearAlgebra
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
    @test f_diff_weight(5, 3) == -10
    @test f_diff_weights(3) == [-1, 3, -3, 1]
    @test f_diff_weights_array(3) ==  [[1], [-1, 1], [1, -2, 1], [-1, 3, -3, 1]]
#@test f_diff_expansion_coeffs_interpolation(2,-1) == [1, -1, 0]
    @test f_diff_expansion_coeffs_interpolation(3,-1) == [1, -1, 0, 0]
    @test f_diff_expansion_coeffs_array_interpolation(2,1) ==  [[1, -2, 1], [1, -1, 0], [1, 0, 0]]
    @test f_diff_expansion_coeffs_extrapolation(5,2) == [1, 2, 3, 4, 5, 6]
    @test f_diff_expansion_weights_extrapolation(5,1) == [-1, 6, -15, 20, -15, 6]
    @test f_diff_expansion_weights(UnitRange(0,5), f_diff_weights_array(5)) == [-5, 29, -69, 85, -55, 15]
#@test f_diff_expansion_weights(UnitRange(0,5), f_diff_weights_array(5)) == [15, -55, 85, -69, 29, -5]
    @test f_diff_expansion_weights_array(7, 3, 1, [[1, 0, 0], [1, -1, 0], [1, -2, 1]]) == [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 1, 0], [1, 0, 0]]
    @test [summation_ranges(7,i,2,1) for i=0:6] == UnitRange{Int64}[1:3, 2:4, 3:5, 4:6, 5:7, 5:7, 5:7]
    @test f_diff_function_sequences([0,1,2,3,4,5,6],2, 1) == [[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [4, 5, 6], [4, 5, 6]]
    @test lagrangian_interpolation([0.0,1,2,3,4,5,6], 0.0..1.0; k=2, i=1) == (0.0:0.08333333333333333:1.0, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0])
    #@test f_diff_expansion_coeffs_differentiation(2,0) == [0.0, 1.0, -1.5]
    #@test f_diff_expansion_coeffs_array_differentiation(2,2) == [[0.0, 1.0, 0.5], [0.0, 1.0, 0.0], [0.0, 1.0, -0.5], [0.0, 1.0, -1.0], [0.0, 1.0, -1.5]]
    #@test f_diff_expansion_coeffs_array_differentiation(2,1) == [[0.0, 1.0, 0.5], [0.0, 1.0, -0.5], [0.0, 1.0, -1.5]]
    #@test lagrangian_differentiation([0.0,1,2,3,4,5], 0.0..5.0; k=2) == (0.0:1.0:5.0, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    @test f_diff_expansion_coeffs_adams_moulton(5) == Rational[1//1, -1//2, -1//12, -1//24, -19//720, -3//160]
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
    @test polynom_deriv_coeffs([1.0 for i=1:6],2) == [-0.0, 0.0, 2.0, 6.0, 12.0, 20.0]
    @test polynom([1.0 for i=1:6],2.0) == 63.0
    @test polynom_multiplication_coeffs([1,1],[1,- 1]) == [1, 0, -1]
    @test polynom_multiplication_coeffs([1,1],[1,- 1,2]) == [1, 0, 1, 2]
    @test polynom_multiplication_coeffs([1,- 1,2],[1,1]) == [1, 0, 1, 2]
    @test permutations_unique_count([[1,2,3],[2,3,1,4,3]],2) == 60
    @test edges(1:5,2.5,2.5) == [-1.25, 1.25, 3.75, 6.25, 8.75]
    @test steps([4,2,6]) == [0, 4, 6, 12]
    @test stepcenters([4,2,6]) == [2.0, 5.0, 9.0]
    @test stepedges([4,2,6]) == [0, 4, 6, 12]
    @test select125([1,2,4,6,8,10,13,16,18,20,40,60,80,100]) == [2, 6, 10, 16, 20, 60, 100]
    @test step125.([5,10,21.3,50,100.1]) == [1, 2, 5, 10, 20]
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
