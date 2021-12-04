using CamiXon
using IntervalSets
using LinearAlgebra
using Test

@testset "CamiXon.jl" begin
    @test sup(-5//2) == "⁻⁵ᐟ²"
    @test sub(-5//2) == "₋₅⸝₂"
    @test frac(-5//2) == "-⁵/₂"
    @test mendeleev(11) == ("Sodium", "Na")
    @test Atom("Helium ion", "⁴Heᐩ", 2, 1, 4.0026, 1//2, 0.0) == Atom("Helium ion", "⁴Heᐩ", 2, 1, 4.0026, 1//2, 0.0)
    @test createAtom(1,0,1.008,1//2,5.59) == Atom("Hydrogen", "¹H", 1, 0, 1.008, 1//2, 5.59)
    @test Term("1s ²S₁⸝₂", 1, 0, 1//2, 0, 1//2) == Term("1s ²S₁⸝₂", 1, 0, 1//2, 0, 1//2)
    @test createTerm(1,0,1//2,0,1//2) == Term("1s ²S₁⸝₂", 1, 0, 1//2, 0, 1//2)
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
    @test f_diff_expansion_coeffs_lagrange(3,-1) == [1, -1, 0, 0]
    @test f_diff_expansion_coeffs_lagrange(5,2) == [1, 2, 3, 4, 5, 6]
    @test f_diff_expansion_weights(UnitRange(0,5), f_diff_weights_array(5)) == [-5, 29, -69, 85, -55, 15]
    @test [summation_range(7,i,2,1) for i=0:6] == UnitRange{Int64}[1:3, 2:4, 3:5, 4:6, 5:7, 5:7, 5:7]
    @test f_diff_function_sequences([0,1,2,3,4,5,6],2) == [[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [4, 5, 6], [4, 5, 6]]
    @test lagrange_interpolation([0.0,1,2,3,4,5,6], 0.0..1.0; k=2, m=2) == (0.0:0.08333333333333333:1.0, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0])
    @test lagrange_extrapolation([0.0,1,2,3,4,5,6,7], 0.0..1.0; k=2, e=1, m=2) == (1.0:0.07142857142857142:1.1428571428571428, [7.0, 7.5, 8.0])
    @test f_diff_expansion_coeffs_differentiation(2,0) == [0.0, 1.0, 0.5]
    @test lagrange_differentiation([0.0,1,2,3,4,5], 0.0..5.0; k=2, m=1) == (0.0:1.0:5.0, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    @test create_lagrange_differentiation_weights(8,0) == [-761//280, 8//1, -14//1, 56//3, -35//2, 56//5, -14//3, 8//7, -1//8]
    @test create_lagrange_differentiation_matrix(3) == [-11//6 3//1 -3//2 1//3; -1//3 -1//2 1//1 -1//6; 1//6 -1//1 1//2 1//3; -1//3 3//2 -3//1 11//6]
    @test trapezoidal_weights(5; rationalize=true) == [95//288, 317//240, 23//30, 793//720, 157//160]
    @test trapezoidal_integration([1.0, 4.0, 15.0, 40.0, 85.0, 156.0], 0.0..5.0, [3//8, 7//6, 23//24]) ≈ 215.4166666
    @test create_adams_moulton_weights(3;rationalize=true) == [1//24, -5//24, 19//24, 3//8]
    @test f_diff_expansion_coeffs_adams_moulton(5) == [1//1, -1//2, -1//12, -1//24, -19//720, -3//160]
    @test f_diff_expansion_coeffs_adams_bashford(5) == [1//1, 1//2, 5//12, 3//8, 251//720, 95//288]
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
    @test bernoulli_numbers(10) == [1//1, -1//2, 1//6, 0//1, -1//30, 0//1, 1//42, 0//1, -1//30, 0//1, 5//66]
    @test faulhaber_polynom(6) == [0//1, 0//1, -1//12, 0//1, 5//12, 1//2, 1//6]
    @test faulhaber_summation(3,5) == 276
    @test harmonic_number(12) == 86021//27720
    @test harmonic_number(12, 3) == 25535765062457//21300003648000
    @test pascal_triangle(5) == [[1], [1, 1], [1, 2, 1], [1, 3, 3, 1], [1, 4, 6, 4, 1], [1, 5, 10, 10, 5, 1]]
    @test pascal_next([1, 4, 6, 4, 1]) == [1, 5, 10, 10, 5, 1]
    @test polynomial([1,1,1,1,1,1],1; deriv=-1) == 49//20
    @test polynom_derivative([1,1,1,1,1]) == [1,2,3,4]
    @test polynom_derivatives([1,1,1,1,1]; deriv=2) == [2, 6, 12]
    @test polynom_derivatives_all([1,1,1,1,1]) == [[1, 2, 3, 4], [2, 6, 12], [6, 24], [24]]
    @test polynom_power([1,1,1],2) == [1, 2, 3, 2, 1]
    @test polynom_powers([1,1,1],3) == [[1, 1, 1], [1, 2, 3, 2, 1], [1, 3, 6, 7, 6, 3, 1]]
    @test polynom_primitive([1,1,1,1,1]) == [0//1, 1//1, 1//2, 1//3, 1//4, 1//5]
    @test polynom_product([1,1],[1,- 1]) == [1, 0, -1]
    @test polynom_product([1,1],[1,- 1,2]) == [1, 0, 1, 2]
    @test polynom_product([1,- 1,2],[1,1]) == [1, 0, 1, 2]
    @test polynom_product_expansion([1,-1,1],[1,1,-1,1,1,1], 5) == [1, 0, -1, 3, -1, 1]
    @test permutations_unique_count([[1,2,3],[2,3,1,4,3]],2) == 60
    @test normalize_VectorRational([2//3,4//5]).num == [10,12]
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
