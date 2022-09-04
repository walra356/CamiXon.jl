using CamiXon
#using IntervalSets
using LinearAlgebra
using Test

@testset "CamiXon.jl" begin
    atom = castAtom(Z=1, A=1, Q=0);
    orbit = castOrbit(n=1, ℓ=0);
    codata = castCodata(2018);
    grid = autoGrid(atom, orbit, codata, Float64);
    def = castDef(grid, atom, orbit);
    @test grid.name == "exponential"
    @test def.atom.element.name == "hydrogen"
    @test sup(-5//2) == "⁻⁵ᐟ²"
    @test sub(-5//2) == "₋₅⸝₂"
    @test frac(-5//2) == "-⁵/₂"
    @test strRational(-5) == "-5"
    @test bohrformula(2, 4) == -1//8
    @test get(dictAtomicNumbers, "Rb", nothing) == 37
    @test get(dictElements, 37, nothing) == ("rubidium", "Rb", 85.468)
    @test castElement(Z=1, msg=false) == Element("hydrogen", "H", 1.008)
    @test castElement("Rb"; msg=false) == castElement(Z=37, msg=false)
    @test listElement("H") == listElement(1)
    @test listElement(1; fmt=String) == "H, hydrogen, Z=1, weight=1.008"
    @test listElements(1,3) == listElements(1:3)
    @test listElements(1:2; fmt=String) == ["H, hydrogen, Z=1, weight=1.008", "He, helium, Z=2, weight=4.0026" ]
    @test castIsotope(Z=1,A=1,msg=false) == Isotope("¹H", "hydrogen", 1, 1, 0, 0.8783, 1.007825032, 1//2, 1, 1.0e100, 2.792847351, 0.0, 99.9855)
    @test castIsotope("Rb"; A=87, msg=false) == castIsotope(Z=37, A=87, msg=false)
    @test listIsotope(2,3; fmt=Latex) == "2 & helium & \$^{3}\$He & 3\\, & 1 & 1.9661 & 3.016029322 & 1/2\$^+\$ & -2.12762531 & 0.0 & 0.0002 \\\\\n"
    @test listIsotope(2,3; fmt=String) == "³He, helium, Z=2, A=3, N=1, R=1.9661, M=3.016029322, I=1/2⁺, μI=-2.12762531, Q=0.0, RA=0.0002%, (stable)"
    @test listAtom("H", 3, 0) == listAtom(1, 3, 0)
    @test listAtom(1,3,1; fmt=String) == "tritium ion, ³Tᐩ, Z=1, A=3, Q=1, Zc=2"
    @test listAtoms(2:2,0; fmt=String) == ["helium, neutral atom, ³He, Z=2, A=3, Q=0, Zc=1", "helium, neutral atom, ⁴He, Z=2, A=4, Q=0, Zc=1"]
    @test castAtom(Z=1, A=1, Q=0, msg=false) == Atom(1, 1, 0, 1, Element("hydrogen", "H", 1.008), Isotope("¹H", "hydrogen", 1, 1, 0, 0.8783, 1.007825032, 1//2, 1, 1.0e100, 2.792847351, 0.0, 99.9855))
    @test castAtom("Rb"; A=87, Q=0, msg=false) == castAtom(Z=37, A=87, Q=0, msg=false)
    @test createSpinOrbit(orbit; msg=false) == SpinOrbit("1s↑", 1, 0, 0, 1//2)
    @test Term("1s ²S₁⸝₂", 1, 0, 0, 1//2, 0, 1//2) == Term("1s ²S₁⸝₂", 1, 0, 0, 1//2, 0, 1//2)
    @test createTerm(1; ℓ=0, S=1//2, L=0, J=1//2, msg=false) == Term("1s ²S₁⸝₂", 1, 0, 0, 1//2, 0, 1//2)
    @test convertUnit(1, codata; unitIn="Hz", unitOut="Joule") == Value(6.62607015e-34, "Joule")
    @test convertUnit(1, codata) == Value(6.579683920501762, "PHz")
    @test strValue(Value(1,"Hz")) == "1 Hz"
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
    @test [fdiff_weight(5,j) for j=0:5] == [1, -5, 10, -10, 5, -1]
    @test isforward(fwd) == true
    @test fdiff_interpolation_expansion_coeffs(-1, 5) == [1, 1, 1, 1, 1, 1]
    @test fdiff_interpolation_expansion_coeffs(-1, 5, bwd) == [1, 1, 1, 1, 1, 1]
    @test fdiff_interpolation_expansion_coeffs(1, 5, fwd) == [1, -1, 1, -1, 1, -1]
    @test [fdiff_interpolation([ν^3 for ν=-5:2], ν; k=3) for ν=1:0.5:8] == [-125.0, -91.125, -64.0, -42.875, -27.0, -15.625, -8.0, -3.375, -1.0, -0.125, 0.0, 0.125, 1.0, 3.375, 8.0]
    @test fdiff_expansion_weights([0,1,2,3,4,5],fwd,reg) ==  [-3, 15, -33, 37, -21, 5]
    @test fdiff_expansion_weights([0,1,2,3,4,5],bwd,rev) == [-5, 29, -69, 85, -55, 15]
    @test fdiff_expansion_weights([0,1,2,3,4,5]) == [-5, 29, -69, 85, -55, 15]
    @test fdiff_expansion([1,-1,1,-1],[1,4,9,16],fwd) == 0
    @test fdiff_expansion([1,1,1,1],[1,4,9,16],bwd) == 25
    @test fdiff_expansion([1,1,1,1],[1,4,9,16]) == 25
    @test fdiff_expansion_coeffs_differentiation(0,2) == [0.0, 1.0, 0.5]
    @test fdiff_differentiation([0.0, 1, 4, 9, 16, 25]; k=3) ≈ [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    @test create_lagrange_differentiation_matrix(3) == [-11//6 3//1 -3//2 1//3; -1//3 -1//2 1//1 -1//6; 1//6 -1//1 1//2 1//3; -1//3 3//2 -3//1 11//6]
    @test trapezoidal_weights(5; rationalize=true) == [95//288, 317//240, 23//30, 793//720, 157//160]
    @test trapezoidal_integration([1.0, 4.0, 15.0, 40.0, 85.0, 156.0], 0.0, 5.0, [3//8, 7//6, 23//24]) ≈ 215.4166666
    @test create_adams_moulton_weights(3;rationalize=true) == [1//24, -5//24, 19//24, 3//8]
    @test fdiff_expansion_coeffs_adams_moulton(5) == [1//1, -1//2, -1//12, -1//24, -19//720, -3//160]
    @test fdiff_expansion_coeffs_adams_bashford(5) == [1//1, 1//2, 5//12, 3//8, 251//720, 95//288]
    @test gridname(2) == "quasi-exponential"
    @test [gridfunction(2, n-1, 0.1; p = 1) for n=1:5] ==  [0.0, 0.10000000000000009, 0.19999999999999996, 0.30000000000000004, 0.3999999999999999]
    @test [gridfunction(1, n-1, 0.1) for n=1:4] == [0.0, 0.10517091807564771, 0.22140275816016985, 0.3498588075760032]
    @test [gridfunction(2, n-1, 0.1; p = 4) for n=1:4] == [0.0, 0.10517083333333321, 0.22140000000000004, 0.3498375]
    @test [gridfunction(3, n-1, 0.1; coords=[0,1,1/2,1/6,1/24]) for n=1:3] == [0.0, 0.10517083333333334, 0.2214]
    @test castGrid(2, 3, Float64; p=1, h=0.1, r0=1.0, msg=false).r == [0.0, 0.10000000000000009, 0.19999999999999996]
    @test castGrid(1, 3, Float64; h=0.1, r0=1.0, msg=false).r == [0.0, 0.10517091807564771, 0.22140275816016985]
    @test castGrid(2, 3, Float64; p=4, h=0.1, r0=1.0, msg=false).r == [0.0, 0.10517083333333321, 0.22140000000000004]
    @test castGrid(3, 3, Float64; coords=[0,1,1/2,1/6,1/24], h=0.1, r0=1.0, msg=false).r == [0.0, 0.10517083333333334, 0.2214]
    grid = castGrid(4, 6, Float64; r0=1.0, h=1.0, msg=false);
    @test grid_differentiation([0.0, 1, 4, 9, 16, 25], grid; k=3) ≈ [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    @test autoRmax(atom, orbit) == 63.0
    @test autoNtot(orbit) == 100
    @test autoPrecision(100.0, orbit) == Float64
    @test autoSteps(1, 100, 100) == (0.1, 0.004540199100968777)
    @test grid_differentiation([0.0, 1, 4, 9, 16, 25], grid) ≈ [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    @test grid_trapezoidal_integral([0.,1.,2.,3.,4.], 1:5, castGrid(2, 5, Float64; p=1, msg=false)) == 0.008
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
    @test pochhammer.([0,1,2,3,4],5) == [0, 120, 720, 2520, 6720]
    @test polynomial([1,1,1,1,1,1],1; deriv=-1) == 49//20
    @test polynom_derivative([1,1,1,1,1]) == [1,2,3,4]
    @test polynom_derivatives([1,1,1,1,1]; deriv=2) == [2, 6, 12]
    @test polynom_derivatives_all([1,1,1,1,1]) == [[1, 2, 3, 4], [2, 6, 12], [6, 24], [24]]
    @test polynom_power([1,1,1],2) == [1, 2, 3, 2, 1]
    @test polynom_powers([1,1,1],3) == [[1, 1, 1], [1, 2, 3, 2, 1], [1, 3, 6, 7, 6, 3, 1]]
    @test polynom_primitive([1,1,1,1,1]) == [0//1, 1//1, 1//2, 1//3, 1//4, 1//5]
    @test polynom_product([1.0,-1],[1,2,-3]) == [1.0, 1.0, -5.0, 3.0]
    @test polynom_product([1//2,-1],[1,2,-3]) == [1//2, 0//1, -7//2, 3//1]
    @test polynom_product([1,1],[1,- 1,2]) == [1, 0, 1, 2]
    @test polynom_product_expansion([1,-1,1],[1,1,-1,1,1,1], 5) == [1, 0, -1, 3, -1, 1]
    @test permutations_unique_count([[1,2,3],[2,3,1,4,3]],2) == 60
    @test normalize_VectorRational([2//3,4//5]).num == [10,12]
    @test edges(1:5,2.5,2.5) == [-1.25, 1.25, 3.75, 6.25, 8.75]
    @test steps([4,2,6]) == [0, 4, 6, 12]
    @test stepcenters([4,2,6]) == [2.0, 5.0, 9.0]
    @test stepedges([4,2,6]) == [0, 4, 6, 12]
    @test select125([1,2,4,6,8,10,13,16,18,20,40,60,80,100]) == [2, 6, 10, 16, 20, 60, 100]
    @test step125.([5,10,21.3,50,100.1]) == [1, 2, 5, 10, 20]
    @test texp(1.0, 0.0, 5) == 2.7166666666666663
    @test texp(1, 0, 5) == 163//60
    @test texp(1//1, 0//1, 5) == 163//60
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
