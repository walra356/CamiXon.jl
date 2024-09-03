# SPDX-License-Identifier: MIT

using CamiXon
import CamiMath

# using IntervalSets
#using BenchmarkTools
using LinearAlgebra
using Test

@testset "CamiXon.jl" begin 
    codata = castCodata(2022)
    codata = castCodata(2018)
    atom = castAtom(Z=1, A=1, Q=0)
    orbit = castOrbit(n=2, ℓ=0; msg=false)
    grid = autoGrid(atom, orbit, Float64; Nboost=1, msg=false)
    def = castDef(grid, atom, orbit, codata)
    RH2s_example = [RH2s(1, grid.r[n]) for n = 1:grid.N]
    ZH2s_example = reduce_wavefunction(RH2s_example, grid)
    ZH2s_generic = hydrogenic_reduced_wavefunction(1, orbit, grid)
    @test ZH2s_example ≈ ZH2s_generic
    orbit = castOrbit(n=2, ℓ=1; msg=false)
    grid = autoGrid(atom, orbit, Float64; Nboost=1, msg=false)
    def = castDef(grid, atom, orbit, codata)
    RH2p_example = [RH2p(1, grid.r[n]) for n = 1:grid.N]
    ZH2p_example = reduce_wavefunction(RH2p_example, grid)
    ZH2p_generic = hydrogenic_reduced_wavefunction(1, orbit, grid)
    @test ZH2p_example ≈ ZH2p_generic
    orbit = castOrbit(n=1, ℓ=0)
    grid = autoGrid(atom, orbit, Float64)
    def = castDef(grid, atom, orbit, codata)
    RH1s_example = [RH1s(atom.Z, grid.r[n]) for n = 1:grid.N]
    ZH1s_generic = hydrogenic_reduced_wavefunction(1, orbit, grid)
    ZH1s_example = reduce_wavefunction(RH1s_example, grid)
    @test ZH1s_example ≈ ZH1s_generic
    E = initE(def)
    #adams = castAdams(E, grid, def)
    E, Z = adams_moulton_master(E, grid, def; Δν=Value(1, "kHz"), imax=25, msg=false)
    @test ((real(ZH1s_example .- Z)) .< [1.0e-7 for i = 1:grid.N]) == ones(Bool, grid.N)
    @test ((imag(ZH1s_example .- Z)) .< [1.0e-7 for i = 1:grid.N]) == ones(Bool, grid.N)
    #Z1 = hydrogenic_reduced_wavefunction(1, orbit, grid);
    #P = real(Z)
    #val = UF(0, P, grid)[1];
    @test round(Int, UF(0, real(Z), grid)[1]) == 1
    @test calibrationReport(1.1, 1.0, codata; unitIn="Hartree", msg=false) == "calibration report (Float64):\nEcal = 1 Hartree \nE = 1.1000000000000001 Hartree \nabsolute accuracy: ΔE = 0.1 Hartree (657.968 THz)\nrelative accuracy: ΔE/E = 0.0909091\n"
    @test listCodata(codata; msg=false) == "∆νCs = 9192631770 Hz      - ¹³³Cs hyperfine transition frequency\n   c = 299792458 m s⁻¹    - speed of light in vacuum\n   h = 6.62607e-34 J Hz⁻¹ - Planck constant\n   ħ = 1.05457e-34 J s    - Planck constant (reduced)\n   e = 1.60218e-19 C      - elementary charge\n  kB = 1.38065e-23 J K⁻¹  - Boltzmann constant\n  NA = 6.02214e23 mol⁻¹   - Avogadro constant\n Kcd = 683 lm W⁻¹         - Luminous efficacy\n  mₑ = 9.10938e-31 Kg     - electron rest mass\n  R∞ = 1.09737e7 m⁻¹      - Rydberg constant\n  Ry = 3.28984e15 Hz      - Rydberg frequency\n  Eₕ = 4.35974e-18 J      - Hartree atomic unit\n   α = 0.00729735         - fine-structure constant\n  μ₀ = 1.25664e-6 N A⁻²   - magnetic permitivity of vacuum\n  ε₀ = 8.85419e-12 F m⁻¹  - electric permitivity of vacuum\n  KJ = 4.83598e14 Hz V⁻¹  - Josephson constant\n  RK = 25812.8 Ω          - Von Klitzing constant\n   R = 8.31446 J mol⁻¹K⁻¹ - Molar gas constant\n   u = 1.66054e-27 kg     - unified atomic mass unit\n" 
    @test grid.name == "exponential"
    @test findIndex(0.0042, grid) == 9
    @test def.atom.element.name == "hydrogen"
    @test def.pos.Na == 8
    @test def.pos.Nb == 103
    @test get_Na(Z, def) == 8
    @test get_Nb(Z, def) == 103
    @test get_Nuctp(E, def) == 76
    @test grid_integration(real(Z) .^ 2, 1, grid.N, grid) ≈ 1.0
    @test grid_integration(real(ZH1s_generic) .^ 2, 1, grid.N, grid) ≈ 1.0
    @test sup(-5 // 2) == "⁻⁵ᐟ²"
    @test sub(-5 // 2) == "₋₅⸝₂"
    @test sub("e") == "ₑ"
    @test frac(-5 // 2) == "-⁵/₂"
    @test strRational(-5) == "-5"
    @test strRational(-5 // 2) == "-5/2"
    @test bohrformula(2, 4) == -1 // 8
    @test get(dictAtomicNumbers, "Rb", nothing) == 37
    @test get(dictElements, 37, nothing) == ("rubidium", "Rb", 85.468)
    @test castElement(Z=1, msg=false) == Element("hydrogen", "H", 1.008)
    @test castElement("Rb"; msg=false) == castElement(Z=37, msg=false)
    @test listElement("H") == listElement(1)
    @test listElement(1; fmt=String) == "H, hydrogen, Z=1, weight=1.008"
    @test listElement(1; fmt=Info, msg=false) == "Element: hydrogen\n  symbol: H\n  atomic number: Z = 1\n  atomic weight (relative atomic mass): 1.008"
    @test listElements(1, 3) == listElements(1:3)
    @test listElements(1:2; fmt=String) == ["H, hydrogen, Z=1, weight=1.008", "He, helium, Z=2, weight=4.0026"]
    @test castIsotope(Z=1, A=1, msg=false) == Isotope("¹H", "hydrogen", 1, 1, 0, 0.8783, 1.007825032, 1 // 2, 1, 1.0e100, 2.792847351, 0.0, 99.9855)
    @test castIsotope("Rb"; A=87, msg=false) == castIsotope(Z=37, A=87, msg=false)
    @test listIsotope(2, 3; fmt=Latex) == "2 & helium & \$^{3}\$He & 3\\, & 1 & 1.9661 & 3.016029322 & 1/2\$^+\$ & -2.12762531 & 0.0 & 0.0002 \\\\\n"
    @test listIsotope(2, 3; fmt=String) == "³He, helium, Z=2, A=3, N=1, R=1.9661, M=3.016029322, I=1/2⁺, μI=-2.12762531, Q=0.0, RA=0.0002%, (stable)"
    @test listAtom("H", 3, 0) == listAtom(1, 3, 0)
    @test listAtom(1, 3, 1; fmt=String) == "tritium ion, ³Tᐩ, Z=1, A=3, Q=1, Zc=2"
    @test listAtoms(2:2, 0; fmt=String) == ["helium, neutral atom, ³He, Z=2, A=3, Q=0, Zc=1", "helium, neutral atom, ⁴He, Z=2, A=4, Q=0, Zc=1"]
    @test castAtom(Z=1, A=1, Q=0, msg=false) == Atom(1, 1, 0, 1, Element("hydrogen", "H", 1.008), Isotope("¹H", "hydrogen", 1, 1, 0, 0.8783, 1.007825032, 1 // 2, 1, 1.0e100, 2.792847351, 0.0, 99.9855))
    @test castAtom("Rb"; A=87, Q=0, msg=false) == castAtom(Z=37, A=87, Q=0, msg=false)
    @test castSpinorbit(n=1, ℓ=0, msg=false) == Spinorbit("1s↑", Orbit("1s", 1, 0, 0, 0), 1//2)
    @test Term("1s ²S₁⸝₂", 1, 0, 0, 1 // 2, 0, 1 // 2) == Term("1s ²S₁⸝₂", 1, 0, 0, 1 // 2, 0, 1 // 2)
    @test createTerm(1; ℓ=0, S=1 // 2, L=0, J=1 // 2, msg=false) == Term("1s ²S₁⸝₂", 1, 0, 0, 1 // 2, 0, 1 // 2)
    @test convertUnit(1, codata; unitIn="Hz", unitOut="J") == Value(6.62607015e-34, "J")
    @test convertUnit(1, codata) == Value(6.579683920501762, "PHz")
    @test strValue(Value(1, "Hz")) == "1 Hz"
    @test lc_eltype(([1 // 2, 1 // 3]; (1 // 4, 1 // 1, 1 // 6))) == Rational{Int} 
    @test lc_primitivetype(([1 // 2, 1 // 3]; (1 // 4, 1 // 1, 1 // 6))) == Int64
    @test primitivetype(Rational{UInt16}) == UInt16
    @test conditionalType(47, 46) == BigInt
    @test typeof(bigconvert([[1 // 1, 1 // 2], [1 // 1, 1 // 2]])) == Vector{Vector{Rational{BigInt}}}
    @test find_all([:📑, :📌, :📢, :📌, :📞]) == [[1], [2, 4], [3], [5]]
    @test find_all([:📑, :📌, :📢, :📌, :📞]; count=true) == [1, 2, 1, 1]
    @test find_all([:📑, :📌, :📢, :📌, :📞], :📌) == [[2, 4]]
    @test find_all([:📑, :📌, :📢, :📌, :📞], :📌; count=true) == [2]
    @test find_all("aβcβd") == [[1], [2, 4], [3], [5]]
    @test find_all("aβcβd"; count=true) == [1, 2, 1, 1]
    @test find_all("aβcβd", 'β') == [[2, 4]]
    @test find_all("aβcβd", 'β'; count=true) == [2]
    @test find_first([:📑, :📌, :📢, :📌, :📞]) == [1, 2, 3, 5]
    @test find_first([:📑, :📌, :📢, :📌, :📞]; dict=true) == [:📑 => 1, :📌 => 2, :📢 => 3, :📞 => 5]
    @test find_first([:📑, :📌, :📢, :📌, :📞], :📌) == [2]
    @test find_first([:📑, :📌, :📢, :📌, :📞], :📌; dict=true) == [:📌 => 2]
    @test find_first([:📑, :📌, :📢, :📌, :📞]) == find_first([1, 2, 3, 2, 5]) == find_first("aβcβd")
    @test find_last([:📑, :📌, :📢, :📌, :📞]) == [1, 4, 3, 5]
    @test find_last([:📑, :📌, :📢, :📌, :📞]; dict=true) == [:📑 => 1, :📌 => 4, :📢 => 3, :📞 => 5]
    @test find_last([:📑, :📌, :📢, :📌, :📞], :📌) == [4]
    @test find_last([:📑, :📌, :📢, :📌, :📞], :📌; dict=true) == [:📌 => 4]
    @test find_last([:📑, :📌, :📢, :📌, :📞]) == find_last([1, 2, 3, 2, 5]) == find_last("aβcβd")
    @test [fdiff_weight(5, j) for j = 0:5] == [1, -5, 10, -10, 5, -1]
    @test isforward(fwd) == true
    @test isregular(reg) == true
    @test a_direct(2, 1, 1, 2, 2) == 2 // 35
    @test b_exchange(1, 1, 1, 2, 2) == 2 // 5
    @test a_direct(6, 3, 2, 3, -1) == -250 // 20449
    @test b_exchange(6, 3, 2, 3, -1) == 1050 // 20449
    @test fdiff_interpolation_expansion_coeffs(-1, 5) == [1, 1, 1, 1, 1, 1]
    @test fdiff_interpolation_expansion_coeffs(-1, 5, bwd) == [1, 1, 1, 1, 1, 1]
    @test fdiff_interpolation_expansion_coeffs(1, 5, fwd) == [1, -1, 1, -1, 1, -1]
    @test [fdiff_interpolation([ν^3 for ν = -5:2], ν; k=3) for ν = 1:0.5:8] == [-125.0, -91.125, -64.0, -42.875, -27.0, -15.625, -8.0, -3.375, -1.0, -0.125, 0.0, 0.125, 1.0, 3.375, 8.0]
    @test fdiff_expansion_weights([0, 1, 2, 3, 4, 5], fwd, reg) == [-3, 15, -33, 37, -21, 5]
    @test fdiff_expansion_weights([0, 1, 2, 3, 4, 5], bwd, rev) == [-5, 29, -69, 85, -55, 15]
    @test fdiff_expansion_weights([0, 1, 2, 3, 4, 5]) == [-5, 29, -69, 85, -55, 15]
    @test fdiff_expansion([1, -1, 1, -1], [1, 4, 9, 16], fwd) == 0
    @test fdiff_expansion([1, 1, 1, 1], [1, 4, 9, 16], bwd) == 25
    @test fdiff_expansion([1, 1, 1, 1], [1, 4, 9, 16]) == 25
    @test fdiff_differentiation_expansion_coeffs(0, 3) == [0 // 1, 1 // 1, 1 // 2, 1 // 3]
    @test fdiff_differentiation_expansion_coeffs(1, 3) == [0 // 1, 1 // 1, -1 // 2, -1 // 6]
    @test [fdiff_differentiation([16, 9, 4, 1, 0, 1, 4, 9, 16], v) for v = 1:9] == [-8 // 1, -6 // 1, -4 // 1, -2 // 1, 0 // 1, 2 // 1, 4 // 1, 6 // 1, 8 // 1]
    @test fdiff_differentiation([16, 9, 4, 1, 0, 1, 4, 9, 16], 5.5) == 1.0
    @test create_lagrange_differentiation_matrix(3) == [-11//6 3//1 -3//2 1//3; -1//3 -1//2 1//1 -1//6; 1//6 -1//1 1//2 1//3; -1//3 3//2 -3//1 11//6]
    @test trapezoidal_epw(5; rationalize=true) == [95 // 288, 317 // 240, 23 // 30, 793 // 720, 157 // 160]
    @test trapezoidal_integration([1.0, 4.0, 15.0, 40.0, 85.0, 156.0], 0.0, 5.0, [3 // 8, 7 // 6, 23 // 24]) ≈ 215.4166666
    @test create_adams_moulton_weights(3; rationalize=true) == [1 // 24, -5 // 24, 19 // 24, 3 // 8]
    @test fdiff_adams_moulton_expansion_coeffs(5) == [1 // 1, -1 // 2, -1 // 12, -1 // 24, -19 // 720, -3 // 160]
    @test fdiff_adams_moulton_expansion_coeff(0) == 1//1
    @test fdiff_adams_moulton_expansion_coeff(20; msg=false) == -12365722323469980029 // 4817145976189747200000
    @test fdiff_adams_bashford_expansion_coeffs(5) == [1 // 1, 1 // 2, 5 // 12, 3 // 8, 251 // 720, 95 // 288]
    @test fdiff_adams_bashford_expansion_coeff(0) == 1//1
    @test fdiff_adams_bashford_expansion_coeff(20; msg=false) == 8136836498467582599787//33720021833328230400000
    @test gridname(2) == "quasi-exponential"
    @test [gridfunction(2, n - 1, 0.1; p=1) for n = 1:5] == [0.0, 0.10000000000000009, 0.19999999999999996, 0.30000000000000004, 0.3999999999999999]
    @test [gridfunction(1, n - 1, 0.1) for n = 1:4] == [0.0, 0.10517091807564771, 0.22140275816016985, 0.3498588075760032]
    @test [gridfunction(2, n - 1, 0.1; p=4) for n = 1:4] == [0.0, 0.10517083333333321, 0.22140000000000004, 0.3498375]
    @test [gridfunction(4, n - 1, 0.1; coords=[0, 1, 1 / 2, 1 / 6, 1 / 24]) for n = 1:3] == [0.0, 0.10517083333333334, 0.2214]
    @test castGrid(2, 3, Float64; p=1, h=0.1, r0=1.0, msg=false).r == [eps(Float64), 0.10000000000000009, 0.19999999999999996]
    @test castGrid(1, 3, Float64; h=0.1, r0=1.0, msg=false).r == [eps(Float64), 0.10517091807564771, 0.22140275816016985]
    @test castGrid(2, 3, Float64; p=4, h=0.1, r0=1.0, msg=false).r == [eps(Float64), 0.10517083333333321, 0.22140000000000004]
    @test castGrid(4, 3, Float64; coords=[0, 1, 1 / 2, 1 / 6, 1 / 24], h=0.1, r0=1.0, msg=false).r == [eps(Float64), 0.10517083333333334, 0.2214]
    grid = castGrid(3, 6, Float64; r0=1.0, h=1.0, msg=false)
    @test grid_differentiation([0.0, 1, 4, 9, 16, 25], grid; k=3) ≈ [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    @test autoRmax(atom, orbit) == 84.0 #63.0
    @test autoNtot(orbit, 2) == 240
    @test autoPrecision(100.0, orbit) == Float64
    @test autoSteps(1, 100, 100) == (0.1, 0.004540199100968777)
    @test grid_differentiation([0.0, 1, 4, 9, 16, 25], grid) ≈ [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    @test grid_integration([0.0, 1.0, 2.0, 3.0, 4.0], 1:5, castGrid(2, 5, Float64; p=1, msg=false)) == 0.008
    @test edges(1:5, 2.5, 2.5) == [-1.25, 1.25, 3.75, 6.25, 8.75]
    @test steps([4, 2, 6]) == [0, 4, 6, 12]
    @test stepcenters([4, 2, 6]) == [2.0, 5.0, 9.0]
    @test stepedges([4, 2, 6]) == [0, 4, 6, 12]
    @test select125([1, 2, 4, 6, 8, 10, 13, 16, 18, 20, 40, 60, 80, 100]) == [2, 6, 10, 16, 20, 60, 100]
    @test step125.([5, 10, 21.3, 50, 100.1]) == [1, 2, 5, 10, 20]
    @test latent_heat_vaporization("Yb", 763) == 24170.448513975916
    @test latent_heat_vaporization("Li", 623) == 18473.64020109123
    @test latent_heat_vaporization("Li", 400) == 19134.482122780522
    @test svp("Yb", 763) == 2.2859295292626745
    @test svp("Li", 623) == 0.0015230367024569058
    @test svp("Li", 400) == 7.901690229445235e-11
    @test melting_point("Li") == 453.65
    @test CGC(3, 0, 4, -1, 5, -1; msg=true) == parse(BigFloat, "-0.36364052611670255269921486774521555203216489725107181148303161368088211274565")
    @test threeJsymbol(3, 0, 4, -1, 5, 1; msg=true) == parse(BigFloat, "-0.1096417439724123565166029917781360897459044055433631161836138910409772907333476")
    @test latexIsotopeTable(1:10) == latexIsotopeTable(1,10)
    @test latexIsotopeTable(11:22; continuation=true) == latexIsotopeTable(11,22; continuation=true) 
    @test silvera_goldman_triplet(10) == -7.71843646003074e-6
    @test silvera_goldman_singlet(10) == -8.696045600341206e-6
    @test silvera_goldman_exchange(10) == 9.776091403104656e-7
    grid = castGrid(3,2000,Float64; h=0.01, r0=1, msg=false);    
    @test silvera_goldman_potential(grid; ℓ=0, S=1)[700] == -1.2954953056510744e-6
    @test silvera_goldman_potential(grid; ℓ=0, S=0)[700] == -0.00020738292434731114

end
