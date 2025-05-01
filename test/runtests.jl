# SPDX-License-Identifier: MIT

using CamiXon
using CamiDiff
using CamiMath

# using IntervalSets
#using BenchmarkTools
using LinearAlgebra
using Test 

println("CamiXon.jl | 167 runtests | runtime 40s (estimated) | start")

@testset "CamiXon.jl" begin 
    @test CamiMath.frac(-5 // 2) == "-⁵/₂"
    @test_throws DomainError castCodata(2016)
    codata = castCodata(2018)
    @test listCodata(codata; msg=false) == "∆νCs = 9192631770 Hz      - ¹³³Cs hyperfine transition frequency\n   c = 299792458 m s⁻¹    - speed of light in vacuum\n   h = 6.62607e-34 J Hz⁻¹ - Planck constant\n   ħ = 1.05457e-34 J s    - Planck constant (reduced)\n   e = 1.60218e-19 C      - elementary charge\n  kB = 1.38065e-23 J K⁻¹  - Boltzmann constant\n  NA = 6.02214e23 mol⁻¹   - Avogadro constant\n Kcd = 683 lm W⁻¹         - Luminous efficacy\n  mₑ = 9.10938e-31 kg     - electron mass\n  mₚ = 1.67262e-27 kg     - proton mass\n  R∞ = 1.09737e7 m⁻¹      - Rydberg constant\n  Ry = 3.28984e15 Hz      - Rydberg frequency\n  Eₕ = 4.35974e-18 J      - Hartree atomic unit\n   α = 0.00729735         - fine-structure constant\n  a0 = 5.29177e-11 m      - Bohr radius\n  μB = 9.27401e-24 J T⁻¹  - Bohr magneton\n  μN = 5.05078e-27 J T⁻¹  - nuclear magneton\n  μ₀ = 1.25664e-6 N A⁻²   - magnetic permitivity of vacuum\n  ε₀ = 8.85419e-12 F m⁻¹  - electric permitivity of vacuum\n  KJ = 4.83598e14 Hz V⁻¹  - Josephson constant\n  RK = 25812.8 Ω          - Von Klitzing constant\n   R = 8.31446 J mol⁻¹K⁻¹ - Molar gas constant\n   u = 1.66054e-27 kg     - unified atomic mass unit\n"
    codata = castCodata(2022)
    @test listCodata(codata; msg=false) == "∆νCs = 9192631770 Hz      - ¹³³Cs hyperfine transition frequency\n   c = 299792458 m s⁻¹    - speed of light in vacuum\n   h = 6.62607e-34 J Hz⁻¹ - Planck constant\n   ħ = 1.05457e-34 J s    - Planck constant (reduced)\n   e = 1.60218e-19 C      - elementary charge\n  kB = 1.38065e-23 J K⁻¹  - Boltzmann constant\n  NA = 6.02214e23 mol⁻¹   - Avogadro constant\n Kcd = 683 lm W⁻¹         - Luminous efficacy\n  mₑ = 9.10938e-31 kg     - electron mass\n  mₚ = 1.67262e-27 kg     - proton mass\n  R∞ = 1.09737e7 m⁻¹      - Rydberg constant\n  Ry = 3.28984e15 Hz      - Rydberg frequency\n  Eₕ = 4.35974e-18 J      - Hartree atomic unit\n   α = 0.00729735         - fine-structure constant\n  a0 = 5.29177e-11 m      - Bohr radius\n  μB = 9.27401e-24 J T⁻¹  - Bohr magneton\n  μN = 5.05078e-27 J T⁻¹  - nuclear magneton\n  μ₀ = 1.25664e-6 N A⁻²   - magnetic permitivity of vacuum\n  ε₀ = 8.85419e-12 F m⁻¹  - electric permitivity of vacuum\n  KJ = 4.83598e14 Hz V⁻¹  - Josephson constant\n  RK = 25812.8 Ω          - Von Klitzing constant\n   R = 8.31446 J mol⁻¹K⁻¹ - Molar gas constant\n   u = 1.66054e-27 kg     - unified atomic mass unit\n"
    @test convertUnit(1, codata) == Value(6.5796839204999635, "PHz") 
    @test convertUnit(1, codata; unitIn="Hz", unitOut="J") == Value(6.62607015e-34, "J")
    @test convertUnit(1, codata; unitIn="Hartree", unitOut="eV") == Value(27.211386245981025, "eV")
    @test convertUnit(1, codata; unitOut="cm-1") == Value(219474.63136314, "cm-1")
    @test convertUnit(1, codata; unitOut="K") == Value(315775.02480398543, "K")
    @test convertUnit(1, codata; unitIn="J", unitOut="eV") == Value(6.241509074460763e18, "eV")
    @test convertUnit(1, codata; unitIn="J", unitOut="cm-1") == Value(5.0341165675427096e22, "cm-1")
    @test convertUnit(1, codata; unitIn="J", unitOut="K") == Value(7.24297051603992e22, "K")
    @test convertUnit(1, codata; unitIn="eV", unitOut="cm-1") == Value(8065.543937349212, "cm-1")
    @test convertUnit(1, codata; unitIn="eV", unitOut="K") == Value(11604.518121550082, "K")
    @test convertUnit(1, codata; unitIn="cm-1", unitOut="K") == Value(1.4387768775039338, "K")
    @test strValue(Value(1, "Hz")) == "1 Hz"
    @test castNamedValue(Value(1.602176634e-19, "C"); name="e") == NamedValue(Value(1.602176634e-19, "C"), "e", " ")
    @test calibrationReport(1.1, 1.0, codata; unitIn="Hartree", msg=false) == "\ncalibration report (Float64):\nEcal = 1 Hartree \nE = 1.1000000000000001 Hartree \nabsolute accuracy: ΔE = 0.1 Hartree (657.968 THz)\nrelative accuracy: ΔE/E = 0.0909091\n"
    @test bohrformula(2, 4) == -1 // 8
    @test get(dictAtomicNumber, "Rb", nothing) == 37
    @test get(dictElement, 37, nothing) == ("rubidium", "Rb", 85.468)
    @test get(dictIsotope, (10, 21), nothing) ==  ("²¹Ne" , "neon"     ,  10,  21,  11, 2.9695,  20.99384669, 3//2, 1, 1.0e100, -0.6617, 0.1016, 0.27)
    @test get(dictAtomicOrbital, "3s", nothing) == (3, 0)
    Z = get(dictAtomicNumber, "Ta", nothing);
    @test get(dictConfiguration, (Z,0), nothing) == "[Yb]5d³"
    @test get(dictCoreConfiguration, "[Yb]", nothing) == "1s²2s²2p⁶3s²3p⁶3d¹⁰4s²4p⁶5s²4d¹⁰5p⁶6s²4f¹⁴"
    
    @test extractCore("[Ar]4s¹") == "1s²2s²2p⁶3s²3p⁶"
    @test extractValence("[Ar]4s¹") == "4s¹"
    @test extractValence("1s²2s²2p⁶3s²3p⁶4s¹") == "1s²2s²2p⁶3s²3p⁶4s¹"

    @test collectConfig("[Be]2p⁵") == ["1s↓0", "1s↑0", "2s↓0", "2s↑0", "2p↓-1", "2p↓0", "2p↓1", "2p↑-1", "2p↑0"]
    @test collectConfig("2p⁵") ==  ["2p↓-1", "2p↓0", "2p↓1", "2p↑-1", "2p↑0"]
    @test collectConfig("[Be]") ==  ["1s↓0", "1s↑0", "2s↓0", "2s↑0"]
    @test collectConfig("1s↑1s↓02p↑-1") == ["1s↑", "1s↓0", "2p↑-1"]

    @test collectSpinorbit("[Be]2p⁵") == collectSpinorbit("1s²2s²2p⁵")

    @test castSpinorbit("1s"; msg=true) == Spinorbit("1s↓", 1, 0, 0, 0, -1//2)
    @test castSpinorbit("1s↓"; msg=false) == Spinorbit("1s↓", 1, 0, 0, 0, -1//2)
    @test castSpinorbit("1s↓0"; msg=false) == Spinorbit("1s↓", 1, 0, 0, 0, -1//2)
    @test castSpinorbit("2p"; msg=false) == Spinorbit("2p↓", 2, 0, 1, 0, -1//2)
    @test castSpinorbit("2p↓"; msg=false) == Spinorbit("2p↓", 2, 0, 1, 0, -1//2)
    @test castSpinorbit("2p↓0"; msg=false) == Spinorbit("2p↓", 2, 0, 1, 0, -1//2)
    @test castSpinorbit("2p↑-1"; msg=false) == Spinorbit("2p↑", 2, 0, 1, -1, 1//2)
    
    @test castElement(Z=1, msg=false) == Element("hydrogen", "H", 1.008)
    @test castElement("Rb"; msg=false) == castElement(Z=37, msg=false)
    @test listElement("H") == listElement(1)
    @test listElement(1; fmt=String) == "H, hydrogen, Z=1, weight=1.008"
    @test listElement(1; fmt=Info, msg=false) == "Element: hydrogen\n  symbol: H\n  atomic number: Z = 1\n  atomic weight (relative atomic mass): 1.008"
    @test listElements(1, 3) == listElements(1:3)
    @test listElements(1:2; fmt=String) == ["H, hydrogen, Z=1, weight=1.008", "He, helium, Z=2, weight=4.0026"]
    @test castIsotope(Z=1, A=1, msg=false) == Isotope("¹H", "hydrogen", 1, 1, 0, 0.8783, 1.007825032, 1 // 2, 1, 1.0e100, 2.792847351, 0.0, 99.9855)
    @test castIsotope("Rb"; A=87, msg=false) == castIsotope(Z=37, A=87, msg=false)
    @test_throws DomainError listIsotope(1,3; fmt=1, msg=false)
    @test listIsotope(2, 3; fmt=Latex) == "2 & helium & \$^{3}\$He & 3\\, & 1 & 1.9661 & 3.016029322 & 1/2\$^+\$ & -2.12762531 & 0.0 & 0.0002 \\\\\n"
    @test listIsotope(2, 3; fmt=String) == "³He, helium, Z=2, A=3, N=1, R=1.9661, M=3.016029322, I=1/2⁺, μI=-2.12762531, Q=0.0, RA=0.0002%, (stable)"
    @test listIsotope(1,3; fmt=Info, msg=false) == "Isotope: tritium-3\n  symbol: ³T\n  element: tritium\n  atomic number: Z = 1\n  atomic mass number: A = 3\n  neutron number: N = 2\n  rms nuclear charge radius: R = 1.7591 fm\n  atomic mass: M = 3.016049281 amu\n  nuclear spin: I = 1/2 ħ\n  parity of nuclear state: π = even\n  nuclear magnetic dipole moment: μI = 2.97896246 μN\n  nuclear electric quadrupole moment: Q = 0.0 barn\n  relative abundance: RA = trace\n  lifetime: 12.33 years"
    @test listIsotopes(1,3) == listIsotopes(1:3)
    @test_throws DomainError listAtom(1,3,0; fmt=Latex)
    @test listAtom("H", 3, 0) == listAtom(1, 3, 0)
    @test listAtom(1, 3, 1; fmt=String) == "tritium ion, ³Tᐩ, Z=1, A=3, Q=1, Zc=2, config=[bare nucleus]"
    @test listAtom(1, 3, 1; fmt=Info, msg=false) == "Atom: tritium ion\n  symbol: ³Tᐩ\n  atomic charge: Z = 1\n  Rydberg charge: Zc = 2\n  electron configuration: [bare nucleus]"
    @test listAtoms(2:2, 0; fmt=String) == Any["helium, neutral atom, ³He, Z=2, A=3, Q=0, Zc=1, config=[He]", "helium, neutral atom, ⁴He, Z=2, A=4, Q=0, Zc=1, config=[He]"]
    @test castAtom("Rb"; A=87, Q=0, msg=false) == castAtom(Z=37, A=87, Q=0, msg=false)
#   ---------------------------------------------------------------------------------
    @test castAtom(Z=1, A=1, Q=0, msg=true) == Atom(1, 1, 0, 1, Element("hydrogen", "H", 1.008), Isotope("¹H", "hydrogen", 1, 1, 0, 0.8783, 1.007825032, 1 // 2, 1, 1.0e100, 2.792847351, 0.0, 99.9855), "1s¹")
    @test castOrbit("2s"; msg=true) == Orbit("2s", 2, 1, 0, 0)

    @test castTerm(1; ℓ=0, S=1 // 2, L=0, J=1 // 2, msg=true) == Term("1s ²S₁⸝₂", 1, 0, 0, 1 // 2, 0, 1 // 2)
    shell = castShell("3s"; msg=false)
    @test shell.spinorbit[1] == Spinorbit("3s↓", 3, 2, 0, 0, -1//2)
    shells = castShells("1s2s",msg=true) 
    @test shells.shell[1].spinorbit[1] == Spinorbit("1s↓", 1, 0, 0, 0, -1//2)

    atom = castAtom(Z=1, A=1, Q=0);
    spinorbit = castSpinorbit(n=2, ℓ=0);
    grid = autoGrid(atom, Float64; N=5000, msg=true);
    castDef(grid, atom, spinorbit, codata, msg=true);
#   ---------------------------------------------------------------------------------
    @test lc_eltype(([1 // 2, 1 // 3]; (1 // 4, 1 // 1, 1 // 6))) == Rational{Int}
    @test lc_eltype(([1//2, 1//3]; (1//4, big(1)//big(5), 1//6))) == Rational
    @test lc_primitivetype(([1 // 2, 1 // 3]; (1 // 4, 1 // 1, 1 // 6))) == Int64
    @test primitivetype(Rational{UInt16}) == UInt16
    @test conditionalType(47, 46) == BigInt 
#   ---------------------------------------------------------------------------------------- 
    atom = castAtom(Z=1, A=1, Q=0)
    spinorbit = castSpinorbit(n=10, ℓ=6)
    grid = autoGrid(atom, Float64; nmax=spinorbit.n, N=5000);
    def = castDef(grid, atom, spinorbit, codata)
    Ecal = convert(grid.T, bohrformula(atom.Z, spinorbit.n))
    E = 0 
    scr = zeros(grid.T,grid.N)
    def, Z = test_adams_moulton(E, scr, grid, def; test=1, msg=false)
    @test real(Z[1:9]) ≈ [0.0, 2.1136553063413117e-29, 7.545007043950154e-28, 5.9779440932204965e-27, 2.594588600478603e-26, 8.107873707845301e-26, 2.0586700404193268e-25, 4.529698568573715e-25, 8.974980054314587e-25] 
    def, Z = test_adams_moulton(E, scr, grid, def; test=2, msg=false)
    @test real(Z[1:9]) ≈ [0.0, 6.283484877494698e-40, 2.2429834003188549e-38, 1.7771261565424394e-37, 7.713205736723988e-37, 2.410312678635134e-36, 7.186046646055796e-36, 2.1536202482561519e-35, 5.474176687868545e-35]
    def, Z = test_adams_moulton(E, scr, grid, def; test=3, msg=false)
    @test real(Z[1:9]) ≈ [0.0, 6.283484877494698e-40, 2.2429834003188549e-38, 1.7771261565424394e-37, 7.713205736723988e-37, 2.410312678635134e-36, 7.186046646055796e-36, 2.1536202482561519e-35, 5.474176687868545e-35] 
    def, Z = test_adams_moulton(E, scr, grid, def; test=4, msg=false)
    @test real(Z[1:9]) ≈ [0.0, 6.283484877494698e-40, 2.2429834003188549e-38, 1.7771261565424394e-37, 7.713205736723988e-37, 2.410312678635134e-36, 7.186046646055796e-36, 2.1536202482561519e-35, 5.474176687868545e-35]
    def, Z = test_adams_moulton(E, scr, grid, def; test=5, msg=false)
    @test real(Z[1:9]) ≈ [0.0, 1.0991926025814348e-40, 3.9237315111140236e-39, 3.1087906841661287e-38, 1.3492988132050849e-37, 4.216446634181015e-37, 1.2570809780159811e-36, 3.767405336060273e-36, 9.576174110134831e-36]
#   ----------------------------------------------------------------------------------------
    atom = castAtom(Z=1, A=1, Q=0; msg=false);
    spinorbit = castSpinorbit(n=2, ℓ=0; msg=false);
    grid = autoGrid(atom, Float64; N=5000);
    RH2s_example = [RH2s(atom.Z, grid.r[n]) for n=1:grid.N];
    ZH2s_example = reduce_wavefunction(RH2s_example, grid);
    ZH2s_generic = hydrogenic_reduced_wavefunction(atom, spinorbit, grid);
    @test ZH2s_example ≈ ZH2s_generic 
    RH2s_generic = restore_wavefunction(ZH2s_generic, atom, spinorbit, grid);  
    @test RH2s_example ≈ RH2s_generic
#   ---------------------------------------------------------------------------------------- 
    println("--- H9f ---" * repeat('-', 39))
    atom = castAtom(Z=1, A=1, Q=0; msg=false);
    spinorbit = castSpinorbit(n=9, ℓ=3; msg=false);
    grid = autoGrid(atom, Float64; nmax=spinorbit.n, msg=true);
    ZH9f_generic = hydrogenic_reduced_wavefunction(atom, spinorbit, grid);
    scr = zeros(grid.T, grid.N);
    scr[1] = 0.1;
    def = castDef(grid, atom, spinorbit, codata);
    def, adams, init, Z = adams_moulton_nodes(0, scr, grid, def; imax=25, msg=true);
    def, adams, init, Z = adams_moulton_iterate!(Z, init, grid, def, adams; imax=25, ϵ=1e-15, msg=true);
    @test isapprox(ZH9f_generic, Z; rtol=1e-4)
#   ---------------------------------------------------------------------------------------- 
    println("--- ¹H:[n=30, ℓ=29] ---" * repeat('-', 39))
    atom = castAtom(Z=1, A=1, Q=0; msg=false);
    spinorbit = castSpinorbit(n=30, ℓ=29; msg=false);
    grid = autoGrid(atom, Float64; nmax=spinorbit.n, msg=true);
    ZH3029_generic = hydrogenic_reduced_wavefunction(atom, spinorbit, grid);
    scr = zeros(grid.T, grid.N);
    def = castDef(grid, atom, spinorbit, codata);
    def, adams, init, Z = adams_moulton_nodes(0, scr, grid, def; imax=25, msg=false);
    def, adams, init, Z = adams_moulton_iterate!(Z, init, grid, def, adams; imax=25, ϵ=1e-16, msg=false);
    @test isapprox(ZH3029_generic, Z; rtol=1e-6)
#   ---------------------------------------------------------------------------------------- 
    println("--- H2p ---" * repeat('-', 39))
    atom = castAtom(Z=1, A=1, Q=0; msg=false);
    spinorbit = castSpinorbit(n=2, ℓ=1; msg=false);
    grid = autoGrid(atom, Float64; p=5, rmax=60.0, nmax=spinorbit.n, N=5000);
    RH2p_example = [RH2p(atom.Z, grid.r[n]) for n=1:grid.N];
    ZH2p_example = reduce_wavefunction(RH2p_example, grid);
    ZH2p_generic = hydrogenic_reduced_wavefunction(atom, spinorbit, grid);
    @test ZH2p_example ≈ ZH2p_generic 
    RH2p_generic = restore_wavefunction(ZH2p_generic, atom, spinorbit, grid);  
    @test RH2p_example ≈ RH2p_generic 
    scr = zeros(grid.T, grid.N);
    def = castDef(grid, atom, spinorbit, codata);
    def, adams, init, Z = adams_moulton_nodes(0, scr, grid, def; imax=25, msg=false);
    def, adams, init, Z = adams_moulton_iterate!(Z, init, grid, def, adams; imax=25, ϵ=1.0e-15, msg=false);
    @test isapprox(ZH2p_generic, Z; rtol=1e-5)
#   ---------------------------------------------------------------------------------------- 
    println("--- H2s ---" * repeat('-', 39))
    atom = castAtom(Z=1, A=1, Q=0; msg=false);
    spinorbit = castSpinorbit(n=2, ℓ=0; msg=false);
    grid = autoGrid(atom, Float64; rmax=60.0, nmax=spinorbit.n, N=10000);
    RH2s_example = [RH2s(atom.Z, grid.r[n]) for n=1:grid.N];
    ZH2s_example = reduce_wavefunction(RH2s_example, grid);
    ZH2s_generic = hydrogenic_reduced_wavefunction(atom, spinorbit, grid);
    @test ZH2s_example ≈ ZH2s_generic 
    RH2s_generic = restore_wavefunction(ZH2s_generic, atom, spinorbit, grid);  
    @test isapprox(RH2s_example, RH2s_generic; rtol=1e-6)
    E=0;
    scr = zeros(grid.T, grid.N);
    def = castDef(grid, atom, spinorbit, codata);
    def, adams, init, Z = adams_moulton_nodes(E, scr, grid, def; imax=25, msg=false);
    def, adams, init, Z = adams_moulton_iterate!(Z, init, grid, def, adams; imax=25, ϵ=1.0e-20, msg=false);
    @test isapprox(ZH2s_generic, Z; rtol=1e-5)
#   ---------------------------------------------------------------------------------------- 
    println("--- H1s ---" * repeat('-', 39))
    atom = castAtom(Z=1, A=1, Q=0; msg=false);
    spinorbit = castSpinorbit(n=1, ℓ=0; msg=false);
    grid = autoGrid(atom, BigFloat; rmax=25.0, nmax=spinorbit.n, N=7500);
    RZH1s_example = [RH1s(atom.Z, grid.r[n]) for n=1:grid.N];
    ZH1s_example = reduce_wavefunction(RZH1s_example, grid);
    ZH1s_generic = hydrogenic_reduced_wavefunction(atom, spinorbit, grid); 
    @test isapprox(ZH1s_example, ZH1s_generic; rtol=1e-6)
    RZH1s_generic = restore_wavefunction(ZH1s_generic, atom, spinorbit, grid);  
    @test isapprox(RZH1s_example, RZH1s_generic; rtol=1e-6)
    scr = zeros(grid.T, grid.N);
    def = castDef(grid, atom, spinorbit, codata);
    def, adams, init, Z = adams_moulton_nodes(0, scr, grid, def; imax=25, msg=false);
    def, adams, init, Z = adams_moulton_iterate!(Z, init, grid, def, adams; imax=5, ϵ=1e-15, msg=false);
    def, adams, init, Z = adams_moulton_iterate!(Z, init, grid, def, adams; imax=50, ϵ=1e-25, msg=false);
    @test isapprox(ZH1s_generic, Z; rtol=1e-5)
#   ---------------------------------------------------------------------------------------- 
    @test a_direct(2, 1, 1, 2, 2) == 2 // 35
    @test b_exchange(1, 1, 1, 2, 2) == 2 // 5
    @test a_direct(6, 3, 2, 3, -1) == -250 // 20449
    @test b_exchange(6, 3, 2, 3, -1) == 1050 // 20449
    l = 2; l′=3; wl=2(2l+1);

    @test [A_direct(k,l) for k=0:2l] == [1,0,0,0,0]
    @test [sum([2*a_direct(k, l, ml, l, 0) for ml=-l:l])//wl for k=0:4] == [1,0,0,0,0]
    @test [B_exchange(k, l, l′) for k=0:5] == [0, 3//70, 0, 2//105, 0, 5//231]
    @test [sum([b_exchange(k, l, ml, l′, 0) for ml=-l:l])//wl for k=0:(l+l′)] == [0, 3//70, 0, 2//105, 0, 5//231]
    @test autoRmax(atom, 1; rmax=84.0) == 84.0 #63.0
    @test autoNtot(1) == 500
    @test autoPrecision(100.0, 0) == Float64
#   ---------------------------------------------------------------------------------------- 
    atom = castAtom(Z=2, A=4, Q=0; msg=false);
    spinorbit1 = castSpinorbit(n=1, ℓ=0; msg=false);
    spinorbit2 = castSpinorbit(n=2, ℓ=0; msg=false);
    grid = autoGrid(atom, Float64; nmax=spinorbit2.n,);
    Z1 = hydrogenic_reduced_wavefunction(atom, spinorbit1, grid);
    Z2 = hydrogenic_reduced_wavefunction(atom, spinorbit2, grid);
    P1 = real(Z1);
    P2 = real(Z2);
    r = grid.r
    U1s = [(1-exp(-2*atom.Z*r[n])*(1+atom.Z*r[n]))/r[n] for n=1:grid.N];
    U1s = CamiDiff.regularize!(U1s);
    @test isapprox(UFk(0, P1, grid), U1s; rtol=1e-6)
    @test isapprox(UF(spinorbit1, spinorbit1, P1, grid), U1s; rtol=1e-6)
    X1s2s = [(2/27)*sqrt(2)*atom.Z*exp(-1.5*atom.Z*r[n])*(2+3atom.Z*r[n]) for n=1:grid.N];
    @test isapprox(UGk(0, P1, P2, grid), X1s2s; rtol=1e-6)
    @test isapprox(UG(spinorbit1, spinorbit1, P1, P2, grid), X1s2s; rtol=1e-6)
    @test isapprox(Fk(0, P1, P2, grid), 34/81; rtol=1e-6)
    @test isapprox(Gk(0, P1, P2, grid), 32/729; rtol=1e-6)
    @test isapprox(𝒥(spinorbit1, spinorbit2, P1, P2, grid), 34/91; rtol=1e6)
    @test isapprox(𝒥(spinorbit1, spinorbit2, P2, P1, grid), 34/91; rtol=1e6)
    @test isapprox(𝒥(spinorbit2, spinorbit1, P1, P2, grid), 34/91; rtol=1e6)
    @test isapprox(𝒦(spinorbit1, spinorbit2, P1, P2, grid), 32/729; rtol=1e6)
    @test isapprox(𝒦(spinorbit1, spinorbit2, P2, P1, grid), 32/729; rtol=1e6)
    @test isapprox(𝒦(spinorbit2, spinorbit1, P1, P2, grid), 32/729; rtol=1e6)
#   ----------------------------------------------------------------------------------------    
    atom = castAtom(Z=2, A=4, Q=0; msg=false)
    spinorbit = castSpinorbit(n=1, ℓ=0; msg=false)
    grid = autoGrid(atom, Float64; nmax=spinorbit.n, msg=false)
    Z = hydrogenic_reduced_wavefunction(atom, spinorbit, grid);
    P = real(Z);
    @test isapprox(Fk(0, P, grid)/atom.Z, 5/8; rtol=1e-10)
    spinorbit = castSpinorbit(n=2, ℓ=1; msg=false)
    grid = autoGrid(atom, Float64; nmax=spinorbit.n, msg=false)
    Z = hydrogenic_reduced_wavefunction(atom, spinorbit, grid);
    P = real(Z);
    @test isapprox(Fk(0, P, grid)/atom.Z, 93/512; rtol=1e-10)
    @test isapprox(Fk(2, P, grid)/atom.Z, 45/512; rtol=1e-10)
    spinorbit = castSpinorbit(n=3, ℓ=2; msg=false)
    grid = autoGrid(atom, Float64; nmax=spinorbit.n, msg=false)
    Z = hydrogenic_reduced_wavefunction(atom, spinorbit, grid);
    P = real(Z);
    @test isapprox(Fk(0, P, grid)/atom.Z, 3965/46080; rtol=1e-10)
    @test isapprox(Fk(2, P, grid)/atom.Z, 2093/46080; rtol=1e-10)
    @test isapprox(Fk(4, P, grid)/atom.Z, 1365/46080; rtol=1e-9)
    spinorbit = castSpinorbit(n=4, ℓ=3; msg=false)
    grid = autoGrid(atom, Float64; nmax=spinorbit.n, msg=false)
    Z = hydrogenic_reduced_wavefunction(atom, spinorbit, grid);
    P = real(Z);
    @test isapprox(Fk(0, P, grid)/atom.Z, 184331/3670016; rtol=1e-10)
    @test isapprox(Fk(2, P, grid)/atom.Z, 103275/3670016; rtol=1e-10)
    @test isapprox(Fk(4, P, grid)/atom.Z, 69003/3670016; rtol=1e-9)
    @test isapprox(Fk(6, P, grid)/atom.Z, 51051/3670016; rtol=1e-9)
#   ----------------------------------------------------------------------------------------    
    f = [-exp(-x^2) for x=-1.0:0.01:1.0];
    f0 = -0.606530659712633;
    @test getNmin(f, 1:201) == 101
    @test getNmax(f, 1:201) == 201
    @test getNcut(f0, f, 1:101) == 30
    @test getNcut(f0, f, 101:201) == 172
    @test_throws DomainError getNcut(1.0, f, 1:201)
    Nlcut = getNcut(f0, f, 1, 101);
    Nucut = getNcut(f0, f, 101, 201);
    ΔNlcut = getΔNcut(f0, f, Nlcut, fwd; ϵ = 1e-6, k=7)
    ΔNucut = getΔNcut(f0, f, Nucut, bwd; ϵ = 1e-6, k=7)
    @test 0.01*(30.0-101.0+getΔNcut(f0, f, Nlcut, fwd; ϵ = 1e-6, k=7)) ≈ -sqrt(1//2)
    @test 0.01*(172.0-101.0+getΔNcut(f0, f, Nucut, bwd; ϵ = 1.0e-6, k=7)) ≈ sqrt(1//2) 
    polynomfwd = lagrange_polynom(f, Nlcut, Nlcut+7, fwd)
    polynombwd = lagrange_polynom(f, Nucut-7, Nucut, bwd)
    @test CamiMath.polynomial(polynomfwd, ΔNlcut) ≈ f0
    @test CamiMath.polynomial(polynombwd, ΔNucut) ≈ f0
#   ----------------------------------------------------------------------------------------    
  
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
#   ------------------------------------------------------------------------------------------------------------
    @test latent_heat_vaporization("Yb", 763) == 24170.448513975916
    @test latent_heat_vaporization("Li", 623) == 18473.64020109123
    @test latent_heat_vaporization("Li", 400) == 19134.482122780522
    @test_throws DomainError latent_heat_vaporization(1, 100)
    @test svp("Yb", 763) == 2.2859295292626745
    @test svp("Li", 623) == 0.0015230367024569058
    @test svp("Li", 400) == 7.901690229445235e-11
    @test_throws DomainError svp("H", 400)
    @test melting_point("Li") == 453.65
    @test latexIsotopeTable(1:10) == latexIsotopeTable(1,10)
    @test latexIsotopeTable(11:22; continuation=true) == latexIsotopeTable(11,22; continuation=true) 
    @test silvera_goldman_triplet(10) == -7.71843646003074e-6
    @test silvera_goldman_singlet(10) == -8.696045600341206e-6
    @test silvera_goldman_exchange(10) == 9.776091403104656e-7
    grid = CamiDiff.castGrid(3,2000,Float64; rmax=20, msg=false); 
    @test silvera_goldman_potential(grid; S=1)[700] ≈ -1.501057905105486e-6
    @test silvera_goldman_potential(grid; S=0)[700] ≈ -0.00020633920967786176
    @test rotbarrier(grid; ℓ=0)[700] ≈ 0.0
    @test rotbarrier(grid; ℓ=1)[700] ≈ 0.040892272017453914

end
