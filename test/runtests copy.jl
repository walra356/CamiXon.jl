# SPDX-License-Identifier: MIT

using CamiXon
using CamiDiff
using CamiMath

# using IntervalSets
#using BenchmarkTools
using LinearAlgebra
using Test

@testset "CamiXon.jl" begin 
#   ---------------------------------------------------------------------------------------- 
    codata = castCodata(2022)
    atom = castAtom(Z=1, A=1, Q=0)
    orbit = castOrbit(n=1, ℓ=0)
    grid = autoGrid(atom, orbit, BigFloat; Ntot=5000);
    def = castDef(grid, atom, orbit, codata)
    Ecal = convert(grid.T, bohrformula(atom.Z, orbit.n))
    E = 0 
    scr = zeros(grid.T,grid.N)
    def, Z = test_adams_moulton(E, scr, grid, def; test=1, msg=false)
    @test real(Z[1:9]) ≈  [0.0, 5.1428402317309706e-27, 1.8358054664173181e-25, 1.4545109486041696e-24, 6.312941767665052e-24, 1.9727333231630704e-23, 5.008945116783607e-23, 1.1021150943733992e-22, 2.1836812980487107e-22]
    def, Z = test_adams_moulton(E, scr, grid, def; test=2, msg=false)
    @test real(Z[1:9]) ≈ [0.0, 1.316529120127285e-36, 4.699526422219377e-35, 3.7234406147142423e-34, 1.6160664722811144e-33, 5.050051623547137e-33, 1.282252959498499e-32, 3.182655514351617e-32, 8.201804186634829e-32]
    def, Z = test_adams_moulton(E, scr, grid, def; test=3, msg=false)
    @test real(Z[1:9]) ≈ [0.0, 1.316529120127285e-36, 4.699526422219377e-35, 3.7234406147142423e-34, 1.6160664722811144e-33, 5.050051623547137e-33, 1.282252959498499e-32, 3.182655514351617e-32, 8.201804186634829e-32]
    def, Z = test_adams_moulton(E, scr, grid, def; test=4, msg=false)
    @test real(Z[1:9]) ≈ [0.0, 1.316529120127285e-36, 4.699526422219377e-35, 3.7234406147142423e-34, 1.6160664722811144e-33, 5.050051623547137e-33, 1.282252959498499e-32, 3.182655514351617e-32, 8.201804186634829e-32]
    def, Z = test_adams_moulton(E, scr, grid, def; test=5, msg=false)
    @test real(Z[1:9]) ≈ [0.0, 2.303207641828953e-37, 8.221607105497416e-36, 6.513989509686655e-35, 2.8272345759443254e-34, 8.834834955917978e-34, 2.2432430623252777e-33, 5.567910644661615e-33, 1.434868229698988e-32]
    println(Z[1:9])
end