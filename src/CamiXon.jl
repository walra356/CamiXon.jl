module CamiXon

import Documenter
import DocumenterInterLinks
import CamiMath
import CamiDiff
import Printf
import Dates           
#using IntervalSets
#using LaTeXStrings
import LinearAlgebra

fwd = CamiMath.fwd
bwd = CamiMath.bwd
reg = CamiMath.reg
rev = CamiMath.rev

sup(i) = CamiMath.sup(i)
sub(i) = CamiMath.sub(i)
undosup(i) = CamiMath.undosup(i)
undosub(i) = CamiMath.undosub(i)
#frac(i) = CamiMath.frac(i)
strRational(i) = CamiMath.strRational(i)

struct Object
end

struct Info
end

struct Latex
end

export Object
export Info
export Latex

export primitivetype
export lc_primitivetype
export lc_eltype
export conditionalType
export find_all
export find_first
export find_last

export convertUnit
export Value
export strValue
export NamedValue
export castNamedValue
export Codata
export castCodata
export listCodata
export calibrationReport

export dictAntoineCoefficient
export dictAtomicNumber
export dictAtomicOrbital
export dictElement
export dictIsotope
export dictMeltingPoint
export dictCoreConfiguration
export dictConfiguration

export extractCore
export extractValence
export collectConfig
export collectSpinorbit

export bohrformula
export Element
export listElement
export listElements
export castElement
export Isotope
export listIsotope
export listIsotopes
export castIsotope
export latexIsotopeTable
export Atom
export listAtom
export listAtoms
export castAtom
export Orbit
export castOrbit
export Spinorbit
export castSpinorbit
export Shell
export castShell
export Shells
export castShells
export Term
export castTerm
export a_direct
export A_direct
export b_exchange
export B_exchange

export UFk
export UGk
export UF
export UG
export Fk
export Gk
export 𝒥
export 𝒦

export svp
export latent_heat_vaporization
export melting_point

export Ein
export castEin
export inE!

export matG
export matσ 
export matMinv
export OUTSCH!
export OUTSCH_WKB!
export OUTSCH_WJ!
export Adams
export castAdams
export updateAdams!
export INSCH!
export INSCH_WKB!
export INSCH_WJ!
export adams_moulton_outward!
export adams_moulton_inward!
export adams_moulton_normalize!
export adams_moulton_patch
export adams_moulton_solve!
export adams_moulton_solve_refine!
export adams_moulton_nodes
export adams_moulton_iterate!
export adams_moulton_report_nodes
export adams_moulton_report_iterate
export test_adams_moulton

export demo_hydrogen
export hydrogenic_reduced_wavefunction
export RH1s
export RH2s
export RH2p
export silvera_goldman_triplet
export silvera_goldman_singlet
export silvera_goldman_exchange
export silvera_goldman_potential
export rotbarrier
export restore_wavefunction
export reduce_wavefunction

export autoRmax
export autoNtot
export autoPrecision
export autoGrid

export getNmin
export getNmax
export getNcut
export getΔNcut

export Pos
export castPos
export updatePos!
export listPos
export Def
export castDef



include("julia.jl")
include("codata.jl")
include("dicts.jl")
include("latex.jl")
include("element.jl")
include("isotope.jl")
include("atom.jl")
include("orbit.jl")
include("thermal-properties.jl")
include("grid_autoset.jl")
include("pos.jl")
include("def.jl")
include("inE.jl")
include("hydrogen.jl")
include("adams.jl")
include("outsch.jl")
include("insch.jl")
include("Coulomb_Integrals.jl")
include("adams-moulton.jl")

end
