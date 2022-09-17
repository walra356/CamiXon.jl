module CamiXon

using Printf
using Dates           # used in fits_private_sector
#using IntervalSets
#using LaTeXStrings
using LinearAlgebra

export sup
export sub
export frac
export strRational

export fwd
export bwd
export reg
export rev
export isforward
export isregular
export Object
export Info
export Latex

#export convert
export convertUnit
export Value
export strValue
export NamedValue
export castNamedValue
export Codata
export castCodata
export listCodata
export calibrationReport

export bohrformula
export Element
export dictAtomicNumbers
export dictElements
export listElement
export listElements
export castElement
export Isotope
export dictIsotopes
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
export SpinOrbit
export createSpinOrbit
export Term
export createTerm

export FITS_HDU
export FITS_header
export FITS_data
export FITS_table
export parse_FITS_TABLE
export FITS_name
export cast_FITS_name

export fdiff_weight
export fdiff_expansion_weights
export fdiff_expansion
export fwd_diff_expansion_weights
export fdiff_interpolation_expansion_coeffs
export fdiff_interpolation
export fdiff_differentiation_expansion_coeffs
export fdiff_differentiation
export create_lagrange_differentiation_matrix
export fdiff_adams_moulton_expansion_coeffs
export create_adams_moulton_weights
export fdiff_adams_bashford_expansion_coeffs
export trapezoidal_weights
export trapezoidal_integration

export matG
export matσ
export matMinv
export OUTSCH
export OUTSCH_WKB
export Adams
export castAdams
export updateAdams!
export INSCH
export adams_moulton_inward
export adams_moulton_outward
export adams_moulton_normalized
export get_nodes
export adams_moulton_solve

export Grid
export gridname
export gridfunction
export castGrid
export findIndex
export autoRmax
export autoNtot
export autoPrecision
export autoSteps
export autoGrid
export grid_differentiation
export grid_trapezoidal_integral

export get_Na
export get_Nb
export get_Nlctp
export get_Nmin
export get_Nuctp
export Pos
export Def
export castDef
export initE

export fits_create
export fits_read
export fits_extend
export fits_info
export fits_copy
export fits_combine
export fits_add_key
export fits_edit_key
export fits_delete_key
export fits_rename_key

#export plot_matrices
#export plot!
export step125
export select125
export edges
export steps
export stepcenters
export stepedges
#export centerticks
#export edgeticks
#export centers
#export edges

export FORTRAN_format
export cast_FORTRAN_format
export cast_FORTRAN_datatype

export find_all
export find_first
export find_last

export bernoulli_numbers
export VectorRational
export normalize_VectorRational
export canonical_partitions
export faulhaber_polynom
export faulhaber_summation
export harmonic_number
export integer_partitions
export log10_characteristic_power
export log10_mantissa
export pascal_triangle
export pascal_next
export pochhammer
export polynomial
export polynom_derivative
export polynom_derivatives
export polynom_derivatives_all
export polynom_primitive
export polynom_product
export polynom_product_expansion
export polynom_power
export polynom_powers
export permutations_unique_count
export texp

include("singleton.jl")
include("codata.jl")
include("dicts.jl")
include("strings.jl")
include("latex.jl")
include("finite_differences.jl")
include("finite_difference_adams.jl")
include("element.jl")
include("isotope.jl")
include("atom.jl")
include("grid.jl")
include("grid_autoset.jl")
include("def.jl")
include("outsch.jl")
include("adams.jl")
include("insch.jl")
include("adams-moulton.jl")
include("fits_pointers.jl")
include("read_io.jl")
include("write_io.jl")
include("fits_objects.jl")
include("fits_private_sector.jl")
include("fits_public_sector.jl")
include("plot_private_sector.jl")
include("plot_public_sector.jl")
include("Header-and-Data-Input.jl")
include("FORTRAN.jl")
include("search_algorithms.jl")
include("mathematics.jl")

end
