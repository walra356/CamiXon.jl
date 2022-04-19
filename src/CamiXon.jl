module CamiXon

using Printf
using Dates           # used in fits_private_sector
using IntervalSets
#using LaTeXStrings
using LinearAlgebra

export sup
export sub
export frac

export myconvert
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
export mendeleev
export Atom
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

export f_diff_weight
export f_diff_weights
export f_diff_weights_array
export f_diff_expansion_coeffs_lagrange
export f_diff_expansion_weights
export summation_range
export f_diff_function_sequences
export lagrange_interpolation
export lagrange_extrapolation
export f_diff_expansion_coeffs_differentiation
export create_lagrange_differentiation_weights
export create_lagrange_differentiation_matrix
export lagrange_differentiation
export f_diff_expansion_coeffs_adams_moulton
export create_adams_moulton_weights
export f_diff_expansion_coeffs_adams_bashford
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
export solve_adams_moulton

export Grid
export gridname
export gridfunction
export castGrid
export autoRmax
export autoNtot
export autoPrecision
export autoSteps
export autoGrid
export grid_lagrange_derivative
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

include("codata.jl")
include("strings.jl")
include("finite_differences.jl")
include("finite_difference_adams.jl")
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
