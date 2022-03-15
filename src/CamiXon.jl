module CamiXon

using Dates
using IntervalSets
#using LaTeXStrings
using LinearAlgebra

export sup
export sub
export frac

export bohrformula
export mendeleev
export Atom
export createAtom
export Orbit
export createOrbit
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

export Grid
export gridfunction
export createGrid

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

include("strings.jl")
include("adams.jl")
include("atom.jl")
include("grid.jl")
include("fits_pointers.jl")
include("read_io.jl")
include("write_io.jl")
include("finite_differences.jl")
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
