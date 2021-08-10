module CamiXon

using Dates
using IntervalSets
using LaTeXStrings
using LinearAlgebra

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
export f_diff_expansion_coeffs_interpolation
export f_diff_expansion_coeffs_extrapolation
export f_diff_expansion_coeffs_array_interpolation
export f_diff_expansion_weights
export f_diff_expansion_weights_extrapolation
export f_diff_expansion_weights_array
export summation_ranges
export f_diff_function_sequences
export lagrangian_interpolation
export f_diff_expansion_coeffs_differentiation
export f_diff_expansion_coeffs_array_differentiation
export lagrangian_differentiation

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
export canonical_partitions
export integer_partitions
export log10_characteristic_power
export log10_mantissa
export polynom
export polynom_deriv_coeffs
export polynom_multiplication_coeffs
export permutations_unique_count
export f_diff_expansion_coeffs_adams_moulton

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
