module CamiXon

using Dates
using AbstractPlotting

export FITS_HDU
export FITS_header
export FITS_data
export FITS_table
export parse_FITS_TABLE
export FITS_name
export cast_FITS_name

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

export Plot2Dset
export cast_Plot2Dset

export plot_matrices

export FORTRAN_format
export cast_FORTRAN_format
export cast_FORTRAN_datatype

export find_all
export find_first
export find_last
export canonical_partitions
export integer_partitions

include("fits_pointers.jl")
include("read_io.jl")
include("write_io.jl")
include("fits_objects.jl")
include("fits_private_sector.jl")
include("fits_public_sector.jl")
include("plot_objects.jl")
include("plot_private_sector.jl")
include("plot_public_sector.jl")
include("Header-and-Data-Input.jl")
include("FORTRAN.jl")
include("search_algorithms.jl")
include("mathematics.jl")

end
