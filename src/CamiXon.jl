module CamiXon

export decompose_filnam
export fits_combine
export fits_copy
export fits_info
export fits_key_copy
export fits_key_create
export fits_key_delete
export fits_key_edit
export fits_key_info
export fits_key_rename
export find_all
export find_first
export find_last
export canonical_partitions
export integer_partitions
export cast_FORTRAN_format
export cast_FORTRAN_datatype
export FORTRAN_format
export FITS_name
export FITS_table
export FITS_header
export FITS_data
export FITS_HDU

include("file_manipulation.jl")
include("fits_pointers.jl")
include("read_io.jl")
include("write_io.jl")
include("fits_objects.jl")
include("FORTRAN.jl")
include("search_algorithms.jl")
include("mathematics.jl") 

end
