module CamiXon


export fits_create
export fits_read
export fits_extend
export fits_info
export fits_copy
export fits_add_key
export fits_edit_key
export fits_delete_key
export FITS_HDU
export FITS_header
export FITS_data
export FITS_table
export parse_FITS_TABLE
export FITS_name
export cast_FITS_name
export FORTRAN_format
export cast_FORTRAN_format
export cast_FORTRAN_datatype

#export decompose_filnam
#export fits_combine
#export fits_copy
#export fits_info
#export fits_key_copy
#export fits_key_create
#export fits_key_delete
#export fits_key_edit
#export fits_key_info
#export fits_key_rename
export find_all
export find_first
export find_last
export canonical_partitions
export integer_partitions

#include("file_manipulation.jl")
include("fits_pointers.jl")
include("read_io.jl")
include("write_io.jl")
include("fits_objects.jl")
include("Header-and-Data-Input.jl")
include("FORTRAN.jl")
include("private_sector.jl")
include("public_sector.jl")
include("search_algorithms.jl")
include("mathematics.jl") 
include("run_tests.jl") 

end
