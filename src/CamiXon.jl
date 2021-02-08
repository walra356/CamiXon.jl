module CamiXon

export decompose_filnam
export fits_info
export fits_copy
export fits_combine
export find_all
export find_first
export find_last
export canonical_partitions
export integer_partitions

include("file_manipulation.jl")
include("search_algorithms.jl")
include("mathematics.jl") 

end
