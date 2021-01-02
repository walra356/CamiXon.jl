using CamiXon
using Test

export find_indices

"""
    find_indices(A::AbstractArray{T,N}, a::T...)  where {T,N}

Find all indices of selected Array elements (default: all elements).  

#### Example 1:
```
julia> A = collect("ahsgh");
julia> find_indices(A,'h')
1-element Array{Array{Int64,1},1}:
 [2, 5]
julia> find_indices(A)
4-element Array{Array{Int64,1},1}:
 [1]
 [2, 5]
 [3]
 [4]
``` 
#### Example 2:
```
julia> A = [1,2,3,4,2];
julia> find_indices(A,2)
1-element Array{Array{Int64,1},1}:
 [2, 5]
julia> find_indices(A)
4-element Array{Array{Int64,1},1}:
 [1]
 [2, 5]
 [3]
 [4]
``` 
""" 
function find_indices(A::AbstractArray{T,N}, a::T...)  where {T,N}
    if a == () # default
        Au = unique(A)
        return [findall(A .== fill(Au[i],length(A))) for i in eachindex(Au)]
    else
        return [findall(A .== fill(a[i],length(A))) for i in eachindex(a)]
    end
end

end
