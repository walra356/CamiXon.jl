module CamiXon

export get_indices
export get_indices_count
export get_permutation_count


"""
    p(n::Int,k::Int)
    the number of partitions of n in k parts
"""
function p(n::Int,k::Int)
    (n<0)|(k<0)|(k>n) ? 0 : (k==n)|(k==1) ? 1 : p(n-k,k) + p(n-1,k-1)
end

"""
    p(n::Int,k::Int)
    the total number of partitions of n
"""
function p(n) #total number of integer partions of n
    c = 1
    for k=2:n c += p(n,k) end
    c
end



"""
    get_indices(A::AbstractArray{T,N}, a::T...)

Find the index (indices) of selected Array elements (default: all elements).

#### Examples:
```
A = collect("ahsgh")
get_indices(A,'h')
1-element Array{Array{Int64,1},1}:
 [2, 5]

get_indices(A)
4-element Array{Array{Int64,1},1}:
 [1]
 [2, 5]
 [3]
 [4]
```
"""
function get_indices(A::AbstractArray{T,N}, a::T...)  where {T,N}
    a == () ? a = unique(A) : false
    [findall(A .== fill(a[i],length(A))) for i in eachindex(a)]
end


"""
    get_indices_count(A::AbstractArray{T,N}, a::T...)  where {T,N}

Count the number of indices of selected Array elements (default: all elements).

#### Examples:
```
A = collect("ahsgh")
get_indices_count(A,'h')
1-element Array{Array{Int64,1},1}:
 2

get_indices_count(A)
4-element Array{Array{Int64,1},1}:
 1
 2
 1
 1
```
"""
function get_indices_count(A::AbstractArray{T,N}, a::T...)  where {T,N}
    a == () ? a = unique(A) : false
    [length(findall(A .== fill(a[i],length(A)))) for i in eachindex(a)]
end


"""
    get_permutation_count(A::AbstractArray{T,N}; unique = false)  where {T,N}

Count the number of permutations (option: unique permutations).

#### Examples:
```
A = collect("ahsgh")
permutation_count(A)
 120

permutation_count(A; unique=true)
 60
```
"""
function get_permutation_count(A::AbstractArray{T,N}; unique = false)  where {T,N}
    if unique
        o = factorial(length(A))
        c = find_count(A)
        for i in eachindex(c) o = o ÷ c[i] end
        return out
    else
        return factorial(length(A)) # default
    end
end

end
