module CamiXon

export indices
export indices_cnt
export partitions_cnt
export indices_cnt
export permutations_cnt


"""
    partitions_cnt(n::Int,k::Int)
    the number of partitions of n in k parts
"""
function p(n::Int,k::Int)
    (n<0)|(k<0)|(k>n) ? 0 : (k==n)|(k==1) ? 1 : partitions_cnt(n-k,k) + partitions_cnt(n-1,k-1)
end


"""
    partitions_cnt(n::Int,k::Int)
    the total number of partitions of n
"""
function p(n) #total number of integer partions of n
    c = 1
    for k=2:n c += partitions_cnt(n,k) end
    c
end


"""
    function indices(A::AbstractArray{T,N}, a::T...)  where {T,N}

Find the index (indices) of selected Array elements (default: all elements).

#### Examples:
```
A = collect("ahsgh")
indices(A,'h')
1-element Array{Array{Int64,1},1}:
 [2, 5]

indices(A)
4-element Array{Array{Int64,1},1}:
 [1]
 [2, 5]
 [3]
 [4]
```
"""
function indices(A::AbstractArray{T,N}, a::T...)  where {T,N}
    a == () ? a = unique(A) : false
    [findall(A .== fill(a[i],length(A))) for i in eachindex(a)]
end


"""
    indices_cnt(A::AbstractArray{T,N}, a::T...)  where {T,N}

Count the number of indices of selected Array elements (default: all elements).

#### Examples:
```
A = collect("ahsgh")
indices_cnt(A,'h')
1-element Array{Array{Int64,1},1}:
 2

indices_cnt(A)
4-element Array{Array{Int64,1},1}:
 1
 2
 1
 1
```
"""
function indices_cnt(A::AbstractArray{T,N}, a::T...)  where {T,N}
    a == () ? a = unique(A) : false
    [length(findall(A .== fill(a[i],length(A)))) for i in eachindex(a)]
end


"""
    permutation_cnt(A::AbstractArray{T,N}; unique = false)  where {T,N}

Count the number of permutations (option: unique permutations).

#### Examples:
```
A = collect("ahsgh")
permutation_cnt(A)
 120

permutation_cnt(A; unique=true)
 60
```
"""
function permutation_cnt(A::AbstractArray{T,N}; unique = false)  where {T,N}
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
