module CamiXon

export indices
export indices_cnt
export partitions_cnt
export permutations_cnt


"""
    partitions_cnt(n::Int,k::Int)

The number of integer partitions of n in k parts
#### Example:
```
partitions_cnt(5,2)
 2
```
"""
function partitions_cnt(n::Int,k::Int)
    (n<0)|(k<0)|(k>n) ? 0 : (k==n)|(k==1) ? 1 : partitions_cnt(n-k,k) + partitions_cnt(n-1,k-1)
end


"""
    partitions_cnt(n::Int)

The total number of integer partitions of n

#### Example
```
partitions_cnt(5)
 7
```
"""
function partitions_cnt(n) 
    c = 1
    for k=2:n c += partitions_cnt(n,k) end
    c
end


"""
    indices(A::AbstractArray{T,1}, a::T...)  where T

The index (indices) of selected Array elements (default: all elements)

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

A = [1,2,3,4,2]
indices(A,2)
1-element Array{Array{Int64,1},1}:
 [2, 5]
```
"""
function indices(A::AbstractArray{T,1}, a::T...)  where T
    a == () ? a = unique(A) : false
    [findall(A .== fill(a[i],length(A))) for i in eachindex(a)]
end


"""
    indices_cnt(A::AbstractArray{T,1}, a::T...)  where T

The number of indices of selected Array elements (default: all elements)

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
function indices_cnt(A::AbstractArray{T,1}, a::T...)  where T
    a == () ? a = unique(A) : false
    [length(findall(A .== fill(a[i],length(A)))) for i in eachindex(a)]
end


"""
    permutations_cnt(A::AbstractArray{T,1}; unique = false)  where T

The number of permutations (option: unique permutations) of the elements of a 1D array

#### Examples:
```
A = collect("ahsgh")
permutations_cnt(A)
 120

permutations_cnt(A; unique=true)
 60
```
"""
function permutations_cnt(A::AbstractArray{T,1}; unique = false)  where T
    if unique
        o = factorial(length(A))
        c = indices_cnt(A)
        for i in eachindex(c) o = o ÷ c[i] end
        return o
    else
        return factorial(length(A)) # default
    end
end

end
