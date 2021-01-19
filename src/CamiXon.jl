module CamiXon

export indices
export indices_cnt
export partitions_cnt
export canonical_partitions
export next_partitions
export all_partitions
export permutations_cnt


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

#### Example:
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
    canonical_partitions(A,i; trailer=false)

The canonical partition in integers of element i of array A

trailer: array A included as last element in the output container

#### Examples:
```
canonical_partitions([6],1; trailer=true)
6-element Array{Array{Int64,1},1}:
 [1, 1, 1, 1, 1, 1]
 [2, 2, 2]
 [3, 3]
 [4, 2]
 [5, 1]
 [6]

canonical_partitions([6],1)
5-element Array{Array{Int64,1},1}:
 [1, 1, 1, 1, 1, 1]
 [2, 2, 2]
 [3, 3]
 [4, 2]
 [5, 1]

canonical_partitions([4,4,4,4],3; trailer=true)
4-element Array{Array{Int64,1},1}:
 [4, 4, 1, 1, 1, 1, 1, 1, 1, 1]
 [4, 4, 2, 2, 2, 2]
 [4, 4, 3, 3, 2]
 [4, 4, 4, 4]

canonical_partitions([4,4,1,1],3)
Error: integer partition of 1 at index 3 in [4, 4, 1, 1] coincides with the trailer
```
"""
function canonical_partitions(A::Array{Int,1},i::Int; trailer=false)
    n = Base.sum(A)                                     # n: partition constant
    A1 = A[1:i-1]                                       # i: partition index (array A[i] to be partitioned)
    ni = n-Base.sum(A1)                                 # ni: sub-partition constant at index i
    m = trailer ? A[i] : A[i]-1                         # m: depth of partition (with or without trailer)
    o = [Base.fill(k,Base.cld(ni,k)) for k=1:m]         # init partitions starting at element i and fill folder o
    for k in eachindex(o)
        o[k][Base.cld(ni,k)]=((ni%k)≠0 ? ni%k : k)      # adjust last value of the partitions
        o[k] = prepend!(o[k],A1)                        # complete partition with A1
    end
    o = o==[] ? println("Error: integer partition of $(A[i]) at index $i in $A coincides with the trailer") : o
end


"""
    next_partitions(o,i)

The canonical partitions of  o[1][:] o[2][:] ... at index i

#### Examples:
```
o = [[4,4,4,4]]
1-element Array{Array{Int64,1},1}:
 [4, 4, 4, 4]

o = next_partitions(o,2)
3-element Array{Array{Int64,1},1}:
 [4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 [4, 2, 2, 2, 2, 2, 2]
 [4, 3, 3, 3, 3]

o = next_partitions(o,3)
3-element Array{Array{Int64,1},1}:
 [4, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 [4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 [4, 3, 2, 2, 2, 2, 1]

o = next_partitions(o,4)
1-element Array{Array{Int64,1},1}:
 [4, 3, 2, 1, 1, 1, 1, 1, 1, 1]

o = next_partitions(o,5)
```
"""
function next_partitions(o::Array{Array{Int,1},1},i)
    c::Array{Array{Int,1},1} = []
    for p ∈ eachindex(o)
        if i<=length(o[p]) && o[p][i]>1
            c = append!(c,canonical_partitions(o[p],i))
        end
    end
    c==[] ? nothing : c
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
