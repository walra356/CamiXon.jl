# ================================== ConditionalType(n::T, nc::T [; msg=false]) ===========
@doc raw"""
    ConditionalType(n::T, nc::T [; msg=true]) where T<:Integer  

Convert type `T` to `BigInt` for `n > nc`.
#### Example:
```
ConditionalType(46, 46)
#  Int64

ConditionalType(47, 46)
#  BigInt
```
"""
function ConditionalType(n::T, nc::Int; msg=true) where {T<:Integer}

    T == Int || return BigInt

    V = n > nc ? BigInt : T
    V ≠ T ? (msg && println("Warning: output converted to BigInt")) : false

    return V

end

@doc raw"""
    bigconvert(o)

Convert type `T` to `BigInt` for `n > nc`.
#### Example:
```
julia> o = [[1//1, 1//2, 1//3],[1//1, 1//2, 1//3]]
2-element Vector{Vector{Rational{Int64}}}:
 [1//1, 1//2, 1//3]
 [1//1, 1//2, 1//3]
 
julia> bigconvert(o)
2-element Vector{Vector{Rational{Int64}}}:
 [1//1, 1//2, 1//3]
 [1//1, 1//2, 1//3]
```
"""
function bigconvert(o)

    T = typeof(o)
    W = get(dictBigConversion, T, "key not implemented")
    o = convert(W, o)

    return o

end

# ========================= find_all(A [,a...]; count=false) ===========

"""
    find_all(A [,a...]; count=false)

A: string/array of elements of the same type

default   : Array containing the index (indices) of selected elements of A (default: all elements)

count=true: The number of indices found for selected elements of A (default: all elements)
#### Examples:
```
A = [:📑,:📌,:📢,:📌,:📞]
B = [1,2,3,2,5]
str = "aβcβd";
find_all(A) == find_all(B) == find_all(str)
true

find_all(A,:📌)
1-element Array{Array{Int64,1},1}:
 [2, 4]

find_all(str)
4-element Array{Array{Int64,1},1}:
 [1]
 [2, 4]
 [3]
 [5]

find_all(A; count=true)
4-element Array{Int64,1}:
 1
 2
 1
 1

str = "📑📌📢📌📞"
find_all(str,'📌')
1-element Array{Array{Int64,1},1}:
 [2, 4]
```
"""
function find_all(A::Union{String,AbstractArray{T,1}}, a::T...; count=false)  where T

    typeof(A) == String ? A = Base.collect(A) : nothing

    a == () ? a = Base.unique(A) : nothing

    o = [Base.findall(x -> x == a[i], A) for i in eachindex(a)]

    return count ? length.(o) : o

end


"""
    find_first(A [,a...]; dict=false)
The first index of selected Array element

A: string/array of elements of the same type

default  : Array containing the first index (indices) of selected elements of A (default: all elements)

dict=true: Dict for the first index (indices) of selected elements of A (default: all elements)
#### Examples:
```
A = [:📑,:📌,:📢,:📌,:📞]
B = [1,2,3,2,5]
str = "aβcβd";

find_first(A) == find_first(B) == find_first(str)
true

find_first(A,:📌)
1-element Array{Array{Int64,1},1}:
 2

find_last(A,:📌; dict=true)
1-element Array{Pair{Symbol,Int64},1}:
 :📌 => 2

find_last(A; dict=true)
4-element Array{Pair{Symbol,Int64},1}:
 :📑 => 1
 :📌 => 2
 :📢 => 3
 :📞 => 5

find_first(str)
4-element Array{Int64,1}:
 1
 2
 3
 5
```
"""
function find_first(A::Union{String,AbstractArray{T,1}}, a::T...; dict=false)  where T

    typeof(A) == String ? A = Base.collect(A) : nothing

    a == () ? a = Base.unique(A) : nothing

    o = [Base.findfirst(x -> x == a[i], A) for i in eachindex(a)]

    return dict ? [a[i] => o[i] for i in eachindex(a)] : o

end

"""
    find_last(A [,a...]; dict=false)

The last index of selected Array element

A: string/array of elements of the same type

default  : Array containing the lasst index (indices) of selected elements of A (default: all elements)

dict=true: Dict for the lasst index (indices) of selected elements of A (default: all elements)
#### Examples:
```
A = [:📑,:📌,:📢,:📌,:📞]
B = [1,2,3,2,5]
str = "aβcβd";
find_last(A) == find_first(B) == find_first(str)
true

find_last(A,:📌)
1-element Array{Array{Int64,1},1}:
 4

find_last(A,:📌; dict=true)
1-element Array{Pair{Symbol,Int64},1}:
 :📌 => 4

find_last(A; dict=true)
4-element Array{Pair{Symbol,Int64},1}:
 :📑 => 1
 :📌 => 4
 :📢 => 3
 :📞 => 5

find_last(str)
4-element Array{Int64,1}:
 1
 4
 3
 5
```
"""
function find_last(A::Union{String,AbstractArray{T,1}}, a::T...; dict=false)  where T

    typeof(A) == String ? A = Base.collect(A) : nothing

    a == () ? a = Base.unique(A) : nothing

    o = [Base.findlast(x -> x == a[i], A) for i ∈ eachindex(a)]

    return dict ? [a[i] => o[i] for i ∈ eachindex(a)] : o

end
