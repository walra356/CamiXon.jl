# SPDX-License-Identifier: MIT

# ====================== primitivetype(T) ==============================

@doc raw"""
    primitivetype(T::Type)

The primitive type of a Type
#### Examples:
```
julia> T = Complex{Float16}};
julia> primitivetype(T)
Float16

julia> T = String
julia> primitivetype(T)
Char
```
"""
function primitivetype(T::Type)
    
    T = eltype(T)
          
    T = T <: Union{Int, Rational{Int}}          ? Int      :  
        T <: Union{Float64, Complex{Float64}}   ? Float64  :  
        T <: Union{BigInt, Rational{BigInt}}    ? BigInt   :
        T <: Union{BigFloat, Complex{BigFloat}} ? BigFloat : 
        T <: Union{Float32, Complex{Float32}}   ? Float32  :
        T <: Union{Float16, Complex{Float16}}   ? Float16  :
        T <: Union{Int8, Rational{Int8}}        ? Int8     :  
        T <: Union{UInt8, Rational{UInt8}}      ? UInt8    :
        T <: Union{Int16, Rational{Int16}}      ? Int16    :
        T <: Union{UInt16, Rational{UInt16}}    ? UInt16   :
        T <: Union{Int32, Rational{Int32}}      ? Int32    :
        T <: Union{UInt32, Rational{UInt32}}    ? UInt32   :
        T <: Union{Int64, Rational{Int64}}      ? Int64    :
        T <: Union{UInt64, Rational{UInt64}}    ? UInt64   :
        T <: Union{Int128, Rational{Int128}}    ? Int128   :
        T <: Union{UInt128, Rational{UInt128}}  ? UInt128  : T
        
    return T
        
end

# ========================== lc_primitivetype(T) ==============================

@doc raw"""
    lc_primitivetype(o::Any)

Lowest comon primitive type of Any Type
#### Examples:
```
julia> o = ([1//2, 1//3]; (1//4, 1//1, 1//6));
julia> lc_primitivetype(o)
Int64
```
"""
function lc_primitivetype(o::Any)

    T = lc_eltype(o)

    T = T <: Union{Int,Rational{Int}} ? Int :
        T <: Union{Float64,Complex{Float64}} ? Float64 :
        T <: Union{BigInt,Rational{BigInt}} ? BigInt :
        T <: Union{BigFloat,Complex{BigFloat}} ? BigFloat :
        T <: Union{Float32,Complex{Float32}} ? Float32 :
        T <: Union{Float16,Complex{Float16}} ? Float16 :
        T <: Union{Int8,Rational{Int8}} ? Int8 :
        T <: Union{UInt8,Rational{UInt8}} ? UInt8 :
        T <: Union{Int16,Rational{Int16}} ? Int16 :
        T <: Union{UInt16,Rational{UInt16}} ? UInt16 :
        T <: Union{Int32,Rational{Int32}} ? Int32 :
        T <: Union{UInt32,Rational{UInt32}} ? UInt32 :
        T <: Union{Int64,Rational{Int64}} ? Int64 :
        T <: Union{UInt64,Rational{UInt64}} ? UInt64 :
        T <: Union{Int128,Rational{Int128}} ? Int128 :
        T <: Union{UInt128,Rational{UInt128}} ? UInt128 : T

    return T

end

# ================================ lc_eltype(o) ==============================

@doc raw"""
    lc_eltype(o)

Lowest common eltype of a collection.
#### Examples:
```
julia> o = ([1//2, 1//3]; (1//4, 1//1, 1//6));
julia> lc_eltype(o)
Rational{Int64}

julia> o = ([1//2, 1//3]; (1//4, big(1)//big(5), 1//6));
julia> lc_eltype(o)
Rational

julia> o = ([1//2, 1//3]; (1//4, [big(1)//big(5)], 1//6));
julia> lc_eltype(o)
Any

julia> o = ([1/2, 1/3]; (1/4, 1/1, 1/6));
julia> lc_eltype(o)
Float64
```
"""
function lc_eltype(o)

    T = eltype(o)
    U = eltype(T)

    while T ≠ U
        T = U
        U = eltype(T)
    end

    return T

end

# ================================== ConditionalType(n::T, nc::T [; msg=false]) ===========
@doc raw"""
    conditionalType(n::T, nc::T [; msg=true]) where T<:Integer  

Convert type `T` to `BigInt` for `n > nc`.
#### Example:
```
julia> conditionalType(46, 46)
Int64

julia> conditionalType(47, 46)
BigInt
```
"""
function conditionalType(n::T, nc::Int; msg=true) where {T<:Integer}

    T == Int || return BigInt

    V = n > nc ? BigInt : T

    return V

end

@doc raw"""
    bigconvert(o)

Convert IntBased types to `BigIntBased` types for `n > nc` in accordance with `dictBigConversion`.
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
    W = get(dictBigConversion, T, "Type not implemented")
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
