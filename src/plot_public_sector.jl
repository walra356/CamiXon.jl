# ==================================== step125(x) ============================================================

"""
    step125(x)

Step used for deviding the number x in steps according to 1-2-5 scheme
#### Examples:
```
step125.([5,10,21.3,50,100.1])
5-element Vector{Int64}:
  1
  2
  5
 10
 20
```
"""
function step125(x::Real)

    m = log10_mantissa(x)
    p = log10_characteristic_power(x)
    v = 10^m
    d = v > 7.9 ? 2.0 : v > 3.9 ? 1.0 : v > 1.49 ? 0.5 : 0.2

    return max(1,round(Int, d *= 10^p))

end


# ==================================== select125(x) ===============================================================

"""
    select125(x)

Select elements of the collection x by index according to 1-2-5 scheme
#### Examples:
```@docs
x1 = [1,2,4,6,8,10,13,16,18,20,40,60,80,100]
x2 = string.(x1)
x3 = Base.OneTo(100)
x4 = (1:100)
println(select125(x1)); println(select125(x2)); println(select125(x3)); println(select125(x4))
 [2, 6, 10, 16, 20, 60, 100]
 ["2", "6", "10", "16", "20", "60", "100"]
 [20, 40, 60, 80, 100]
 [20, 40, 60, 80, 100]
```
"""
select125(x) = (n = length(x); return [x[i] for i=step125(n):step125(n):n])

# ==================================== edges(x) ===============================================================

"""
    edges(x)

Heatmap range transformation from `center` (default) to `edge` format.
For ClosedInterval ranging the heatmap dimension has to be supplied manually.
#### Examples:
```@docs
edges(Base.OneTo(10))
 0.5:1.0:9.5

edges(UnitRange(-8:0))
 -8.5:1.0:-0.5

edges(range(-21.1, 0, length=6))
 [-23.21, -18.99, -14.77, -10.549999999999999, -6.33, -2.11]

edges(LinRange(-21.1,0,6))
 [-23.21, -18.990000000000002, -14.77, -10.55, -6.33, -2.1100000000000003]

edges((-21.1)..0; dim=6)
 (-23.21)..(-2.1100000000000003)
```
"""
function edges(x; dim = 0)

    T = typeof(x)
    E = eltype(x)

    strErr = "RangeType $T not implemented"

    T <: Base.OneTo   ? Δx = 1 :
    T <: UnitRange    ? Δx = 1 :
    T <: StepRange    ? Δx = x.step :
    T <: StepRangeLen ? Δx = x.step :
    T <: LinRange     ? Δx = (x.stop-x.start)/(x.len-1) :
    T <: IntervalSets.ClosedInterval ? Δx = (x.right-x.left)/(dim-1) : error(strErr)

    return dim > 0 ? (x.left-0.5Δx)..(x.right-0.5Δx) : Δx == 1 ? x .- 0.5Δx : x .- E[0.5Δx]

end

# =================================== steps(x) ===============================================================

"""
    steps(x)

Heatmap range transformation for steplength specification vector x
#### Examples:
```@docs
steps([4,2,6])
 [0, 4, 6, 12]
```
"""
function steps(x::Vector{T} where T<:Real)

    sum(x .< 0) == 0 || error("Error: $x - nagative step length not allowed")

    return (s = append!(eltype(x)[0],x); [Base.sum(s[1:i]) for i ∈ Base.eachindex(s)])

end

# =================================== stepcenters(x) ===============================================================
"""
    stepcenters(x)

Stepcenter positions for steplength specification vector x
#### Examples:
```@docs
stepcenters([4,2,6])
 [2.0, 5.0, 9.0]
```
"""
function stepcenters(x::Vector{T} where T<:Real)

    s = append!(eltype(x)[0],x)

    return [Base.sum(s[1:i]) + 0.5x[i] for i ∈ Base.eachindex(x)]

end

# =================================== stepedges(x) ===============================================================

"""
    stepedges(x)

Stepedge positions for steplength specification vector x
#### Examples:
```@docs
stepedges([4,2,6])
 [0, 4, 6, 12]
```
"""
function stepedges(x::Vector{T} where T<:Real)

    return steps(x)

end
