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


# ==================================== ticks125(x) ================================================================

"""
    ticks(x)

Tickvalues for x according to 1-2-5 scheme
#### Examples:
```
ticks.([5,10,21.3,50,100.1])
5-element Vector{StepRange{Int64, Int64}}:
 -5:1:5
 -10:2:10
 -20:5:20
 -50:10:50
 -100:20:100
```
"""
function ticks(x)

    max = typeof(x) <: IntervalSets.ClosedInterval ? x.right : x[end]

    Δx = step125(max)
    xm = (round(Int, max) ÷ Δx) * Δx
    return (-xm:Δx:xm)

end

# ==================================== edges(itr) ============================================================

"""
    edges(itr)

Iterators defining the pixel edges of a heatmap (1-2-5 scheme).
#### Examples:
```
edges(1:100)
 0.5:1.0:99.5

edges(-200:100)
 -200.5:1.0:99.5

edges(-200:10:100)
 -200.5:10.0:99.5

edges(-200..100)
 -200.5..99.5

edges([0, 1, 5, 7, 12])
 [0, 1, 5, 7, 12]
```
"""
function edges(itr)

    typeof(itr) <: Base.OneTo && return (0.5):(itr[end] - 0.5)
    typeof(itr) <: UnitRange  && return (itr[1]-0.5):(itr[end] - 0.5)
    typeof(itr) <: StepRange  && return (itr[1]-0.5):itr.step:(itr[end] - 0.5)
    typeof(itr) <: LinRange   && return (itr[1]-0.5):itr.len:(itr[end] - 0.5)
    typeof(itr) <: Vector     && return itr
    typeof(itr) <: IntervalSets.ClosedInterval && return (itr.left-0.5)..(itr.right-0.5)

    return error("Edges: use different iterator type")

end
