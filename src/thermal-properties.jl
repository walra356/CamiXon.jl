function svp(atomicnumber::Int, temp::Real, dictAntoineCoefficients::Dict)

    T = float(temp)
    i = atomicnumber
    d = dictAntoineCoefficients

    groupA = [3, 4, 11, 13, 19, 21, 22, 23, 26, 27, 28, 29, 30, 31, 37, 39, 40]
    append!(groupA, [45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 64, 65])
    append!(groupA, [68, 71, 78, 79, 80, 81, 82, 90, 91, 92, 93, 94, 96])
    groupB = [12, 20, 24, 25, 38, 41, 42, 44, 62, 63, 66, 67, 69, 70, 72, 73]
    append!(groupB, [74, 75, 76, 77, 95])
    groupC = [1, 2, 5, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35]
    append!(groupC, [36, 43, 51, 52, 53, 54, 61, 83, 84, 85, 86, 87, 88])
    append!(groupC, [89, 97, 98, 99, 100, 101, 102])
    groupD = (6)

    (Tmin, mp, Tmax) = d[i][3]

    if i ∈ groupA
        (A, B, C, D) = (T ≤ mp) ? d[i][1] : d[i][2]
    elseif i ∈ groupB
        (A, B, C, D) = (T ≤ mp) ? d[i][1] :
            println("Antoine coefficients for T > m.p. = $(mp) not found")
    else
        println("Antoine coefficients not found")
    end

    lp = A + B / T + C * log10(T) + D * T / 1000.0

    p = exp(lp) # saturated vapor pressure

    return p

end
function svp(element::String, temp::Float64, dictAntoineCoefficients::Dict)

    A = dictAtomicNumbers[element]
    p = svp(A, temp, dictAntoineCoefficients)

    return p

end