var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = CamiXon","category":"page"},{"location":"#CamiXon.jl","page":"Home","title":"CamiXon.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A package for image analysis of backscattered light","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Table-of-contents","page":"Home","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Functions","page":"Home","title":"Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"find_all(A::Union{String,AbstractArray{T,1}}, a::T...; count=false)  where T\nfind_first(A::Union{String,AbstractArray{T,1}}, a::T...; dict=false)  where T\nfind_last(A::Union{String,AbstractArray{T,1}}, a::T...; dict=false)  where T\ncanonical_partitions(n::Int, m=0; header=true, reverse=true)\ninteger_partitions(n::Int, m=0; transpose=false, count=false)\npermutations_cnt(A::AbstractArray{Any,1}; unique = false)","category":"page"},{"location":"#CamiXon.find_all-Union{Tuple{T}, Tuple{Union{AbstractArray{T,1}, String},Vararg{T,N} where N}} where T","page":"Home","title":"CamiXon.find_all","text":"find_all(A [,a...]; count=false)\n\nA: string/array of elements of the same type\n\ndefault   : Array containing the index (indices) of selected elements of A (default: all elements)  count=true: The number of indices found for selected elements of A (default: all elements)  \n\nExamples:\n\nA = [:📑,:📌,:📢,:📌,:📞]\nB = [1,2,3,2,5]\nstr = \"aβcβd\";\n\nfind_all(A) == find_all(B) == find_all(str)\ntrue\n\nfind_all(A,:📌)\n1-element Array{Array{Int64,1},1}:\n [2, 4]\n\nfind_all(str)\n4-element Array{Array{Int64,1},1}:\n [1]\n [2, 4]\n [3]\n [5]\n\nfind_all(A; count=true)\n4-element Array{Int64,1}:\n 1\n 2\n 1\n 1\n\nstr = \"📑📌📢📌📞\"\nfind_all(str,'📌')\n1-element Array{Array{Int64,1},1}:\n [2, 4]\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.find_first-Union{Tuple{T}, Tuple{Union{AbstractArray{T,1}, String},Vararg{T,N} where N}} where T","page":"Home","title":"CamiXon.find_first","text":"find_first(A [,a...]; dict=false)\n\nThe first index of selected Array element\n\nA: string/array of elements of the same type\n\ndefault  : Array containing the first index (indices) of selected elements of A (default: all elements)  dict=true: Dict for the first index (indices) of selected elements of A (default: all elements) \n\nExamples:\n\nA = [:📑,:📌,:📢,:📌,:📞]\nB = [1,2,3,2,5]\nstr = \"aβcβd\";\n\nfind_first(A) == find_first(B) == find_first(str)\ntrue\n\nfind_first(A,:📌)\n1-element Array{Array{Int64,1},1}:\n 2\n\nfind_last(A,:📌; dict=true)\n1-element Array{Pair{Symbol,Int64},1}:\n :📌 => 2\n\nfind_last(A; dict=true)\n4-element Array{Pair{Symbol,Int64},1}:\n :📑 => 1\n :📌 => 2\n :📢 => 3\n :📞 => 5\n\nfind_first(str)\n4-element Array{Int64,1}:\n 1\n 2\n 3\n 5\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.find_last-Union{Tuple{T}, Tuple{Union{AbstractArray{T,1}, String},Vararg{T,N} where N}} where T","page":"Home","title":"CamiXon.find_last","text":"find_last(A [,a...]; dict=false)\n\nThe last index of selected Array element\n\nA: string/array of elements of the same type\n\ndefault  : Array containing the lasst index (indices) of selected elements of A (default: all elements)  dict=true: Dict for the lasst index (indices) of selected elements of A (default: all elements)\n\nExamples:\n\nA = [:📑,:📌,:📢,:📌,:📞]\nB = [1,2,3,2,5]\nstr = \"aβcβd\";\n\nfind_last(A) == find_first(B) == find_first(str)\ntrue\n\nfind_last(A,:📌)\n1-element Array{Array{Int64,1},1}:\n 4\n\nfind_last(A,:📌; dict=true)\n1-element Array{Pair{Symbol,Int64},1}:\n :📌 => 4\n\nfind_last(A; dict=true)\n4-element Array{Pair{Symbol,Int64},1}:\n :📑 => 1\n :📌 => 4\n :📢 => 3\n :📞 => 5\n\nfind_last(str)\n4-element Array{Int64,1}:\n 1\n 4\n 3\n 5\n\n\n\n\n\n","category":"method"},{"location":"#CamiXon.canonical_partitions","page":"Home","title":"CamiXon.canonical_partitions","text":"canonical_partitions(n; header=false, reverse=true)\n\nThe canonical partition in integers of the integer n header=true : unit patition included in output\n\nExamples:\n\ncanonical_partitions(6; header=true, reverse=false)\n6-element Array{Array{Int64,1},1}:\n [6]\n [5, 1]\n [4, 2]\n [3, 3]\n [2, 2, 2]\n [1, 1, 1, 1, 1, 1]\n\ncanonical_partitions(6; header=true)\n6-element Array{Array{Int64,1},1}:\n [1, 1, 1, 1, 1, 1]\n [2, 2, 2]\n [3, 3]\n [4, 2]\n [5, 1]\n [6]\n\ncanonical_partitions(6)\n5-element Array{Array{Int64,1},1}:\n [1, 1, 1, 1, 1, 1]\n [2, 2, 2]\n [3, 3]\n [4, 2]\n [5, 1]\n\n\n\n\n\n","category":"function"},{"location":"#CamiXon.integer_partitions","page":"Home","title":"CamiXon.integer_partitions","text":"integer_partitions(n [,m]; transpose=false, count=false)\n\ndefault              : The integer partitions of n  count=true           : The number of integer partitions of n  transpose=false/true : for m>0 restricted to partitions with maximum part/length m\n\ndefinitions:  The integer partition of the positive integer n is a nonincreasing sequence of positive integers p1, p2,... pk whose sum is n. The elements of the sequence are called the parts of the partition. \n\nExamples:\n\ninteger_partitions(7)\n15-element Array{Array{Int64,1},1}:\n [1, 1, 1, 1, 1, 1, 1]\n [2, 2, 2, 1]\n [3, 3, 1]\n [4, 3]\n [5, 2]\n [6, 1]\n [7]\n [2, 2, 1, 1, 1]\n [3, 2, 2]\n [4, 2, 1]\n [5, 1, 1]\n [2, 1, 1, 1, 1, 1]\n [3, 2, 1, 1]\n [4, 1, 1, 1]\n [3, 1, 1, 1, 1]\n\ninteger_partitions(7; count=true)\n15\n\ninteger_partitions(7,4; count=true)\n3\n\ninteger_partitions(7,4)\n3-element Array{Array{Int64,1},1}:\n [4, 3]\n [4, 2, 1]\n [4, 1, 1, 1]\n\ninteger_partitions(7,4; transpose=true)\n3-element Array{Array{Int64,1},1}:\n [2, 2, 2, 1]\n [3, 2, 1, 1]\n [4, 1, 1, 1]\n\n\n\n\n\n\n","category":"function"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}