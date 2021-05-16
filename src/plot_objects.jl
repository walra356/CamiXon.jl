"""
    Plotset2D

Type to hold plot parameters

Fields:
* `.dims::Tuple{Int,Int}`: dimensions (size) of the matrix to be plotted
* `.gain::Real`: scaling factor to adjust the comstrast of the plot
* `.aspect::Real`: aspect ratio
* `.center::Tuple{Bool,Bool}`: centering true/false of ticks on x/y axes
* `.ticks::Tuple{AbstractRange,AbstractRange}`: ticks range on x/y axes
* `.title::String`:  title of plot
* `.labels::Tuple{String,String}`: labels on x/y axes
"""
struct Plotset2D

    dims::Tuple{Int,Int}
    gain::Real
    aspect::Real
    center::Tuple{Bool,Bool}
    range::Tuple{AbstractRange,AbstractRange}
    ticks::Tuple{AbstractRange,AbstractRange}
    labels::Tuple{String,String}

end
