"""
    Plot2Dset

Plot object to decompose the names of .fits files.
Fields:
* `.aspect::Real`: aspect ratio
* `.center::Tuple{Bool,Bool}`: centering true/false of ticks on x/y axes
* `.ticks::Tuple{AbstractRange,AbstractRange}`: ticks range on x/y axes
* `.title::String`:  title of plot
* `.labels::Tuple{String,String}`: labels on x/y axes
"""
struct Plot2Dset

    aspect::Real
    center::Tuple{Bool,Bool}
    ticks::Tuple{AbstractRange,AbstractRange}
    title::String
    labels::Tuple{String,String}

end
