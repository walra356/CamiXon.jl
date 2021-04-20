# ............................................ FITS private sector ................................................................

function _validate_FITS_name(filename::String)
     
    return cast_FITS_name(filename)
    
end

function _isavailable(filename::String, protect::Bool)
     
    str = "FitsWarning: '$filename': filename in use (set ';protect=false' to overwrite)"
     
    available = (!Base.Filesystem.isfile(filename) | !protect)
    available || println(str)
        
    return  available ? true : false
    
end

function _rm_blanks(records::Array{String,1})               # remove blank records
    
    record_B = repeat(' ',length(records[1]))
    
    return [records[i] for i ∈ findall(records .≠ record_B)]
    
end

function _fits_parse(str::String) # routine 
    
    T = Float32    
    s = Base.strip(str)
    c = Base.collect(s)
    l = Base.length(s)    
    
    Base.length(s) == 0 && return 0 
     
    d = [Base.Unicode.isdigit(c[i]) for i ∈ Base.eachindex(c)]  # d: digits
    p = [Base.Unicode.ispunct(c[i]) for i ∈ Base.eachindex(c)]  # p: punctuation
   
    s[1] == '\'' && return str      
    s[1] == '-'  && (d[1] = true) && (p[1] = false)     # change leading sign into digit (for type parsing only)   
    s[1] == '+'  && (d[1] = true) && (p[1] = false)     # change leading sign into digit (for type parsing only)
    
    a = p .== d                                         # a: other non-digit or punctuation characters
    
    ia = [i for i ∈ Base.eachindex(a) if a[i]==1]       # ia: indices of nonzero elements of a
    ip = [i for i ∈ Base.eachindex(p) if p[i]==1]       # ip: indices of nonzero elements of p
    
    sd = Base.sum(d)                                    # sd: number of digits
    sp = Base.sum(p)                                    # sd: number of punctuation characters
    sa = Base.sum(a)                                    # sa: number of other non-digit or punctuation characters
    
    E = ['E','D','e','p']
    
    sd == l && return Base.parse(Int,s)    
    sa >= 2 && return str
    sp >= 3 && return str
    sa == 1 && s == "T" && return true
    sa == 1 && s == "F" && return false
    sa == 0 && sp == 1 && s[ip[1]] == '.' && return Base.parse(T,s)
    sp == 0 && sa == 1 && s[ia[1]] ∈ E && return Base.parse(T,s)
    sp == 1 && s[ip[1]] == '-' && s[ip[1]-1] ∈ E && return Base.parse(T,s)
    sp == 1 && s[ip[1]] == '.' && s[ia[1]] ∈ E && ip[1] < ia[1]  && return Base.parse(T,s)
    sp == 2 && s[ip[1]] == '.' && s[ip[2]] == '-' && s[ip[2]-1] ∈ E && ip[1] < ia[1]  &&  return Base.parse(T,s)
    
    return error("FitsError: $str: parsing error")

end
