# ............................................ FITS private sector ................................................................

using Dates

function _fits_bzero(E::DataType)

    return E ∉ [Int8, UInt16, UInt32, UInt64, UInt128] ? 0.0 : E == Int8 ? -128 : 2^(8E.size-1)

end

function _fits_eltype(nbits::Int, bzero::Number)

    E = nbits ==   8 ? UInt8 :
        nbits ==  16 ? Int16 :
        nbits ==  32 ? Int32 :
        nbits ==  64 ? Int64 : Int64

    E = nbits == -32 ? Float32 :
        nbits == -64 ? Float64 : E

    E = bzero == 0.0 ? E : E == UInt8 ? Int8 :
                           E == Int16 ? UInt16 :
                           E == Int32 ? UInt32 :
                           E == Int64 ? UInt64 : E

    return E

end

function _is_recordvalue_charstring(val::Any)

    _isascii_printable(val) || error("FitsError: string not within ASCII range 32-126")

    v = collect(strip(val))

    length(v) > 1 || error("FitsError: string value not delimited by single quotes")

    (v[1] == '\'') & (v[end] == '\'') || error("FitsError: string value not delimited by single quotes")

    return true

end

function _format_recordvalue(val::Any)

    typeof(val) <: Real   && return _format_recordvalue_numeric(val)
    typeof(val) <: String && _is_recordvalue_charstring(val) && return _format_recordvalue_charstring(val)
    typeof(val) <: DateTime && return _format_recordvalue_datetime(val)

    return error("FitsError: invalid record value type")

end

function _format_recordvalue_charstring(val::String)

    v = strip(val)

    return length(v) == 10 ? rpad(v,20) : length(v) < 21 ? rpad(v[1:end-1],19) * "'" : val

end

function _format_recordvalue_datetime(val::Dates.DateTime)

    return "'" * string(val) * "'"

end

function _format_recordvalue_numeric(val::Real)

    val = typeof(val) != Bool ? string(val) : val == true ? "T" : "F"

    return lpad(val,20)

end

function _format_recordkey(key::String)

    strkey = "FitsError: Key truncated at 8 characters (FITS standard)"

    length(Base.strip(key)) < 9 ? key = Base.strip(key) : (println(strkey); key = key[1:8])

    return rpad(key,8)

end

function _format_recordcomment(com::String)

    c = strip(com)

    return length(c) < 48 ? rpad(c,47) : com

end

function _format_longstringrecord(key::String, val::String, com::String)

    lcom = length(com)
    last = lcom > 0 ? "&'" : "'"
    lval = length(val)

    records = lval > 70 ? [key * "= " * val[1:68] * "&'"] : [key * "= " * val[1:end]]
    comment = lcom > 61 ? ("CONTINUE  '&' / " * com[1:61]) : lcom > 0 ? rpad(("CONTINUE  '' / " * com),80) : ""

    if  (lval < 70)
        records[1] = records[1][1:end-1] * last
    elseif (lval == 70)
        chr = records[1][end-1]
        records[1] = records[1][1:end-2] * last
        push!(records,rpad(("CONTINUE  '" * chr * last),80))
    elseif lval > 70
         val = val[1:end-1]
        lval = length(val)
        nrec = (length(val)-68)÷67 +1
        nrec > 1 ? [push!(records,("CONTINUE  '" * val[68i+2-i:68(i+1)-i] * "&'")) for i=1:nrec-1] : 0
         val = val[68(nrec)+2-nrec:end]; lval = length(val)
        (lcom == 0) & (lval == 1) ? (records[end] = records[end][1:end-2] * val * "'"; lval  =0) : 0
        lval > 0 ? push!(records,rpad(("CONTINUE  '" * val * last),80)) : 0
    end
    lcom > 0 ? push!(records,comment) : 0

    return records

end

function _fits_new_records(key::String, val::Any, com::String)

    isnumeric = typeof(val) <: Real ? true : false

    key = _format_recordkey(key)
    val = _format_recordvalue(val)
    com = _format_recordcomment(com::String)
    com = length(com) < 48 ? rpad(strip(com),67-length(val)) : com

    records = length(val * com) < 68 ? [(key * "= " * val * " / " * com)] : _format_longstringrecord(key, val, com)

    return records

end

function _validate_FITS_name(filename::String)

    return cast_FITS_name(filename)

end

function _isavailable(filename::String, protect::Bool)

    str = "FitsWarning: '$filename': filename in use (set ';protect=false' to overwrite)"

    available = (!Base.Filesystem.isfile(filename) | !protect)
    available || println(str)

    return  available ? true : false

end

function _isascii_printable(str::Union{Char,String})

    return isascii(str) ? !convert(Bool,sum(.!(31 .< Int.(collect(str)) .< 127))) : false

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
