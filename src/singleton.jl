# =========================== singletons =======================================

struct fwd
    end

struct bwd
    end

struct reg
    end

struct rev
    end

struct Object
    end

struct Info
    end

struct Latex
    end

# ============================= isforward(notation) ===========================

"""
    function isforward(val)

Boolean status of `val`, with options: `fwd` (forward) and `bwd` (backward).
#### Example:
```
isforward(fwd)
  true
```
"""
function isforward(val)

    strErr = "Error: invalid value of val (options: fwd, bwd)"

    return val === fwd ? true : val === bwd ? false : error(strErr)

end
# ============================= End ============================================

# ============================= isregular(ordering) ============================
"""
    function isregular(val)

Boolean status of `val`, with options: `reg` (regular) and `rev` (reversed).
#### Example:
```
isregular(reg)
  true
```
"""
function isregular(val)

    val === reg && return true

    val === rev && return false

    error("Error: invalid value of val (options: reg, rev)")

end
# ============================= End ===========================
