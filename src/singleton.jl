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
    function isforward(notation)

Boolean status of notation, with options: `fwd` and `bwd`.
#### Example:
isforward(fwd)
  true
"""
function isforward(notation)

    notation === fwd && return true

    notation === bwd && return false

    error("Error: unknown fdiff notation type")

end
# ============================= End ============================================

# ============================= isregular(ordering) ============================
"""

function isregular(ordering)

Boolean status of ordering, with options: `reg` and `rev`.
#### Example:
isregular(reg)
  true
"""
function isregular(ordering)

    ordering === reg && return true

    ordering === rev && return false

    error("Error: unknown fdiff ordering type")

end
# ============================= End ===========================
