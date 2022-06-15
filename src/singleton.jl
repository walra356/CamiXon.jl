# =========================== singletons =======================================

struct fwd
end

struct bwd
end

struct Object
    end

struct Info
    end

struct Latex
    end

# ============================= isforward(notation) ===========================

function isforward(notation)

    notation === fwd && return true

    notation === bwd && return false

    error("Error: unknown fdiff notation type")

end
# ============================= End ===========================

