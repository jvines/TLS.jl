module TransitLeastSquaresMetalExt

using TransitLeastSquares
using Metal

# Placeholder: Metal-backed hot loop lives here. Wiring into the public
# API (a `backend = :metal` keyword on `tls`) is deferred to a later
# minor version. Keeping this file present and loadable ensures the
# extension hook exists in the compat matrix from v0.1 onward.

function __init__()
    @debug "TransitLeastSquaresMetalExt loaded (stub; CPU backend still used)"
end

end # module TransitLeastSquaresMetalExt
