"""
    Symlog10([lower=-upper,] upper; linscale=1)

An axis scaling which is linear for inputs in the interval `[lower, upper]` (where `lower < 0 < upper`) and logarithmic outside, thus representing both positive and negative values.

The parameter `linscale` (default: 1) controls how much space should be used for the linear region in the output, relative to decades in the logarithmic region.
Specifically, the linear region `[lower, upper]` will occupy the same space as `2 * linscale` decades in the output.

If only one argument is given, `lower` is set to `-upper`, and the linear region is symmetric around zero.

WARNING: The gradient of this transformation is discontinuous at `lower` and `upper`, which may lead to visual artifacts in the data. Other scales such as `AsinhScale` or `pseudolog10` are smooth and do not have this issue.
"""
Symlog10(upper; kwargs...) = Symlog10(-upper, upper; kwargs...)
function Symlog10(lower, upper; linscale = 1)

    lower >= 0 && throw(ArgumentError("Argument `lower` must be < 0. Got: $lower"))
    upper <= 0 && throw(ArgumentError("Argument `upper` must be > 0. Got: $upper"))
    linscale <= 0 && throw(ArgumentError("Argument `linscale` must be > 0. Got: $linscale"))

    function forward(x)
        if lower < x < upper
            x = ((x - lower) / (upper - lower) * 2 - 1) * linscale
        else
            x = sign(x) * (linscale + log10(abs(x) / (x > 0 ? upper : abs(lower))))
        end
        return x - (-lower / (upper - lower) * 2 - 1) * linscale  # Shifts so that 0 maps to 0
    end
    function inverse(x)
        x += (-lower / (upper - lower) * 2 - 1) * linscale  # Undo the shift
        if abs(x) < linscale
            x = (x / linscale + 1) / 2 * (upper - lower) + lower
        else
            x = sign(x) * exp10(abs(x) - linscale) * (x > 0 ? upper : abs(lower))
        end
        return x
    end

    return ReversibleScale(forward, inverse; limits = (-3.0f0, 3.0f0), name = :Symlog10)
end
