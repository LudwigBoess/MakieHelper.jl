function round_to_next_N(x; N=100)
    return ceil(x / N) * N
end

function get_integer_ticks(x_lim; N_spacing=100)
    round_to_next_N(x_l[1], N=N_spacing):N_spacing:round_to_next_N(x_l[2], N=N_spacing)
end

"""
    DecadeTicks(max_labels=8)

Major tick locator for logarithmic axes. On a `log10` scale the ticks sit at
whole decades and are labelled `10ⁿ`; an axis spanning fewer than two decades
gets `1, 2, 5 × 10ⁿ` labelled `m×10ⁿ`, and one spanning more than `max_labels`
decades labels only every few decades. Every other scale falls back to Makie's
automatic ticks, so it can be the theme default for all axes and colorbars.
"""
struct DecadeTicks
    max_labels::Int
end

DecadeTicks() = DecadeTicks(8)


"""
    DecadeMinorTicks(fallback=IntervalsBetween(2))

Minor tick locator for logarithmic axes: `m × 10ⁿ`, `m = 1…9`, on a `log10`
scale, excluding the major ticks; `fallback` on every other scale.
"""
struct DecadeMinorTicks{T}
    fallback::T
end

DecadeMinorTicks() = DecadeMinorTicks(IntervalsBetween(2))


"""
    decade_tickvalues(vmin, vmax, mantissas)

All values `m × 10ⁿ` with `m` in `mantissas` between `vmin` and `vmax`, sorted.
"""
function decade_tickvalues(vmin, vmax, mantissas)

    exponents = floor(Int, log10(vmin)):ceil(Int, log10(vmax))
    values = [m * exp10(n) for n in exponents for m in mantissas]

    return sort(filter(v -> vmin <= v <= vmax, values))
end


"""
    major_decade_ticks(t::DecadeTicks, vmin, vmax)

Major tick values of a `log10` axis from `vmin` to `vmax`, see [`DecadeTicks`](@ref).
"""
function major_decade_ticks(t::DecadeTicks, vmin, vmax)

    major = decade_tickvalues(vmin, vmax, (1,))
    length(major) < 2 && return decade_tickvalues(vmin, vmax, (1, 2, 5))

    # label every `step`-th decade, aligned to exponents divisible by `step`
    step = cld(length(major), t.max_labels)
    step > 1 && filter!(v -> mod(round(Int, log10(v)), step) == 0, major)

    return major
end


"""
    decade_ticklabel(value)

Tick label `10ⁿ` of a whole decade, `m×10ⁿ` otherwise.
"""
function decade_ticklabel(value)

    n = floor(Int, log10(value) + 1.0e-9)
    m = round(Int, value / exp10(n))

    return m == 1 ? rich("10", superscript("$n")) : rich("$m×10", superscript("$n"))
end


function Makie.get_ticks(t::DecadeTicks, scale, formatter, vmin, vmax)

    scale === log10 || return Makie.get_ticks(Makie.automatic, scale, formatter, vmin, vmax)

    values = major_decade_ticks(t, vmin, vmax)
    labels = formatter isa Makie.Automatic ? decade_ticklabel.(values) :
             Makie.get_ticklabels(formatter, values)

    return values, labels
end

# unit-aware axes ask for the tick values alone
Makie.get_tickvalues(t::DecadeTicks, scale, vmin, vmax) =
    scale === log10 ? major_decade_ticks(t, vmin, vmax) :
    Makie.get_tickvalues(Makie.automatic, scale, vmin, vmax)

# categorical and date axes pick their own ticks only when they are `automatic`
Makie.get_ticks(c::Makie.CategoricalConversion, ::DecadeTicks, scale, formatter, vmin, vmax) =
    Makie.get_ticks(c, Makie.automatic, scale, formatter, vmin, vmax)
Makie.get_ticks(c::Makie.DateTimeConversion, ::DecadeTicks, scale, formatter, vmin, vmax) =
    Makie.get_ticks(c, Makie.automatic, scale, formatter, vmin, vmax)


function Makie.get_minor_tickvalues(t::DecadeMinorTicks, scale, tickvalues, vmin, vmax)

    scale === log10 || return Makie.get_minor_tickvalues(t.fallback, scale, tickvalues, vmin, vmax)

    return filter(v -> all(abs(v - major) > 1.0e-6v for major in tickvalues),
        decade_tickvalues(vmin, vmax, 1:9))
end
