"""
    get_theme(darkmode=false; fontsize::Integer=15)

Get the default theme. Logarithmic axes and colorbars get whole-decade tick
labels and `m × 10ⁿ` minor ticks ([`DecadeTicks`](@ref),
[`DecadeMinorTicks`](@ref)); all other axes keep Makie's automatic ticks.
"""
function get_theme(darkmode::Bool=false; fontsize::Integer=15)

    if darkmode
        main_theme = Theme(
            backgroundcolor = :black,
            textcolor = :white,
            linecolor = :white,
            fontsize=fontsize,
            Axis = ( xticksmirrored=true,
                yticksmirrored=true,
                xminorticksvisible=true,
                yminorticksvisible=true,
                xminortickalign=1,
                yminortickalign=1,
                xtickalign=1,
                ytickalign=1,
                backgroundcolor = :transparent,
                bottomspinecolor = :white,
                topspinecolor = :white,
                leftspinecolor = :white,
                rightspinecolor = :white,
                xtickcolor = :white,
                ytickcolor = :white,
                xminortickcolor = :white,
                yminortickcolor = :white,
                xticks = DecadeTicks(),
                yticks = DecadeTicks(),
                xminorticks = DecadeMinorTicks(IntervalsBetween(2)),
                yminorticks = DecadeMinorTicks(IntervalsBetween(2))),
            Colorbar = (
                tickcolor = :white,
                minortickcolor = :white,
                spinecolor = :white,
                topspinecolor = :white,
                bottomspinecolor = :white,
                leftspinecolor = :white,
                rightspinecolor = :white,
                minorticksvisible=true,
                tickalign=1,
                minortickalign=1,
                ticks = DecadeTicks(),
                minorticks = DecadeMinorTicks(IntervalsBetween(5))
                )
            )
    else

        main_theme = Theme(
                fontsize=fontsize,
                Axis = ( xticksmirrored=true,
                    yticksmirrored=true,
                    xminorticksvisible=true,
                    yminorticksvisible=true,
                    xminortickalign=1,
                    yminortickalign=1,
                    xtickalign=1,
                    ytickalign=1,
                    xgridvisible = false,
                    ygridvisible = false,
                    xticks = DecadeTicks(),
                    yticks = DecadeTicks(),
                    xminorticks = DecadeMinorTicks(IntervalsBetween(2)),
                    yminorticks = DecadeMinorTicks(IntervalsBetween(2))),
                Colorbar = (minorticksvisible=true,
                tickalign=1,
                minortickalign=1,
                ticks = DecadeTicks(),
                minorticks = DecadeMinorTicks(IntervalsBetween(5)))
                )
    end

    return merge(main_theme, theme_latexfonts())
end

function set_dark_theme!(; fontsize::Integer=15)
    my_theme = get_theme(true; fontsize)
    set_theme!(my_theme)
end

function set_light_theme!(; fontsize::Integer=15)
    my_theme = get_theme(false; fontsize)
    set_theme!(my_theme)
end