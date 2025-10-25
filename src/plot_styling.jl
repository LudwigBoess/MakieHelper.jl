"""
    get_theme(;fontsize::Integer=15)

Get the default theme
"""
function get_theme(darkmode::Bool=false; fontsize::Integer=15)

    if darkmode
        main_theme = Theme(
            backgroundcolor = :black,
            textcolor = :white,
            linecolor = :white,
            fontsize=15,
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
                yminortickcolor = :white,),
            Colorbar = (
                tickcolor = :white,
                spinecolor = :white,
                topspinecolor = :white,
                bottomspinecolor = :white,
                leftspinecolor = :white,
                rightspinecolor = :white,
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
                    ytickalign=1),
                Colorbar = (tickalign=1)
                )
    end

    return merge(main_theme, theme_latexfonts())
end

function set_dark_theme!(; fontsize::Integer=15)
    my_theme = get_theme(true; fontsize)
    set_theme!(merged_theme)
end

function set_light_theme!(; fontsize::Integer=15)
    my_theme = get_theme(false; fontsize)
    set_theme!(merged_theme)
end