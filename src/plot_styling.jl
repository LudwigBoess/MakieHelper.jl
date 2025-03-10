"""
    get_theme(;fontsize::Integer=15)

Get the default theme
"""
function get_theme(;fontsize::Integer=15)
    main_theme = Theme(
            fontsize=fontsize,
            Axis = ( xticksmirrored=true,
            yticksmirrored=true,
            xminorticksvisible=true,
            yminorticksvisible=true,
            xminortickalign=1,
            yminortickalign=1,
            xtickalign=1,
            ytickalign=1)
            )

    return merge(main_theme, theme_latexfonts())
end