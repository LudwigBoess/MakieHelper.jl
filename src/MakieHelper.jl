module MakieHelper

    using Makie
    using Colors, ColorTypes, ColorSchemes
    using GeometryBasics
    using LaTeXStrings
    using ImageFiltering
    using SPHtoGrid
    using Printf

    include("binning.jl")
    include("plot_styling.jl")
    include("colorschemes.jl")
    include("symlog10.jl")
    include("ticks.jl")
    include("utility.jl")
    include("image_grid.jl")

    export get_theme,
            set_dark_theme!,
            set_light_theme!,
            bin_1D, bin_1D!,
            bin_1D_log, bin_1D_log!,
            bin_1D_loglog, bin_1D_loglog!,
            σ_1D_quantity,
            bin_2D, bin_2D!,
            bin_2D_quantity!,
            bin_2D_log, bin_2D_log!,
            bin_2D_quantity_log!,
            Symlog10,
            spezi,
            round_to_next_N,
            get_integer_ticks,
            plot_image_grid

end
