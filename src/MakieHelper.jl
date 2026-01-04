module MakieHelper

    using Makie

    include("binning.jl")
    include("plot_styling.jl")
    include("symlog10.jl")

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
            Symlog10

end
