using CairoMakie
using GeometryBasics
using LaTeXStrings
using Printf

# Helper to map standard Matplotlib string colors to Makie
_mcolor(c) = c == "k" ? :black : (c == "w" ? :white : Symbol(c))

"""
    propaganda_plot_double_row(Ncols, files, im_cmap, cb_labels, vmin_arr, vmax_arr, plot_name; ...)

Creates a 2-row image grid plot with `Ncols` columns, with colorbars exclusively on the top (for row 1) and bottom (for row 2).
"""
function propaganda_plot_double_row(Ncols, files, im_cmap, cb_labels, vmin_arr, vmax_arr, plot_name;
                                map_arr = nothing, par_arr = nothing,
                                contour_arr = nothing, contour_par_arr = nothing,
                                log_map = trues(2Ncols),
                                smooth_file = falses(2Ncols),
                                smooth_sizes = 0.0,
                                annotate_smoothing = falses(2Ncols),
                                streamline_files = nothing,
                                streamlines = falses(2Ncols),
                                contour_files = nothing,
                                contours = falses(2Ncols),
                                contour_level_values=nothing,
                                contour_color = "white",
                                smooth_contour_file = falses(2Ncols),
                                alpha_contours = ones(2Ncols),
                                cutoffs=vmin_arr,
                                mask_bad = trues(2Ncols),
                                bad_colors = ["k" for _ = 1:2Ncols],
                                annotate_time = falses(2Ncols),
                                time_labels = nothing,
                                annotate_text = falses(2Ncols),
                                text_labels = nothing,
                                annotate_scale = trues(2Ncols),
                                scale_label = L"1 \: h^{-1} c" * "Mpc",
                                scale_kpc = 1000.0,
                                r_circles = nothing,
                                shift_colorbar_labels_inward = trues(2Ncols),
                                upscale = Ncols,
                                aspect_ratio = 1.42,
                                read_mode = 1,
                                image_num = ones(Int64, 2Ncols),
                                ticks_color = "k",
                                N_ticks = zeros(Int64, 2Ncols),
                                annotation_color = "w",
                                dpi=400,
                                cb_label_offset=0.0,
                                overplotting_functions=nothing
                            )

    # Added figure padding to protect outer tick overhangs
    fig = Figure(size = (600 * upscale * Ncols, 600 * upscale * 2), figure_padding = 30)

    Nfile = 1
    Ncontour = 1
    
    # Loop across rows, then columns
    for row = 1:2
        for col = 1:Ncols

            @info "Column $col, Plot $row"

            # Row 1 goes into layout row 2. Row 2 goes into layout row 3.
            # Layout rows 1 and 4 are reserved for Colorbars
            image_row = row == 1 ? 2 : 3
            cb_row = row == 1 ? 1 : 4

            ax = Axis(fig[image_row, col], aspect=DataAspect())
            hidedecorations!(ax)
            hidespines!(ax)

            # read map and parameters
            map, par = read_map_par(read_mode, Nfile, files, map_arr, par_arr)

            println("Maximum value of map: ", maximum(map))
            println("Minimum value of map: ", minimum(map))

            if smooth_file[Nfile]
                map = smooth_map!(map, smooth_sizes[Nfile], par)
            end

            # Handle cutoffs by converting to NaNs so mask_bad can catch them
            if !isnothing(cutoffs)
                map[map .< cutoffs[Nfile]] .= NaN
            end

            im_kwargs = Dict{Symbol, Any}()
            im_kwargs[:colormap] = Symbol(im_cmap[Nfile])
            im_kwargs[:colorrange] = (vmin_arr[Nfile], vmax_arr[Nfile])

            if mask_bad[Nfile]
                im_kwargs[:nan_color] = _mcolor(bad_colors[Nfile])
                map[isnan.(map)] .= NaN
                map[isinf.(map)] .= NaN
            else
                map[isnan.(map)] .= vmin_arr[Nfile]
                map[isinf.(map)] .= vmin_arr[Nfile]
            end

            if log_map[Nfile]
                im_kwargs[:colorscale] = log10
                # ensure strict positivity for log
                map = max.(map, vmin_arr[Nfile]) 
            end

            im = image!(ax, map; im_kwargs...)

            # Contours
            if contours[Nfile]
                cmap_contour, cpar = read_map_par(read_mode, Nfile, contour_files, contour_arr, contour_par_arr)

                if smooth_contour_file[Nfile]
                    cmap_contour = smooth_map!(cmap_contour, smooth_sizes[Nfile], cpar)
                end

                if isnothing(contour_level_values)
                    contour!(ax, cmap_contour, color=_mcolor(contour_color), linewidth=1.2, linestyle=:dash, alpha=alpha_contours[col])
                else
                    contour!(ax, cmap_contour, levels=contour_level_values, color=_mcolor(contour_color), linewidth=1.2, linestyle=:dash, alpha=alpha_contours[col])
                end
                Ncontour += 1
            end

            # Streamlines
            if streamlines[Nfile]
                @info "streamlines"
                vx, cpar, snap_num, units = read_fits_image(streamline_files[Nfile])
                vy, cpar, snap_num, units = read_fits_image(streamline_files[Nfile+1])

                x_grid = 1:size(vx, 1)
                y_grid = 1:size(vx, 2)

                streamplot!(ax, x_grid, y_grid, vx, vy, density=2.0, color=:white, linewidth=0.5)
                Ncontour += 1
            end

            # Additional overplotting
            if !isnothing(overplotting_functions)
                overplotting_functions[Nfile](ax)
            end

            map_x_pixels = size(map, 2)
            pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / map_x_pixels

            # Annotate Scale (Only on the leftmost column)
            if col == 1 && annotate_scale[row]
                length_x = scale_kpc / pixelSideLength
                
                # Dynamically calculate right-aligned placement
                padding_x = map_x_pixels / 14
                padding_y = map_x_pixels / 14 
                end_x = map_x_pixels - padding_x
                start_x = end_x - length_x
                
                lines!(ax, [start_x, end_x], [padding_y, padding_y], color=_mcolor(annotation_color), linewidth=2)
                text!(ax, scale_label, position=Point2f(start_x + length_x/2, padding_y * 1.5), 
                      color=_mcolor(annotation_color), align=(:center, :bottom))
            end

            if annotate_time[Nfile]
                text!(ax, time_labels[Nfile], position=Point2f(0.075 * par.Npixels[1], par.Npixels[1] * 0.925), 
                      color=_mcolor(annotation_color), align=(:left, :top))
            end

            if annotate_text[Nfile]
                text!(ax, text_labels[Nfile], position=Point2f(par.Npixels[1] * 0.925, par.Npixels[1] * 0.925), 
                      color=_mcolor(annotation_color), align=(:right, :top))
            end

            # Circles
            if !isnothing(r_circles)
                add_circle_contours!(ax, r_circles, ["$r" for r in r_circles], 0.5, [":", "--"], par)
            end

            # Smoothing beam (Rectangle + Ellipse)
            if smooth_file[Nfile] || smooth_contour_file[Nfile] || annotate_smoothing[Nfile]
                smooth_pixel = smooth_sizes[Nfile] ./ pixelSideLength
                
                poly!(ax, Rect(0.1*par.Npixels[1] - 1.5*smooth_pixel[1], 0.1*par.Npixels[2] - 1.5*smooth_pixel[2], 
                               3*smooth_pixel[1], 3*smooth_pixel[2]), 
                      color=:gray, strokecolor=:gray, strokewidth=1)
                      
                scatter!(ax, [Point2f(0.1 * par.Npixels[1], 0.1 * par.Npixels[2])], 
                         marker=:circle, markersize=Vec2f(smooth_pixel[1], smooth_pixel[2]), 
                         color=:white, markerspace=:data)
            end

            # Ticks Processing
            vmin_val = vmin_arr[Nfile]
            vmax_val = vmax_arr[Nfile]
            n_ticks_user = N_ticks[Nfile]

            if log_map[Nfile]
                minortick_intervals = IntervalsBetween(9)
                if n_ticks_user > 0
                    cb_ticks = 10.0 .^ Float64.(collect(range(log10(vmin_val), log10(vmax_val), length=n_ticks_user)))
                else
                    p_min = ceil(log10(vmin_val))
                    p_max = floor(log10(vmax_val))
                    cb_ticks = [10.0^p for p in p_min:p_max]
                    pushfirst!(cb_ticks, vmin_val)
                    push!(cb_ticks, vmax_val)
                    cb_ticks = Float64.(unique(sort(cb_ticks)))
                end
                
                # Top row (row=1) gets flipaxis=true so ticks/labels go above the bar
                Colorbar(fig[cb_row, col], im, label=cb_labels[Nfile], vertical=false, flipaxis=(row==1),
                         ticks=cb_ticks, minorticks=minortick_intervals, alignmode=Outside())
            else
                minortick_intervals = IntervalsBetween(10)
                num_ticks = n_ticks_user > 0 ? n_ticks_user : 5
                cb_ticks = Float64.(collect(range(vmin_val, vmax_val, length=num_ticks)))
                
                # Apply linear tick format with auto-abbreviation
                cbtickformat_linear = values -> [@sprintf("%g", value) for value in values]

                Colorbar(fig[cb_row, col], im, label=cb_labels[Nfile], vertical=false, flipaxis=(row==1),
                         ticks=cb_ticks, tickformat=cbtickformat_linear, minorticks=minortick_intervals, alignmode=Outside())
            end

            Nfile += 1
        end
    end

    # Make the two image rows flush vertically (0 gap between layout row 2 and 3)
    rowgap!(fig.layout, 2, 0)
    
    # Make all columns fully flush horizontally
    colgap!(fig.layout, 0)

    # Force square aspect ratio on grid cells
    rowsize!(fig.layout, 2, Auto(1.0))
    rowsize!(fig.layout, 3, Auto(1.0))
    for c in 1:Ncols
        # Ties column width perfectly to image row height
        colsize!(fig.layout, c, Aspect(2, 1.0)) 
    end

    # Shrink wrap the figure to the constrained layout
    resize_to_layout!(fig)

    @info "saving $plot_name"
    save(plot_name, fig)
end