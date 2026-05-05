# Helper to map standard Matplotlib string colors to Makie
_mcolor(c) = c == "k" ? :black : (c == "w" ? :white : Symbol(c))

"""
    plot_image_grid(Nrows, Ncols, files, im_cmap, cb_labels, vmin_arr, vmax_arr, plot_name; ...)

Creates an `image_grid` plot with `Ncols` and `Nrows` with colorbars automatically handled by CairoMakie layouting.
"""
function plot_image_grid(Nrows, Ncols, files, im_cmap, cb_labels, vmin_arr, vmax_arr, plot_name;
                                map_arr = nothing, par_arr = nothing,
                                contour_arr = nothing, contour_par_arr = nothing,
                                log_map = trues(Nrows*Ncols),
                                colorbar_bottom=false,
                                smooth_file = falses(Nrows*Ncols),
                                smooth_sizes = 0.0,
                                annotate_smoothing = falses(Nrows*Ncols),
                                smooth_col=falses(Nrows*Ncols),
                                streamline_files = nothing,
                                streamlines = falses(Nrows*Ncols),
                                contour_files = nothing,
                                contours = falses(Nrows*Ncols),
                                contour_levels = nothing,
                                contour_color = "white",
                                smooth_contour_file = falses(Nrows*Ncols),
                                alpha_contours = ones(Nrows*Ncols),
                                cutoffs = nothing,
                                mask_bad = falses(Nrows*Ncols),
                                bad_colors = ["k" for _ = 1:Nrows*Ncols],
                                annotate_time = falses(Nrows*Ncols),
                                time_labels = nothing,
                                annotate_text = falses(Nrows*Ncols),
                                text_labels = nothing,
                                annotate_scale = trues(Nrows*Ncols),
                                scale_label = L"1 \: h^{-1} c" * "Mpc",
                                scale_kpc = 1000.0,
                                r_circles = nothing,
                                circle_color = "w",
                                circle_alpha = 0.5,
                                circle_lines = [":", "--", "-.", "_"],
                                smooth_contour_col=falses(Nrows*Ncols),
                                shift_colorbar_labels_inward = trues(Nrows*Ncols),
                                upscale = Ncols,
                                read_mode = 1,
                                ticks_color = "k",
                                annotation_color = "w",
                                colorbar_location="top",
                                colorbar_mode="edge",
                                grid_direction="column",
                                dpi=400,
                                overplotting_functions=nothing,
                                colorbar_size = 25
                            )

    # Added figure_padding = 30 to prevent edge tick labels from clipping
    fig = Figure(size = (600 * upscale * Ncols, 650 * upscale * Nrows), figure_padding = 30)

    if grid_direction == "column"
        Ncols, Nrows = Nrows, Ncols
    end

    Nfile = 1
    Ncontour = 1

    for row = 1:Nrows
        for col = 1:Ncols

            @info "Plotting Row $row, Column $col"

            ax = Axis(fig[row, col], aspect=DataAspect())
            hidedecorations!(ax)
            hidespines!(ax)

            # read map and parameters
            map, par = read_map_par(read_mode, Nfile, files, map_arr, par_arr)

            println("Maximum value of map: ", maximum(map))
            println("Minimum value of map: ", minimum(map))

            if smooth_file[Nfile]
                map = smooth_map!(map, smooth_sizes[Nfile], par)
            end

            if smooth_col[col]
                map = smooth_map!(map, smooth_sizes[col], par)
            end

            if !isnothing(cutoffs)
                map[map .< cutoffs[col]] .= cutoffs[col]
                map[isnan.(map)] .= cutoffs[col]
                map[isinf.(map)] .= cutoffs[col]
            end

            selected = colorbar_mode == "single" ? 1 : col

            # Image Attributes
            im_kwargs = Dict{Symbol, Any}()
            im_kwargs[:colormap] = Symbol(im_cmap[selected])
            im_kwargs[:colorrange] = (vmin_arr[selected], vmax_arr[selected])
            im_kwargs[:rasterize] = true

            if log_map[selected]
                im_kwargs[:colorscale] = log10
                # ensure strict positivity for log
                map = max.(map, vmin_arr[selected]) 
            end

            if mask_bad[selected]
                im_kwargs[:nan_color] = _mcolor(bad_colors[selected])
            end

            im = image!(ax, map; im_kwargs...)

            # Contours
            if contours[selected]
                cmap_contour, cpar = read_map_par(read_mode, Nfile, contour_files, contour_arr, contour_par_arr)

                if smooth_contour_file[Nfile]
                    cmap_contour = smooth_map!(cmap_contour, smooth_sizes[Nfile], cpar)
                end

                if !isnothing(contour_levels)
                    cmap_contour[cmap_contour .< contour_levels[1]] .= contour_levels[1]
                    cmap_contour[isnan.(cmap_contour)] .= contour_levels[1]
                    cmap_contour[isinf.(cmap_contour)] .= contour_levels[1]
                    contour!(ax, cmap_contour, levels=contour_levels, color=_mcolor(contour_color), linewidth=1.2, linestyle=:dash, alpha=alpha_contours[selected])
                else
                    contour!(ax, cmap_contour, color=_mcolor(contour_color), linewidth=1.2, linestyle=:dash, alpha=alpha_contours[selected])
                end
                Ncontour += 1
            end

            # Streamlines
            if streamlines[selected]
                @info "streamlines"
                vx, cpar, snap_num, units = read_fits_image(streamline_files[Nfile])
                vy, cpar, snap_num, units = read_fits_image(streamline_files[Nfile+1])

                x_grid = 1:size(vx, 1)
                y_grid = 1:size(vx, 2)

                streamplot!(ax, x_grid, y_grid, vx, vy, density=1.5, color=(RGBf(209/255, 209/255, 224/255), 0.7), linewidth=0.5)
                Ncontour += 1
            end

            # Additional overplotting
            if !isnothing(overplotting_functions)
                overplotting_functions[Nfile](ax)
            end

            map_x_pixels = size(map, 2)
            pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / map_x_pixels

            # Annotations
            if annotate_scale[Nfile]
                length_x = scale_kpc / pixelSideLength
                
                # Define padding from the edges (adjust the divisor to move it closer/further from the edge)
                padding_x = map_x_pixels / 10
                padding_y = map_x_pixels / 10 
                
                # Calculate X coordinates for the lower-right
                end_x = map_x_pixels - padding_x
                start_x = end_x - length_x
                
                # Draw the scale line
                lines!(ax, [start_x, end_x], [padding_y, padding_y], color=_mcolor(annotation_color), linewidth=2)
                
                # Place the label centered above the line
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
                add_circle_contours!(ax, r_circles, ["$r" for r in r_circles], circle_alpha, circle_lines, par)
            end

            # Smoothing beam (Ellipse)
            if smooth_col[Nfile] || smooth_contour_col[col] || annotate_smoothing[col]
                smooth_pixel = smooth_sizes[Nfile] ./ pixelSideLength
                scatter!(ax, [Point2f(0.1 * par.Npixels[1], 0.1 * par.Npixels[2])], 
                         marker=:circle, markersize=Vec2f(smooth_pixel[1], smooth_pixel[2]), 
                         color=:white, markerspace=:data)
            end

            # Logic to force vmin and vmax ticks to be displayed
            vmin_val = vmin_arr[selected]
            vmax_val = vmax_arr[selected]

            if log_map[selected]
                minortick_intervals = IntervalsBetween(9)
                
                ## Get integer powers of 10 and explicitly bookend them with vmin and vmax
                p_min = floor(log10(vmin_val))
                p_max = ceil(log10(vmax_val))
                Nticks = Int(p_max - p_min) + 1
                if Nticks > 10
                    Nticks /= 2
                end
                cb_ticks = 10.0.^LinRange(p_min, p_max, Nticks)
                cbtickformat = values -> [iszero(value) ? "0" : rich("$(value < 0 ? "-" : "")10", superscript("$(@sprintf("%i", log10(abs(value))))")) for value in values]

            else
                minortick_intervals = IntervalsBetween(10)
                cb_ticks = Makie.automatic
                cbtickformat = values -> [@sprintf("%g", value) for value in values]
            end

            # Natively handled Colorbars with custom explicitly defined ticks
            if colorbar_mode == "single"
                if col == Ncols && row == Nrows
                    if colorbar_location == "top" || (colorbar_location == "single" && !colorbar_bottom)
                        Colorbar(fig[0, 1:Ncols], im, label=cb_labels[selected], vertical=false,
                                ticks=cb_ticks, minorticks=minortick_intervals,
                                size=colorbar_size, alignmode=Outside(),
                                tickformat=cbtickformat)
                    elseif colorbar_location == "bottom" || colorbar_bottom
                        Colorbar(fig[Nrows+1, 1:Ncols], im, label=cb_labels[selected], vertical=false, flipaxis=false,
                                ticks=cb_ticks, minorticks=minortick_intervals,
                                size=colorbar_size, alignmode=Outside(),
                                tickformat=cbtickformat)
                    elseif colorbar_location == "right"
                        Colorbar(fig[1:Nrows, Ncols+1], im, label=cb_labels[selected], vertical=true,
                                ticks=cb_ticks, minorticks=minortick_intervals,
                                size=colorbar_size, alignmode=Outside(),
                                tickformat=cbtickformat )
                    end
                end
            else
                if colorbar_location == "top" && row == 1
                    Colorbar(fig[0, col], im, label=cb_labels[selected], vertical=false,
                            ticks=cb_ticks, minorticks=minortick_intervals,
                            size=colorbar_size, alignmode=Outside(),
                            tickformat=cbtickformat)
                elseif colorbar_location == "bottom" && row == Nrows
                    Colorbar(fig[Nrows+1, col], im, label=cb_labels[selected], vertical=false, flipaxis=false,
                            ticks=cb_ticks, minorticks=minortick_intervals,
                            size=colorbar_size, alignmode=Outside(),
                            tickformat=cbtickformat)
                elseif colorbar_location == "right" && col == Ncols
                    Colorbar(fig[row, Ncols+1], im, label=cb_labels[selected], vertical=true,
                            ticks=cb_ticks, minorticks=minortick_intervals,
                            size=colorbar_size, alignmode=Outside(),
                            tickformat=cbtickformat)
                end
            end

            Nfile += 1
        end
    end

    # Remove layout paddings to keep images flush
    rowgap!(fig.layout, 0)
    colgap!(fig.layout, 0)

    # Force square aspect ratio on grid cells to prevent DataAspect() whitespace
    for r in 1:Nrows
        rowsize!(fig.layout, r, Auto(1.0))
    end
    for c in 1:Ncols
        colsize!(fig.layout, c, Aspect(1, 1.0)) 
    end

    # Shrink wrap the figure to the constrained layout
    resize_to_layout!(fig)

    @info "saving $plot_name"
    save(plot_name, fig)
end