# Helper to map standard Matplotlib string colors to Makie
_mcolor(c) = c == "k" ? :black : (c == "w" ? :white : Symbol(c))

"""
    plot_image_grid(Nrows, Ncols, files, im_cmap, cb_labels, vmin_arr, vmax_arr, plot_name; kwargs...)

Plot a grid of 2D maps with colorbars and save it to `plot_name`.

Panels sit flush in square cells without axes. Each panel can carry a scale bar, text labels,
contours, streamlines, circles and a smoothing beam.

# Panel order and indexing
Panels are filled row by row: the panel index `Nfile` runs across all columns of row 1 first.
`files`, `map_arr`, `par_arr` and every keyword indexed by `Nfile` follow this order; a matrix is
read in its column-major linear order.

Properties indexed by `selected` (colormap, colour range, log scale, NaN masking, contours,
streamlines) take one entry per colorbar:
- `colorbar_mode="edge"`: one colorbar per column, `selected = col`. With
  `colorbar_location="right"` there is one colorbar per row instead, `selected = row`.
- `colorbar_mode="single"`: one colorbar for the whole grid, `selected = 1`.

`grid_direction="column"` (the default) swaps `Nrows` and `Ncols` before the grid is built; pass
`"row"` to use them as given.

# Arguments
- `Nrows`, `Ncols`: grid shape.
- `files`: map filenames in panel order, read according to `read_mode`.
- `im_cmap`: colormap names, indexed by `selected`.
- `cb_labels`: colorbar labels, indexed by `selected`.
- `vmin_arr`, `vmax_arr`: colour range, indexed by `selected`. Log maps are clamped to `vmin_arr`.
- `plot_name`: output path; its extension sets the file format.

# Keywords
Map input:
- `read_mode = 1`: format of `files`: 1 SPHtoGrid FITS, 2 Smac1 binary, 3 Smac2, 4 Smac1 FITS.
- `map_arr`, `par_arr`: preloaded maps and their `mappingParameters` in panel order, used instead
  of `files` when `map_arr` is given.
- `cutoffs`: floor for each map; NaN and Inf pixels are set to it too. Indexed by `col`, or by
  `row` for right-hand edge colorbars.
- `log_map`: log10 colour scale, indexed by `selected`.
- `mask_bad`, `bad_colors`: colour NaN pixels with `bad_colors`, indexed by `selected`. No effect
  together with `cutoffs`, which removes NaNs.

Smoothing:
- `smooth_sizes`: Gaussian smoothing width per entry, an `(x, y)` pair in map length units; the
  kernel σ is half the width.
- `smooth_file`: smooth panel `Nfile` with `smooth_sizes[Nfile]`.
- `smooth_col`: smooth every panel in column `col` with `smooth_sizes[col]`.
- `annotate_smoothing`, `smooth_contour_col`: draw a white ellipse of diameter
  `smooth_sizes[Nfile]` in the lower left of each panel in column `col`. `smooth_col[Nfile]` also
  draws it.

Contours and streamlines:
- `contours`: overlay dashed contours, indexed by `selected`. The contour map of panel `Nfile` is
  read from `contour_files`, or taken from `contour_arr` and `contour_par_arr`.
- `contour_levels`: contour levels; values below the first level, NaN and Inf are set to it.
  `nothing` picks levels automatically.
- `contour_color`, `alpha_contours`: contour colour and opacity (`alpha_contours[selected]`).
- `smooth_contour_file`: smooth the contour map of panel `Nfile` with `smooth_sizes[Nfile]`.
- `streamlines`: overlay streamlines, indexed by `selected`, with x and y components read from
  the SPHtoGrid FITS files `streamline_files[Nfile]` and `streamline_files[Nfile+1]`.

Annotations (indexed by `Nfile` unless noted):
- `annotate_scale`, `scale_kpc`, `scale_label`: scale bar of length `scale_kpc` in map length
  units in the lower right, labelled `scale_label`.
- `annotate_time`, `time_labels`: label in the upper left.
- `annotate_text`, `text_labels`: label in the upper right.
- `annotation_color = "w"`: colour of scale bars and labels.
- `r_circles`, `circle_alpha`, `circle_lines`: circles of radius `r_circles` in map length units
  around the map centre on every panel, with linestyles `circle_lines` from the largest radius
  down.
- `overplotting_functions`: callbacks `f(ax)` in panel order, called after the image, contours
  and streamlines are drawn.

Layout:
- `colorbar_mode = "edge"`: `"edge"` or `"single"`, see above.
- `colorbar_location = "top"`: `"top"`, `"bottom"` or `"right"`.
- `colorbar_bottom = false`: with `colorbar_mode="single"`, put the colorbar at the bottom unless
  `colorbar_location="top"`.
- `colorbar_size = 25`: colorbar thickness in pixels. Log colorbars are ticked at whole decades.
- `grid_direction = "column"`: `"column"` or `"row"`, see above.
- `upscale = Ncols`: sets the panel size relative to the font size; the figure starts at
  `600upscale × 650upscale` per panel and is then shrunk to the layout.
- `rasterize = true`: rasterize the images in vector output.

`dpi`, `ticks_color`, `circle_color` and `shift_colorbar_labels_inward` are accepted for
compatibility with `PyPlotUtility.plot_image_grid` but unused.

# Example
B in the top row and divB in the bottom row for three runs, one colorbar per row on the right:
```julia
files = [folder * "/maps/snap_045.\$q.yz.fits" for folder ∈ folders, q ∈ ["B", "divB"]]

plot_image_grid(2, 3, files, ["magma", "inferno"], [L"B \\: [\\mu \\mathrm{G}]", L"h |\\nabla \\cdot B| / |B|"],
                [1.e-3, 1.e-2], [10.0, 20.0], "BdivB.pdf",
                cutoffs = [1.e-3, 1.e-2], grid_direction = "row", colorbar_location = "right")
```
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
                                colorbar_size = 25,
                                rasterize = true
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

            # right-hand edge colorbars are one per row, so their properties are selected by row
            per_row = colorbar_mode == "edge" && colorbar_location == "right"
            selected = colorbar_mode == "single" ? 1 : (per_row ? row : col)

            if !isnothing(cutoffs)
                icut = per_row ? row : col
                map[map .< cutoffs[icut]] .= cutoffs[icut]
                map[isnan.(map)] .= cutoffs[icut]
                map[isinf.(map)] .= cutoffs[icut]
            end

            # Image Attributes
            im_kwargs = Dict{Symbol, Any}()
            im_kwargs[:colormap] = Symbol(im_cmap[selected])
            im_kwargs[:colorrange] = (vmin_arr[selected], vmax_arr[selected])
            im_kwargs[:rasterize] = rasterize

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