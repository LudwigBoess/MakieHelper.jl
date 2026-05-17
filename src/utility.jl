
function read_map_par(read_mode, Nfile, files, map_arr, par_arr)
    if isnothing(map_arr)
        image_name = files[Nfile]

        # read SPHtoGrid image
        if read_mode == 1
            map, par, snap_num, units = read_fits_image(image_name)

        # read Smac1 binary image
        elseif read_mode == 2
            map = read_smac1_binary_image(image_name)
            smac1_info = read_smac1_binary_info(image_name)
            smac1_center = [smac1_info.xcm, smac1_info.ycm, smac1_info.zcm] ./ 3.085678e21
            par = mappingParameters(center = smac1_center,
                x_size = smac1_info.boxsize_kpc,
                y_size = smac1_info.boxsize_kpc,
                z_size = smac1_info.boxsize_kpc,
                Npixels = smac1_info.boxsize_pix)

        elseif read_mode == 3
            map = read_smac2_image(image_name, 1)
            par = read_smac2_info(image_name)
            
        elseif read_mode == 4
            map = read_smac1_fits_image(image_name, 1)
            par = read_smac1_fits_info(image_name)
        end
    else
        map = map_arr[Nfile]
        par = par_arr[Nfile]
    end

    # transpose map to match expected Makie / image mapping behavior (rows/cols)
    return map, par
end

function smooth_map!(map, smooth_size, par)
    pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]
    smooth_pixel = smooth_size ./ pixelSideLength ./ 2
    map = imfilter(map, reflect(Kernel.gaussian((smooth_pixel[1], smooth_pixel[2]), (par.Npixels[1] - 1, par.Npixels[1] - 1))))
    return map
end

# Replaces `get_circle_contours` with standard geometric plotting
function add_circle_contours!(ax, r_circles, circle_labels, circle_alpha, circle_lines, par)
    pixelSideLength = (par.x_lim[2] - par.x_lim[1]) / par.Npixels[1]
    center = Point2f(0.5par.Npixels[1], 0.5par.Npixels[2])

    sorted = reverse(sortperm(r_circles))
    
    for (idx, i) in enumerate(sorted)
        # circle radius in pixel units
        r = Float32(r_circles[i] / pixelSideLength)
        circ = Circle(center, r)
        
        # Draw the circle natively
        lines!(ax, circ, color=(:white, circle_alpha), linestyle=Symbol(circle_lines[idx]), linewidth=1.5)
        
        # Add label text inline
        text!(ax, string(circle_labels[i]), position = center .+ Point2f(0, r), 
              color=(:white, circle_alpha), align=(:center, :bottom), fontsize=14)
    end
    return ax
end