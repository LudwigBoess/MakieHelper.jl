# Orthographic CPU volume renderer for 3D scalar fields: emission-absorption
# ray-march with front-to-back compositing over a black background. The result is
# an RGB image that CairoMakie plots with `image!`, standing in for the volume
# rendering CairoMakie lacks. Two integrators:
# - default: per-sample opacity as in entity's in-situ renderer
#   (src/output/render/raymarch.hpp), so the range and opacity settings of a scene
#   carry over to the in-situ [render] deck;
# - `l_ref` given: the biased ray march of Blender Cycles
#   (intern/cycles/kernel/integrator/shade_volume.h, "Biased" volume option) with
#   extinction and emission per unit length, independent of the step size.

_rt_vec3(v::Number) = (Float64(v), Float64(v), Float64(v))
_rt_vec3(v)         = (Float64(v[1]), Float64(v[2]), Float64(v[3]))
_rt_dot(a, b)       = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
_rt_cross(a, b)     = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
_rt_normalize(a)    = a ./ sqrt(_rt_dot(a, a))

# eye position and orthonormal basis (fwd, right, up) of cam = (eye, look_at, up, ortho_height)
function _rt_camera_basis(cam)
    eye   = _rt_vec3(cam[1])
    fwd   = _rt_normalize(_rt_vec3(cam[2]) .- eye)
    right = _rt_normalize(_rt_cross(fwd, _rt_vec3(cam[3])))
    upc   = _rt_cross(right, fwd)
    return eye, fwd, right, upc
end

# half width and half height of the image plane in world units
function _rt_half_extent(cam, W, H)
    half_h = Float64(cam[4]) / 2
    return half_h * W / H, half_h
end

# opacity at normalised value u: linear between (u, alpha) control points, alpha = u without points
function _rt_alpha_at(pts, u)
    isempty(pts) && return u
    u <= pts[1][1] && return pts[1][2]
    u >= pts[end][1] && return pts[end][2]
    for i = 1:length(pts)-1
        if pts[i][1] <= u <= pts[i+1][1]
            span = pts[i+1][1] - pts[i][1]
            t = span > 0 ? (u - pts[i][1]) / span : 0.0
            return pts[i][2] * (1 - t) + pts[i+1][2] * t
        end
    end
    return pts[end][2]
end

# premultiplied (r a, g a, b a, a) lookup table; colours from the Makie colormap, its own alpha ignored
function _rt_lut(cmap, alpha_pts, n_lut)
    cm  = Makie.to_colormap(cmap)
    lut = Vector{NTuple{4,Float32}}(undef, n_lut)
    for i = 1:n_lut
        u = (i - 1) / (n_lut - 1)
        c = Makie.interpolated_getindex(cm, u)
        a = Float32(_rt_alpha_at(alpha_pts, u))
        lut[i] = (red(c) * a, green(c) * a, blue(c) * a, a)
    end
    return lut
end

# (r, g, b, σ) lookup table for the Cycles integrator: straight colormap colour and the
# extinction per world unit of a slab of thickness l_ref whose opacity is alpha(u)
function _rt_extinction_lut(cmap, alpha_pts, n_lut, l_ref)
    cm  = Makie.to_colormap(cmap)
    lut = Vector{NTuple{4,Float64}}(undef, n_lut)
    for i = 1:n_lut
        u = (i - 1) / (n_lut - 1)
        c = Makie.interpolated_getindex(cm, u)
        a = Float64(_rt_alpha_at(alpha_pts, u))
        lut[i] = (Float64(red(c)), Float64(green(c)), Float64(blue(c)), -log1p(-a) / l_ref)
    end
    return lut
end

# trilinear sample of the cell-centred field F at world point p; positions beyond
# the outermost cell centres take the edge value
@inline function _rt_sample(F, x0, inv_dx, p)
    n1, n2, n3 = size(F)
    g1 = clamp((p[1] - x0[1]) * inv_dx[1] - 0.5, 0.0, n1 - 1.0)
    g2 = clamp((p[2] - x0[2]) * inv_dx[2] - 0.5, 0.0, n2 - 1.0)
    g3 = clamp((p[3] - x0[3]) * inv_dx[3] - 0.5, 0.0, n3 - 1.0)
    # g >= 0 after the clamp, so truncation is floor
    b1 = min(unsafe_trunc(Int, g1), n1 - 2)
    b2 = min(unsafe_trunc(Int, g2), n2 - 2)
    b3 = min(unsafe_trunc(Int, g3), n3 - 2)
    t1 = g1 - b1
    t2 = g2 - b2
    t3 = g3 - b3
    # linear index of the lower corner and the strides of dimensions 2 and 3
    s2 = n1
    s3 = n1 * n2
    i  = 1 + b1 + s2 * b2 + s3 * b3
    @inbounds begin
        c00 = F[i] * (1 - t1) + F[i+1] * t1
        c10 = F[i+s2] * (1 - t1) + F[i+s2+1] * t1
        c01 = F[i+s3] * (1 - t1) + F[i+s3+1] * t1
        c11 = F[i+s2+s3] * (1 - t1) + F[i+s2+s3+1] * t1
    end
    c0 = c00 * (1 - t2) + c10 * t2
    c1 = c01 * (1 - t2) + c11 * t2
    return c0 * (1 - t3) + c1 * t3
end

# slab intersection (t_enter, t_exit) of the ray o + t d with the box [lo, hi]
function _rt_ray_box(o, d, lo, hi)
    tmin, tmax = -Inf, Inf
    for k = 1:3
        if abs(d[k]) < 1e-12
            (o[k] < lo[k] || o[k] > hi[k]) && return (Inf, -Inf)
        else
            t1 = (lo[k] - o[k]) / d[k]
            t2 = (hi[k] - o[k]) / d[k]
            tmin = max(tmin, min(t1, t2))
            tmax = min(tmax, max(t1, t2))
        end
    end
    return (tmin, tmax)
end

# front-to-back compositing along one ray at the global sample positions t_k = k ds
function _rt_march(F, x0, inv_dx, o, d, te, tx, ds, lut, logscale, vlo, inv_range, early)
    nl = length(lut)
    t  = ceil(te / ds) * ds
    ar = ag = ab = aa = 0.0
    while t < tx
        p = (o[1] + t * d[1], o[2] + t * d[2], o[3] + t * d[3])
        s = _rt_sample(F, x0, inv_dx, p)
        u = logscale ? (s > 0 ? (log10(s) - vlo) * inv_range : -1.0) : (s - vlo) * inv_range
        u = isnan(u) ? 0.0 : clamp(u, 0.0, 1.0)
        @inbounds cr, cg, cb, ca = lut[unsafe_trunc(Int, u * (nl - 1) + 0.5) + 1]
        if ca > 0
            w   = 1 - aa
            ar += w * cr
            ag += w * cg
            ab += w * cb
            aa += w * ca
            aa >= early && break
        end
        t += ds
    end
    return ar, ag, ab
end

# Cycles' biased ray march along one ray: segments end at te + (i + offset) ds, the
# first one shortened by offset ∈ [0, 1]; each is shaded once at the fraction shade of
# its length, with transmittance exp(-σ dt) and emission σ c integrated exactly over dt
function _rt_march_segments(F, x0, inv_dx, o, d, te, tx, ds, tab, logscale, vlo, inv_range, offset, shade, early)
    nl = length(tab)
    L1 = L2 = L3 = 0.0
    Tr = 1.0
    lo = te
    i  = 0
    while lo < tx
        hi = min(tx, te + (i + offset) * ds)
        i += 1
        hi > lo || continue
        ts = lo + shade * (hi - lo)
        s  = _rt_sample(F, x0, inv_dx, (o[1] + ts * d[1], o[2] + ts * d[2], o[3] + ts * d[3]))
        u  = logscale ? (s > 0 ? (log10(s) - vlo) * inv_range : -1.0) : (s - vlo) * inv_range
        u  = isnan(u) ? 0.0 : clamp(u, 0.0, 1.0)
        @inbounds cr, cg, cb, σ = tab[unsafe_trunc(Int, u * (nl - 1) + 0.5) + 1]
        if σ > 0
            # segment opacity 1 - exp(-σ dt) = ∫ σ exp(-σ t) dt over the segment
            a   = -expm1(-σ * (hi - lo))
            w   = Tr * a
            L1 += w * cr
            L2 += w * cg
            L3 += w * cb
            Tr *= 1 - a
            1 - Tr >= early && break
        end
        lo = hi
    end
    return L1, L2, L3
end

# splitmix64 finaliser
@inline function _rt_mix(x::UInt64)
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end

# uniform number in [0, 1) for (seed, pixel, ray sample, dimension); a pure function of
# its arguments, so the image does not depend on thread scheduling
@inline function _rt_rand(seed, px, py, s, dim)
    h = _rt_mix((seed % UInt64) + 0x9e3779b97f4a7c15)
    h = _rt_mix(h ⊻ (px % UInt64))
    h = _rt_mix(h ⊻ (py % UInt64))
    h = _rt_mix(h ⊻ (((s % UInt64) << 8) | (dim % UInt64)))
    return (h >> 11) * 0x1.0p-53
end

# casts the ray of every pixel through the box [lo, hi] and stores f(o, te, tx, px, py) -> (r, g, b)
function _rt_pixels!(f, img, eye, fwd, right, upc, half_w, half_h, lo, hi)
    H, W = size(img, 1), size(img, 2)
    Threads.@threads for py = 1:H
        sy = (1 - 2 * (py - 0.5) / H) * half_h
        for px = 1:W
            sx = (2 * (px - 0.5) / W - 1) * half_w
            o  = eye .+ sx .* right .+ sy .* upc
            te, tx = _rt_ray_box(o, fwd, lo, hi)
            te >= tx && continue
            r, g, b = f(o, te, tx, px, py)
            img[py, px, 1] = r
            img[py, px, 2] = g
            img[py, px, 3] = b
        end
    end
    return img
end


"""
    render_volume(F, x0, dx, lo, hi, cam, W, H, scene; ds = 2.0, early = 0.99, n_lut = 256,
                  l_ref = nothing, spp = 1, jitter = true, seed = 0, max_steps = typemax(Int))

Volume-renders the cell-centred 3D scalar field `F` with an orthographic camera
and returns an `H × W × 3` `Float32` RGB image in [0, 1] (row 1 at the top),
composited over black. Plot it with [`rt_image!`](@ref) and [`rt_colorbar`](@ref).

- `x0`: world position of the lower corner of cell `F[1, 1, 1]`; `dx`: cell size,
  scalar or 3-tuple. Samples beyond the outermost cell centres take the edge value.
- `lo`, `hi`: corners of the rendered box in world coordinates.
- `cam = (eye, look_at, up, ortho_height)`: `ortho_height` is the world height
  covered by the image, the width follows from `W / H`.
- `scene = (vmin, vmax, log, cmap, alpha)`: colour range, log10 normalisation,
  any Makie colormap (`:inferno`, `Reverse(:viridis)`, `cgrad(...)`, ...) and the
  opacity transfer function as `(u, alpha)` control points over the normalised
  value `u ∈ [0, 1]`, linear in between (an empty list gives alpha = u).
- `ds`: step along the ray in world units.
- `early`: accumulated opacity at which a ray stops.
- `n_lut`: entries of the colour/opacity lookup table.

Integrators:

- `l_ref = nothing` (default), entity's in-situ renderer: samples at `t_k = k ds`
  from the eye, each with opacity `alpha(u)`, so the opacity of a structure grows
  as `ds` shrinks.
- `l_ref` set, the biased ray march of Blender Cycles: `alpha(u)` is the opacity of
  a slab of thickness `l_ref`, i.e. extinction `σ(u) = -ln(1 - alpha(u)) / l_ref`
  per world unit with emission `σ(u)` times the colormap colour. Segments of length
  `ds` start at the ray's entry into the box, the first one shortened by a random
  offset; each is shaded once at a random fraction of its length, with transmittance
  `exp(-σ dt)` and the emission integrated exactly over `dt`. The image converges
  as `ds` shrinks; at `ds = l_ref` its opacities match the default integrator.
  - `spp`: ray samples per pixel, averaged, with stratified segment offsets.
  - `jitter`: random offset and shading fraction per ray sample, which turns
    banding at sharp features into noise; `false` uses full segments from the
    entry point shaded at their midpoints (then `spp` must be 1).
  - `seed`: seed of the per-pixel random numbers (images do not depend on threading).
  - `max_steps`: segment count per ray above which `ds` is raised to
    `path length / max_steps` (Cycles' default: 1024).

Threaded over image rows.
"""
function render_volume(F::AbstractArray{<:Real,3}, x0, dx, lo, hi, cam, W::Integer, H::Integer, scene;
                       ds = 2.0, early = 0.99, n_lut = 256,
                       l_ref = nothing, spp::Integer = 1, jitter::Bool = true, seed::Integer = 0,
                       max_steps::Integer = typemax(Int))
    Base.require_one_based_indexing(F)
    minimum(size(F)) >= 2 || throw(ArgumentError("F needs at least 2 cells per dimension, got size $(size(F))"))
    ds > 0 || throw(ArgumentError("ds must be positive, got $ds"))
    if isnothing(l_ref)
        spp == 1 || throw(ArgumentError("spp > 1 needs the Cycles integrator (l_ref)"))
    else
        l_ref > 0 || throw(ArgumentError("l_ref must be positive, got $l_ref"))
        spp >= 1 || throw(ArgumentError("spp must be at least 1, got $spp"))
        (jitter || spp == 1) || throw(ArgumentError("spp > 1 needs jitter = true"))
        max_steps >= 1 || throw(ArgumentError("max_steps must be at least 1, got $max_steps"))
    end

    eye, fwd, right, upc = _rt_camera_basis(cam)
    half_w, half_h = _rt_half_extent(cam, W, H)
    x0v    = _rt_vec3(x0)
    inv_dx = 1 ./ _rt_vec3(dx)
    lov    = _rt_vec3(lo)
    hiv    = _rt_vec3(hi)

    logscale  = Bool(scene.log)
    vlo       = logscale ? log10(Float64(scene.vmin)) : Float64(scene.vmin)
    span      = (logscale ? log10(Float64(scene.vmax)) : Float64(scene.vmax)) - vlo
    inv_range = span > 0 ? 1 / span : 0.0
    dsf       = Float64(ds)
    earlyf    = Float64(early)

    img = zeros(Float32, H, W, 3)
    if isnothing(l_ref)
        lut = _rt_lut(scene.cmap, scene.alpha, n_lut)
        _rt_pixels!(img, eye, fwd, right, upc, half_w, half_h, lov, hiv) do o, te, tx, px, py
            _rt_march(F, x0v, inv_dx, o, fwd, te, tx, dsf, lut, logscale, vlo, inv_range, earlyf)
        end
    else
        tab = _rt_extinction_lut(scene.cmap, scene.alpha, n_lut, Float64(l_ref))
        _rt_pixels!(img, eye, fwd, right, upc, half_w, half_h, lov, hiv) do o, te, tx, px, py
            ds_ray = tx - te > max_steps * dsf ? (tx - te) / max_steps : dsf
            r = g = b = 0.0
            for s = 1:spp
                offset = jitter ? (s - 1 + _rt_rand(seed, px, py, s, 1)) / spp : 1.0
                shade  = jitter ? _rt_rand(seed, px, py, s, 2) : 0.5
                cr, cg, cb = _rt_march_segments(F, x0v, inv_dx, o, fwd, te, tx, ds_ray, tab, logscale, vlo,
                                                inv_range, offset, shade, earlyf)
                r += cr
                g += cg
                b += cb
            end
            return r / spp, g / spp, b / spp
        end
    end
    return img
end


"""
    rt_image(img)

Converts the `H × W × 3` output of [`render_volume`](@ref) (row 1 at the top) into
a `W × H` `Matrix{RGBf}` in Makie's `image` orientation (first index along x, y up).
"""
rt_image(img) = rotr90(RGBf.(view(img, :, :, 1), view(img, :, :, 2), view(img, :, :, 3)))


"""
    rt_image!(ax, img, cam; interpolate = false, kwargs...)

Plots a [`render_volume`](@ref) image into `ax`, spanning the image plane in
world units centred on the view axis: `x ∈ ±ortho_height W / 2H`,
`y ∈ ±ortho_height / 2`. World points map into these coordinates with
[`rt_project`](@ref); the orthographic scale is uniform, so a scale bar of
length `L` is `L` axis units long anywhere in the image. `kwargs` go to `image!`.
"""
function rt_image!(ax, img, cam; interpolate = false, kwargs...)
    half_w, half_h = _rt_half_extent(cam, size(img, 2), size(img, 1))
    return image!(ax, -half_w .. half_w, -half_h .. half_h, rt_image(img); interpolate, kwargs...)
end


"""
    rt_project(cam, p)

Position of the world point `p` in the axis coordinates of [`rt_image!`](@ref),
as a `Point2f`.
"""
function rt_project(cam, p)
    eye, _, right, upc = _rt_camera_basis(cam)
    d = _rt_vec3(p) .- eye
    return Point2f(_rt_dot(d, right), _rt_dot(d, upc))
end


"""
    rt_colorbar(pos, scene; kwargs...)

`Colorbar` at layout position `pos` with the colormap, range and log scaling of a
[`render_volume`](@ref) scene. It maps hue to field value; pixel brightness in
the rendering also carries the opacity accumulated along the ray. `kwargs` go to
`Colorbar`.
"""
rt_colorbar(pos, scene; kwargs...) = Colorbar(pos; colormap = scene.cmap, limits = (scene.vmin, scene.vmax),
                                              scale = scene.log ? log10 : identity, kwargs...)
