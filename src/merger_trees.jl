"""
Merger-tree visualisation.

Ports `galaxybox/visualization/tree.py` (`TemporalTreePlotter`) to Julia
and renders with Makie. Operates on a `GadgetIO.MergerTree`, which holds
the depth-first halo array produced by L-HaloTrees (or the Julia
re-implementation in `GadgetIO`).

The data plumbing mirrors galaxybox:

| galaxybox column | source in `MergerTree`                                       |
| :--------------- | :----------------------------------------------------------- |
| `Scale`          | user-supplied `snap → scale` map (default: linear)           |
| `Desc_ID`        | `Descendant` (0-based pointer; converted to 1-based)         |
| `MMP`            | `is_mmp(tree, i)`                                            |
| `Num_prog`       | `num_progenitors(tree, i)`                                   |
| `Coprog_ID`      | `NextProgenitor` of the current node                         |
| `Stellar_mass`   | user-supplied `mass_of(halo)` (default: `log10(Len)`)        |
"""

"""
    struct MergerTreeLayout

Result of running the galaxybox-style layout on a `MergerTree`.

* `x[i]`    — horizontal coordinate (scale factor) of halo `i`
* `y[i]`    — vertical coordinate set by the recursive walk; `NaN` when
              the halo did not pass the `min_scale`/`max_mu` cuts
* `mass[i]` — colour value, defaults to `log10(Len)`
* `edges`   — vector of `(child_idx, parent_idx)` pairs (1-based) to draw
"""
struct MergerTreeLayout
    x::Vector{Float64}
    y::Vector{Float64}
    mass::Vector{Float64}
    edges::Vector{Tuple{Int,Int}}
end


"""
    tree_layout(tree::MergerTree;
                snap_to_scale = nothing,
                mass_of = h -> log10(max(Int(h.Len), 1)),
                min_scale = 0.0, max_mu = 10_000.0,
                plot_style = :centered) -> MergerTreeLayout

Compute `(x, y, mass, edges)` for plotting a single merger tree, using
the algorithm from `galaxybox.visualization.tree.TemporalTreePlotter`:

1. From the root, recurse along the main (MMP) branch first, setting
   `y` to the current level.
2. At each level, walk co-progenitors whose mass ratio
   `mu = 10^(M_MMP - M_coprog)` does not exceed `max_mu`, advancing the
   y cursor.
3. When `plot_style = :centered`, on the way back up, descendants of
   non-MMPs sit at the midpoint between the current y and the MMP's y;
   along an MMP chain the descendant simply inherits its progenitor's y.

`snap_to_scale` converts a halo's `SnapNum` to a scale factor. If
`nothing`, defaults to `(snap+1)/(last_snap+1)` where `last_snap` is the
largest SnapNum in the tree.

Halos with `scale ≤ min_scale` are gated out and end up with `y = NaN`;
plotting filters them.
"""
function tree_layout(tree::GadgetIO.MergerTree;
                     snap_to_scale = nothing,
                     mass_of = h -> log10(max(Int(h.Len), 1)),
                     min_scale::Real = 0.0,
                     max_mu::Real = 10_000.0,
                     plot_style::Symbol = :centered)

    plot_style ∈ (:centered, :linear) ||
        throw(ArgumentError("plot_style must be :centered or :linear"))

    n = length(tree)

    if snap_to_scale === nothing
        last_snap = maximum(Int(h.SnapNum) for h in tree.halos)
        snap_to_scale = s -> GadgetIO.snap_to_scale_default(s, last_snap)
    end

    x = [Float64(snap_to_scale(Int(h.SnapNum))) for h in tree.halos]
    mass = [Float64(mass_of(h)) for h in tree.halos]
    y = fill(NaN, n)

    root = GadgetIO.root_index(tree)
    current_y = Ref(0)
    node = Ref(root)

    function set_pos!()
        i = node[]
        y[i] = current_y[]
        scale_now = x[i]

        # MMP of self.node (in galaxybox terms): our FirstProgenitor
        immp = 0
        if GadgetIO.num_progenitors(tree, i) > 0
            immp = Int(tree[i].FirstProgenitor) + 1
        end

        # Coprog_ID of self.node: our NextProgenitor (subject to mu cut)
        icoprog = _coprog_for_layout(tree, i, mass, max_mu)

        if scale_now > min_scale
            if immp > 0
                node[] = immp
                set_pos!()
            end
            if icoprog > 0
                current_y[] += 1
                node[] = icoprog
                set_pos!()
            end

            if plot_style === :centered
                _update_desc_pos!(tree, y, node)
            end
        end
    end

    set_pos!()

    edges = Tuple{Int,Int}[]
    for i in 1:n
        isnan(y[i]) && continue
        d = tree[i].Descendant
        d < 0 && continue
        push!(edges, (i, Int(d) + 1))
    end

    return MergerTreeLayout(x, y, mass, edges)
end


function _coprog_for_layout(tree::GadgetIO.MergerTree, i_julia::Integer,
                            mass::Vector{Float64}, max_mu::Real)
    # galaxybox's `Coprog_ID` at row i = next halo in i's progenitor
    # linked list, i.e. our `NextProgenitor`. The mu cut compares the
    # coprog against the MMP of i's descendant (i.e. the head of the
    # same chain).
    h = tree[i_julia]
    h.NextProgenitor < 0 && return 0
    h.Descendant < 0 && return Int(h.NextProgenitor) + 1

    # `set_progenitor_pointers!` only guarantees the MMP sits at the *head*
    # of each progenitor list; the remaining co-progenitors are NOT sorted
    # by mass. So we must not stop at the first co-progenitor that fails the
    # mass-ratio cut — a large secondary (a real merger) can sit further down
    # the chain behind tiny accreted halos, and stopping early would drop the
    # whole branch. Walk the sibling chain and return the first co-progenitor
    # that passes; only an exhausted chain yields 0. For already mass-sorted
    # trees (L-HaloTrees) this is equivalent to the old early stop.
    immp_of_desc = Int(tree[h.Descendant + 1].FirstProgenitor) + 1
    icop_julia = Int(h.NextProgenitor) + 1
    while true
        mu = 10.0 ^ (mass[immp_of_desc] - mass[icop_julia])
        mu <= max_mu && return icop_julia
        nxt = tree[icop_julia].NextProgenitor
        nxt < 0 && return 0
        icop_julia = Int(nxt) + 1
    end
end


function _update_desc_pos!(tree::GadgetIO.MergerTree, y::Vector{Float64}, node::Ref{Int})
    # Direct port of galaxybox `_update_desc_pos` with max_depth = 1:
    # at most one non-MMP centering, but unlimited MMP-following.
    max_depth = 1
    depth = 0
    while true
        desc_ptr = tree[node[]].Descendant
        desc_ptr < 0 && break
        depth < max_depth || break

        desc_julia = Int(desc_ptr) + 1
        if !GadgetIO.is_mmp(tree, node[])
            immp_of_desc = Int(tree[desc_julia].FirstProgenitor) + 1
            mmp_y = y[immp_of_desc]
            now_y = y[node[]]
            if !isnan(mmp_y) && !isnan(now_y)
                y[desc_julia] = (now_y + mmp_y) / 2
            end
            depth += 1
        else
            y[desc_julia] = y[node[]]
        end
        node[] = desc_julia
    end
end


"""
    plot_merger_tree(tree::MergerTree; ax=nothing, kwargs...) -> (fig, ax)

Draw a merger tree in galaxybox style. If `ax` is given, draw into that
existing Makie axis (no figure or colorbar is created). Otherwise a new
`Figure` with axis and colorbar is returned.

Keyword arguments forwarded to [`tree_layout`](@ref):
`snap_to_scale`, `mass_of`, `min_scale`, `max_mu`, `plot_style`.

Plot-only arguments:

| Name           | Default                                                  |
| :------------- | :------------------------------------------------------- |
| `colormap`     | `:jet`                                                   |
| `vmin`/`vmax`  | autoscale from valid masses                              |
| `markersize`   | `8`                                                      |
| `edge_color`   | `:black`                                                 |
| `cb_label`     | `L"\\log_{10}(M_{\\ast}/M_\\odot)"`                      |
| `figsize`      | `(1200, 600)`                                            |
| `xlabel`       | `\"Scale factor\"`                                       |
| `redshift_axis` | `true` — overlay a secondary top axis showing redshift  |
| `redshift_ticks` | `[10, 4, 2, 1, 0.5, 0.2, 0.1, 0.03, 0.0]`              |
| `redshift_label` | `\"Redshift z\"`                                       |

The redshift axis is only added when the function owns the figure
(`ax = nothing`). For user-supplied axes call [`add_redshift_axis!`](@ref).
When the redshift axis is active, set the title via `ax_z.title = …` on
the returned secondary axis — `ax.title` would overlap its tick labels.

# Example (L-HaloTrees binary)
```julia
using GadgetIO, MakieHelper, CairoMakie
f = read_merger_tree_file("trees_063.0")
plot_merger_tree(f[1]; mass_of = h -> log10(h.M_Crit200), max_mu = 1e3)
```

# Example (consistent-trees ASCII)
`GadgetIO.read_consistent_trees` returns a `Vector{MergerTree}` built with
exactly the same pointer convention as L-HaloTrees, so it plots with no
changes here. Use `GadgetIO.scale_factor_map` to put the x-axis on the
real per-snapshot scale factor instead of the `(snap+1)/(last_snap+1)`
heuristic:
```julia
using GadgetIO, MakieHelper, CairoMakie
trees     = read_consistent_trees("rockstar_run/trees/tree_0_0_0.dat")
scaledict = scale_factor_map("rockstar_run/trees/tree_0_0_0.dat")
plot_merger_tree(trees[1];
                 mass_of       = h -> log10(h.M_Crit200),
                 snap_to_scale = s -> scaledict[s],
                 max_mu        = 1e3)
```
"""
function plot_merger_tree(tree::GadgetIO.MergerTree;
                          ax = nothing,
                          snap_to_scale = nothing,
                          mass_of = h -> log10(max(Int(h.Len), 1)),
                          min_scale::Real = 0.0,
                          max_mu::Real = 10_000.0,
                          plot_style::Symbol = :centered,
                          colormap = :jet,
                          vmin = nothing, vmax = nothing,
                          markersize::Real = 8,
                          edge_color = :black,
                          cb_label = L"\log_{10}(M_{\ast}/M_\odot)",
                          figsize::Tuple{Int,Int} = (1200, 600),
                          xlabel::AbstractString = "Scale factor",
                          redshift_axis::Bool = true,
                          redshift_ticks::AbstractVector{<:Real} =
                              [10, 4, 2, 1, 0.5, 0.2, 0.1, 0.03, 0.0],
                          redshift_label::AbstractString = "Redshift z")

    layout = tree_layout(tree;
                         snap_to_scale, mass_of,
                         min_scale, max_mu, plot_style)
    return plot_merger_tree(layout;
                            ax, colormap, vmin, vmax, markersize,
                            edge_color, cb_label, figsize, xlabel,
                            redshift_axis, redshift_ticks, redshift_label)
end


"""
    plot_merger_tree(layout::MergerTreeLayout; ax=nothing, kwargs...) -> (fig, ax)

Variant that takes an already-computed `MergerTreeLayout`. Useful when
the same layout is reused across multiple renderings.
"""
function plot_merger_tree(layout::MergerTreeLayout;
                          ax = nothing,
                          colormap = :jet,
                          vmin = nothing, vmax = nothing,
                          markersize::Real = 8,
                          edge_color = :black,
                          cb_label = L"\log_{10}(M_{\ast}/M_\odot)",
                          figsize::Tuple{Int,Int} = (1200, 600),
                          xlabel::AbstractString = "Scale factor",
                          redshift_axis::Bool = true,
                          redshift_ticks::AbstractVector{<:Real} =
                              [10, 4, 2, 1, 0.5, 0.2, 0.1, 0.03, 0.0],
                          redshift_label::AbstractString = "Redshift z")

    valid = .!isnan.(layout.y)
    mvals = layout.mass[valid]

    vmin_eff = vmin === nothing ? minimum(mvals; init = 0.0) : Float64(vmin)
    vmax_eff = vmax === nothing ? maximum(mvals; init = 1.0) : Float64(vmax)

    own_fig = ax === nothing
    if own_fig
        fig = Figure(size = figsize)
        ax = Axis(fig[1, 1]; xlabel = xlabel)
        hideydecorations!(ax; grid = false)
    else
        fig = nothing
    end

    # edges first so the scatter sits on top
    edge_x = Float64[]
    edge_y = Float64[]
    for (c, p) in layout.edges
        push!(edge_x, layout.x[c], layout.x[p], NaN)
        push!(edge_y, layout.y[c], layout.y[p], NaN)
    end
    lines!(ax, edge_x, edge_y; color = edge_color)

    # nodes coloured by mass
    sc = scatter!(ax,
                  layout.x[valid], layout.y[valid];
                  color = mvals,
                  colormap = colormap,
                  colorrange = (vmin_eff, vmax_eff),
                  markersize = markersize,
                  strokewidth = 0)

    if own_fig
        if redshift_axis
            add_redshift_axis!(ax;
                               redshift_ticks = redshift_ticks,
                               label = redshift_label)
        end
        Colorbar(fig[1, 2], sc; label = cb_label)
        return fig, ax
    else
        return nothing, ax
    end
end


"""
    add_redshift_axis!(ax::Axis;
                       redshift_ticks = [10, 4, 2, 1, 0.5, 0.2, 0.1, 0.03, 0.0],
                       label = "Redshift z") -> Axis

Overlay a secondary x-axis at the top of `ax` showing redshift `z`,
assuming the bottom axis is scale factor `a = 1/(1+z)`. Each `z` in
`redshift_ticks` is placed at `a = 1/(1+z)`. The new axis is linked to
`ax`'s x-range so it tracks zoom/pan. Returns the new top axis — set its
`.title` if a title is needed (using `ax.title` would clash with the
top-axis tick labels).

`ax` must be placed in a `GridLayout` (the standard `fig[i, j]` case);
that placement is what we read to overlay the twin axis on the same cell.
"""
function add_redshift_axis!(ax::Axis;
                            redshift_ticks::AbstractVector{<:Real} =
                                [10, 4, 2, 1, 0.5, 0.2, 0.1, 0.03, 0.0],
                            label::AbstractString = "Redshift z")

    gc = ax.layoutobservables.gridcontent[]
    gc === nothing && throw(ArgumentError(
        "add_redshift_axis! requires ax to be placed in a GridLayout"))

    a_positions = Float64[1 / (1 + z) for z in redshift_ticks]
    z_labels    = String[_z_tick_label(z) for z in redshift_ticks]

    ax_z = Axis(gc.parent[gc.span.rows, gc.span.cols];
                xaxisposition   = :top,
                xlabel          = label,
                xticks          = (a_positions, z_labels),
                xgridvisible    = false,
                ygridvisible    = false,
                xticksmirrored  = false)
    hidespines!(ax_z)
    hideydecorations!(ax_z)
    # the parent theme may mirror ticks; that draws the a-ticks across
    # the top axis and the z-ticks across the bottom — kill both sides
    ax.xticksmirrored = false
    linkxaxes!(ax, ax_z)
    return ax_z
end


_z_tick_label(z::Real) =
    z == round(z) ? string(Int(round(z))) : string(z)


# ===========================================================================
# Forest plot
# ===========================================================================

_forest_trees(f::GadgetIO.MergerForest)                       = f.trees
_forest_trees(f::AbstractVector{<:GadgetIO.MergerTree})       = collect(f)

"""
    plot_merger_forest(forest; ax=nothing, kwargs...) -> (fig, ax)

Render every tree of a consistent-trees forest in a single panel, each tree
laid out with [`tree_layout`](@ref) and stacked vertically (separated by
`gap`). `forest` may be a `GadgetIO.MergerForest` or any vector of
`GadgetIO.MergerTree`.

This is the view that exposes an **ongoing, not-yet-coalesced major merger**:
the main halo and the infalling companion are distinct z=0 roots in the same
forest, so they show up as two massive branches converging toward the present
day rather than as a single merger node (which would only appear once the
companion is fully disrupted). Colour encodes `mass_of` on a colour scale
shared across all trees.

Keyword arguments forwarded to [`tree_layout`](@ref) (applied per tree):
`snap_to_scale`, `mass_of`, `min_scale`, `max_mu`, `plot_style`. If
`snap_to_scale === nothing`, a single default mapping is built from the
largest `SnapNum` across the whole forest so all trees share one x-axis.

Forest-specific arguments:

| Name          | Default | Meaning                                              |
| :------------ | :------ | :--------------------------------------------------- |
| `gap`         | `2.0`   | vertical spacing (layout units) between stacked trees |
| `max_trees`   | `nothing` | plot only the `max_trees` largest trees (by halos)  |
| `label_trees` | `true`  | annotate each tree's z=0 root with `mass_of(root)`   |
| `label_fmt`   | `m -> @sprintf("%.1f", m)` | formatter for the root labels        |

All remaining plot-only arguments (`colormap`, `vmin`, `vmax`, `markersize`,
`edge_color`, `cb_label`, `figsize`, `xlabel`, `redshift_axis`,
`redshift_ticks`, `redshift_label`) match [`plot_merger_tree`](@ref).

# Example
```julia
using GadgetIO, MakieHelper, CairoMakie
treedir = "A2256/cr_mhd2/rockstar/trees"
forest  = build_merger_forest(treedir; ntrees_max = 6)
amap    = scale_factor_map(joinpath(treedir, "tree_0_0_0.dat"))
fig, ax = plot_merger_forest(forest;
                             snap_to_scale = s -> amap[s],
                             mass_of       = h -> log10(h.M_TopHat),
                             max_mu        = 1e4)
```
"""
function plot_merger_forest(forest;
                            snap_to_scale = nothing,
                            mass_of = h -> log10(max(Int(h.Len), 1)),
                            min_scale::Real = 0.0,
                            max_mu::Real = 10_000.0,
                            plot_style::Symbol = :centered,
                            gap::Real = 2.0,
                            max_trees::Union{Integer,Nothing} = nothing,
                            label_trees::Bool = true,
                            label_fmt = m -> @sprintf("%.1f", m),
                            colormap = :jet,
                            vmin = nothing, vmax = nothing,
                            markersize::Real = 8,
                            edge_color = :black,
                            cb_label = L"\log_{10}(M/M_\odot)",
                            figsize::Tuple{Int,Int} = (1600, 800),
                            xlabel::AbstractString = "Scale factor",
                            redshift_axis::Bool = true,
                            redshift_ticks::AbstractVector{<:Real} =
                                [10, 4, 2, 1, 0.5, 0.2, 0.1, 0.03, 0.0],
                            redshift_label::AbstractString = "Redshift z")

    trees = _forest_trees(forest)
    if max_trees !== nothing && length(trees) > max_trees
        ord   = sortperm(trees; by = length, rev = true)
        trees = trees[ord[1:max_trees]]
    end
    isempty(trees) && throw(ArgumentError("forest has no trees to plot"))

    # one shared snap → scale mapping so every tree lands on the same x-axis
    if snap_to_scale === nothing
        last_snap = maximum(maximum(Int(h.SnapNum) for h in t.halos) for t in trees)
        snap_to_scale = s -> GadgetIO.snap_to_scale_default(s, last_snap)
    end

    layouts = [tree_layout(t; snap_to_scale, mass_of,
                           min_scale, max_mu, plot_style) for t in trees]

    X = Float64[]; Y = Float64[]; M = Float64[]
    E = Tuple{Int,Int}[]
    lbl_x = Float64[]; lbl_y = Float64[]; lbl_s = String[]

    running = 0.0
    idxoff  = 0
    for (t, L) in zip(trees, layouts)
        valid = .!isnan.(L.y)
        ylo, yhi = any(valid) ? extrema(L.y[valid]) : (0.0, 0.0)
        shift = running - ylo

        append!(X, L.x)
        append!(Y, [isnan(y) ? NaN : y + shift for y in L.y])
        append!(M, L.mass)
        for (c, p) in L.edges
            push!(E, (c + idxoff, p + idxoff))
        end

        if label_trees
            r = GadgetIO.root_index(t)
            if !isnan(L.y[r])
                push!(lbl_x, L.x[r])
                push!(lbl_y, L.y[r] + shift)
                push!(lbl_s, label_fmt(L.mass[r]))
            end
        end

        running += (yhi - ylo) + gap
        idxoff  += length(L.y)
    end

    combined = MergerTreeLayout(X, Y, M, E)
    fig, ax = plot_merger_tree(combined;
                               colormap, vmin, vmax, markersize, edge_color,
                               cb_label, figsize, xlabel,
                               redshift_axis, redshift_ticks, redshift_label)

    if label_trees && !isempty(lbl_x)
        text!(ax, lbl_x, lbl_y; text = lbl_s,
              align = (:left, :center), offset = (6, 0),
              fontsize = 14, color = :black)
        # pad the x-range on the right so the z=0 labels are not clipped
        xv = [X[i] for i in eachindex(X) if !isnan(Y[i])]
        if !isempty(xv)
            xmin, xmax = extrema(xv)
            span = xmax - xmin
            span > 0 && xlims!(ax, xmin - 0.02span, xmax + 0.12span)
        end
    end

    return fig, ax
end
