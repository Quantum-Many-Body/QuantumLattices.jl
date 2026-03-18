module QuantumLatticesMakieExt

using Makie: Makie, Axis, Axis3, Colorbar, Figure, Label, lines!, scatter!, heatmap!, text!, linesegments!
using Makie: axislegend, DataAspect
using QuantumLattices: AbstractLattice, Assignment, Bond, BrillouinZone, Data, FractionalReciprocalSpace, Neighbors, Point, ReciprocalCurve, ReciprocalPath, ReciprocalScatter, ReciprocalZone
using QuantumLattices: atol, bonds, dimension, distance, fractionals, getcontent, label, rank, scalartype, str, ticks

# Helper function to extract the main axis from a Figure
function get_axis(fig::Figure)
    # Try to find an existing Axis or Axis3 in the figure
    for content in fig.content
        if isa(content, Axis) || isa(content, Axis3)
            return content
        end
    end
    error("No axis found in figure")
end

# Utility functions

function get_lattice_name(lattice::AbstractLattice)
    return String(getcontent(lattice, :name))
end

function coordinate_tuple(p, dim)
    if dim == 1
        return (p[1],)
    elseif dim == 2
        return (p[1], p[2])
    else
        return (p[1], p[2], p[3])
    end
end

function get_reciprocal_labels(path::ReciprocalPath)
    x, y, text = [], [], []
    for i = 1:length(path.contents)
        if length(path.contents[i])>0
            push!(x, path.contents[i][1][1])
            push!(y, path.contents[i][1][2])
            push!(text, string(path.labels[i].first))
            j = i%length(path.labels)+1
            if length(path.contents[i])>1 && path.labels[i].second!=path.labels[j].first
                push!(x, path.contents[i][end][1])
                push!(y, path.contents[i][end][2])
                push!(text, string(path.labels[i].second))
            end
        end
    end
    return x, y, text
end

# 1. Lattice plotting - plot! method for in-place plotting
function Makie.plot!(ax::Union{Axis, Axis3}, lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; siteon=false)
    dim = dimension(lattice)

    # Plot sites as scatter
    if isa(neighbors, Int) || 0 in keys(neighbors)
        site_coords = NTuple{dim, scalartype(lattice)}[]
        site_labels = String[]
        for i in eachindex(lattice)
            bond = Bond(Point(i, lattice[i], zero(lattice[i])))
            if filter(bond)
                push!(site_coords, coordinate_tuple(lattice[i], dim))
                push!(site_labels, siteon ? string(i) : "")
            end
        end

        if length(site_coords) > 0
            if dim == 1
                x_coords = [c[1] for c in site_coords]
                scatter!(ax, x_coords; marker=:circle, strokewidth=0)
            elseif dim == 2
                x_coords = [c[1] for c in site_coords]
                y_coords = [c[2] for c in site_coords]
                scatter!(ax, x_coords, y_coords; marker=:circle, strokewidth=0)
                # Add annotations if siteon
                if siteon
                    for (i, label) in enumerate(site_labels)
                        if !isempty(label)
                            text!(ax, x_coords[i], y_coords[i]; text=label, fontsize=8)
                        end
                    end
                end
            else
                x_coords = [c[1] for c in site_coords]
                y_coords = [c[2] for c in site_coords]
                z_coords = [c[3] for c in site_coords]
                scatter!(ax, x_coords, y_coords, z_coords; marker=:circle, strokewidth=0)
            end
        end
    end

    # Plot bonds as lines
    for bond in bonds(lattice, neighbors)
        if length(bond)==2 && filter(bond)
            if dim == 1
                x1, x2 = bond[1].rcoordinate[1], bond[2].rcoordinate[1]
                lines!(ax, [x1, x2], [0.0, 0.0])
            elseif dim == 2
                x1, x2 = bond[1].rcoordinate[1], bond[2].rcoordinate[1]
                y1, y2 = bond[1].rcoordinate[2], bond[2].rcoordinate[2]
                lines!(ax, [x1, x2], [y1, y2])
            else
                x1, x2 = bond[1].rcoordinate[1], bond[2].rcoordinate[1]
                y1, y2 = bond[1].rcoordinate[2], bond[2].rcoordinate[2]
                z1, z2 = bond[1].rcoordinate[3], bond[2].rcoordinate[3]
                lines!(ax, [x1, x2], [y1, y2], [z1, z2])
            end
        end
    end

    return ax
end

function Makie.plot(lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; siteon=false, figure_kw=(;), axis_kw=(;))
    dim = dimension(lattice)
    fig = Figure(; figure_kw...)

    if dim == 1
        ax = Axis(fig[1, 1]; title=get_lattice_name(lattice), titlesize=10, axis_kw...)
    elseif dim == 2
        ax = Axis(fig[1, 1]; title=get_lattice_name(lattice), titlesize=10, aspect=DataAspect(), axis_kw...)
    else
        ax = Axis3(fig[1, 1]; title=get_lattice_name(lattice), titlesize=10)
    end

    Makie.plot!(ax, lattice, neighbors, filter; siteon=siteon)
    return fig
end

# 2. FractionalReciprocalSpace plotting (BrillouinZone, ReciprocalZone, ReciprocalScatter)
function Makie.plot!(ax::Union{Axis, Axis3}, reciprocalspace::FractionalReciprocalSpace; fractional=false, autolims=true)
    N = rank(reciprocalspace)
    dim = dimension(reciprocalspace)

    if fractional && N <= 2
        coords_list = fractionals(reciprocalspace)
        coordinates = [tuple(c...) for c in coords_list]
    else
        coordinates = [tuple(c...) for c in reciprocalspace]
    end

    if dim == 1 || (fractional && N == 1)
        isa(ax, Axis) || error("Axis required for 1D reciprocal space")
        x_vals = [c[1] for c in coordinates]
        if !autolims
            Makie.xlims!(ax, extrema(x_vals)...)
        end
        scatter!(ax, x_vals; marker=:circle, strokewidth=0)
    elseif dim == 2 || (fractional && N == 2)
        isa(ax, Axis) || error("Axis required for 2D reciprocal space")
        x_vals = [c[1] for c in coordinates]
        y_vals = [c[2] for c in coordinates]
        if !autolims
            Makie.xlims!(ax, extrema(x_vals)...)
            Makie.ylims!(ax, extrema(y_vals)...)
        end
        scatter!(ax, x_vals, y_vals; marker=:circle, strokewidth=0)
        ax.xgridvisible = true
        ax.ygridvisible = true
        ax.xminorgridvisible = true
        ax.yminorgridvisible = true
    else
        isa(ax, Axis3) || error("Axis3 required for 3D reciprocal space")
        x_vals = [c[1] for c in coordinates]
        y_vals = [c[2] for c in coordinates]
        z_vals = [c[3] for c in coordinates]
        if !autolims
            Makie.limits!(ax, extrema(x_vals)..., extrema(y_vals)..., extrema(z_vals)...)
        end
        scatter!(ax, x_vals, y_vals, z_vals; marker=:circle, strokewidth=0)
    end

    return ax
end

function Makie.plot(reciprocalspace::FractionalReciprocalSpace; fractional=false, autolims=true, figure_kw=(;), axis_kw=(;))
    fig = Figure(; figure_kw...)

    N = rank(reciprocalspace)
    dim = dimension(reciprocalspace)

    if dim == 1 || (fractional && N == 1)
        ax = Axis(fig[1, 1]; titlesize=10, axis_kw...)
        ax.xlabel = label(reciprocalspace, 1)
    elseif dim == 2 || (fractional && N == 2)
        ax = Axis(fig[1, 1]; titlesize=10, aspect=DataAspect(), axis_kw...)
        if N > 0
            ax.xlabel = label(reciprocalspace, 1)
        end
        if N > 1
            ax.ylabel = label(reciprocalspace, 2)
        end
    else
        ax = Axis3(fig[1, 1]; titlesize=10, axis_kw...)
        if N > 0
            ax.xlabel = label(reciprocalspace, 1)
        end
        if N > 1
            ax.ylabel = label(reciprocalspace, 2)
        end
        if N > 2
            ax.zlabel = label(reciprocalspace, 3)
        end
    end

    Makie.plot!(ax, reciprocalspace; fractional=fractional, autolims=autolims)
    return fig
end

# 3. ReciprocalPath plotting
function Makie.plot!(ax::Union{Axis, Axis3}, path::ReciprocalPath)
    dim = dimension(path)

    x, y, text_labels = get_reciprocal_labels(path)

    if dim == 2
        isa(ax, Axis) || error("Axis required for 2D path")
        coordinates = [tuple(c...) for c in path]
        x_vals = [c[1] for c in coordinates]
        y_vals = [c[2] for c in coordinates]
        scatter!(ax, x_vals, y_vals; marker=:circle, strokewidth=0)

        # Add annotations
        for (xi, yi, ti) in zip(x, y, text_labels)
            text!(ax, xi, yi; text=ti, fontsize=8)
        end
    else
        isa(ax, Axis3) || error("Axis3 required for 3D path")
        coordinates = [tuple(c...) for c in path]
        x_vals = [c[1] for c in coordinates]
        y_vals = [c[2] for c in coordinates]
        z_vals = [c[3] for c in coordinates]
        scatter!(ax, x_vals, y_vals, z_vals; marker=:circle, strokewidth=0)
    end

    return ax
end

function Makie.plot(path::ReciprocalPath; figure_kw=(;), axis_kw=(;))
    fig = Figure(; figure_kw...)
    dim = dimension(path)

    if dim == 2
        ax = Axis(fig[1, 1]; aspect=DataAspect(), titlesize=10,
                  xlabel=label(path, 1), ylabel=label(path, 2), axis_kw...)
    else
        ax = Axis3(fig[1, 1]; titlesize=10,
                   xlabel=label(path, 1), ylabel=label(path, 2), zlabel=label(path, 3), axis_kw...)
    end

    Makie.plot!(ax, path)
    return fig
end

# 4. ReciprocalCurve plotting
function Makie.plot!(ax::Union{Axis, Axis3}, curve::ReciprocalCurve)
    dim = dimension(curve)

    coordinates = [tuple(c...) for c in curve]

    if dim == 2
        isa(ax, Axis) || error("Axis required for 2D curve")
        x_vals = [c[1] for c in coordinates]
        y_vals = [c[2] for c in coordinates]
        scatter!(ax, x_vals, y_vals; marker=:circle, strokewidth=0)
    else
        isa(ax, Axis3) || error("Axis3 required for 3D curve")
        x_vals = [c[1] for c in coordinates]
        y_vals = [c[2] for c in coordinates]
        z_vals = [c[3] for c in coordinates]
        scatter!(ax, x_vals, y_vals, z_vals; marker=:circle, strokewidth=0)
    end

    return ax
end

function Makie.plot(curve::ReciprocalCurve; figure_kw=(;), axis_kw=(;))
    fig = Figure(; figure_kw...)
    dim = dimension(curve)

    if dim == 2
        ax = Axis(fig[1, 1]; aspect=DataAspect(), titlesize=10,
                  xlabel=label(curve, 1), ylabel=label(curve, 2), axis_kw...)
    else
        ax = Axis3(fig[1, 1]; titlesize=10,
                   xlabel=label(curve, 1), ylabel=label(curve, 2), zlabel=label(curve, 3), axis_kw...)
    end

    Makie.plot!(ax, curve)
    return fig
end

# 5. Heatmap for BrillouinZone/ReciprocalZone with data
function Makie.plot(reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractMatrix{<:Number}; figure_kw=(;), axis_kw=(;))
    @assert length(reciprocalspace.reciprocals)==2 "plot error: only two dimensional reciprocal spaces are supported."

    fig = Figure(; figure_kw...)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), titlesize=10, axis_kw...)

    x = collect(range(reciprocalspace, 1))
    y = collect(range(reciprocalspace, 2))

    ax.xlabel = label(reciprocalspace, 1)
    ax.ylabel = label(reciprocalspace, 2)

    hm = heatmap!(ax, x, y, transpose(data); colorrange=extrema(data))
    Colorbar(fig[1, 2], hm)

    return fig
end

# 6. ReciprocalScatter with weights
function Makie.plot(reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}; fractional=true, weightmultiplier=1.0, weightcolors=nothing, weightlabels=nothing, figure_kw=(;), axis_kw=(;))
    fig = Figure(; figure_kw...)

    if fractional
        N = rank(reciprocalscatter)
        coords_list = fractionals(reciprocalscatter)
        coordinates = [tuple(c...) for c in coords_list]
    else
        coordinates = [tuple(c...) for c in reciprocalscatter]
    end

    dim = length(first(coordinates))

    if dim == 2
        ax = Axis(fig[1, 1]; aspect=DataAspect(), titlesize=10, axis_kw...)
        x_vals = [c[1] for c in coordinates]
        y_vals = [c[2] for c in coordinates]

        # Legend entries for weight colors
        show_legend = !isnothing(weightlabels)

        for (i, index) in enumerate(axes(weights, 2))
            color = isnothing(weightcolors) ? Makie.wong_colors()[mod1(i, 7)] : weightcolors[i]
            sizes = weights[:, index] .* weightmultiplier .+ atol

            # Add legend entry
            if show_legend
                scatter!(ax, [x_vals[1]], [y_vals[1]]; markersize=5, color=color, label=weightlabels[i])
            end

            scatter!(ax, x_vals, y_vals; markersize=sizes, color=color, strokewidth=0)
        end

        if show_legend
            axislegend(ax)
        end
    end

    return fig
end

# 7. Line plot for ReciprocalPath with data
function Makie.plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}; figure_kw=(;), axis_kw=(;))
    fig = Figure(; figure_kw...)
    ax = Axis(fig[1, 1]; titlesize=10, xminorticksvisible=true, yminorticksvisible=true,
              xminorticks=Makie.IntervalsBetween(10), yminorticks=Makie.IntervalsBetween(10),
              xgridvisible=true, ygridvisible=true, xminorgridvisible=true, yminorgridvisible=true,
              axis_kw...)

    positions, labels = ticks(path)
    ax.xlabel = string(label(path))
    ax.xticks = (positions, labels)
    Makie.xlims!(ax, 0, distance(path))

    x_vals = [distance(path, i) for i=1:length(path)]

    for j in 1:size(data, 2)
        lines!(ax, x_vals, data[:, j])
    end

    return fig
end

# 8. Scatter plot for ReciprocalPath with data and weights
function Makie.plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3};
                    weightmultiplier=5.0, weightcolors=nothing, weightlabels=nothing, figure_kw=(;), axis_kw=(;))
    fig = Figure(; figure_kw...)
    ax = Axis(fig[1, 1]; titlesize=10,
              xminorticksvisible=true, yminorticksvisible=true,
              xminorticks=Makie.IntervalsBetween(10), yminorticks=Makie.IntervalsBetween(10),
              axis_kw...)

    positions, labels = ticks(path)
    ax.xlabel = string(label(path))
    ax.xticks = (positions, labels)
    Makie.xlims!(ax, 0, distance(path))

    x_vals = [distance(path, i) for i=1:length(path)]
    show_legend = !isnothing(weightlabels)

    # Add legend entries
    if show_legend
        for i in 1:size(weights, 3)
            color = isnothing(weightcolors) ? Makie.wong_colors()[mod1(i, 7)] : weightcolors[i]
            scatter!(ax, [0.0], [0.0]; markersize=5, color=color, label=weightlabels[i])
        end
    end

    # Plot data points with weights
    for i in 1:size(weights, 3)
        color = isnothing(weightcolors) ? Makie.wong_colors()[mod1(i, 7)] : weightcolors[i]
        for j in 1:size(data, 2)
            y_vals = data[:, j]
            w_vals = weights[:, j, i] .* weightmultiplier .+ atol
            scatter!(ax, x_vals, y_vals; markersize=w_vals, color=color, alpha=0.8, strokewidth=0)
        end
    end

    if show_legend
        axislegend(ax)
    end

    return fig
end

# 9. Heatmap for ReciprocalPath with y data
function Makie.plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number}; figure_kw=(;), axis_kw=(;))
    fig = Figure(; figure_kw...)
    ax = Axis(fig[1, 1]; titlesize=10, axis_kw...)

    positions, labels = ticks(path)
    ax.xlabel = string(label(path))
    ax.xticks = (positions, labels)
    Makie.xlims!(ax, 0, distance(path))
    Makie.ylims!(ax, minimum(y), maximum(y))

    x_vals = [distance(path, i) for i=1:length(path)]

    hm = heatmap!(ax, x_vals, y, transpose(data); colorrange=extrema(data))
    Colorbar(fig[1, 2], hm)

    return fig
end

# 10. Multiple heatmaps for 3D data
function Makie.plot(reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractArray{<:Number, 3};
                    subtitles=nothing, subtitlesize=8, nrow=nothing, ncol=nothing, clims=nothing, figure_kw=(;), axis_kw=(;))
    if isnothing(clims)
        clims = extrema(data)
    end
    if isnothing(nrow)
        nrow = round(Int, sqrt(size(data, 3)))
    end
    if isnothing(ncol)
        ncol = ceil(Int, size(data, 3)/nrow)
    end

    fig = Figure(; figure_kw...)

    x = collect(range(reciprocalspace, 1))
    y = collect(range(reciprocalspace, 2))

    for i = 1:size(data, 3)
        row = div(i-1, ncol) + 1
        col = (i-1) % ncol + 1

        ax = Axis(fig[row, col]; aspect=DataAspect(), axis_kw...)
        ax.xlabel = label(reciprocalspace, 1)
        ax.ylabel = label(reciprocalspace, 2)

        if !isnothing(subtitles)
            ax.title = subtitles[i]
            ax.titlesize = subtitlesize
        end

        heatmap!(ax, x, y, transpose(data[:, :, i]); colorrange=clims)
    end

    # Add shared colorbar
    Colorbar(fig[nrow+1, 1:ncol], limits=clims, vertical=false)

    return fig
end

function Makie.plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractArray{<:Number, 3};
                    subtitles=nothing, subtitlesize=8, nrow=nothing, ncol=nothing, clims=nothing, figure_kw=(;), axis_kw=(;))
    if isnothing(clims)
        clims = extrema(data)
    end
    if isnothing(nrow)
        nrow = round(Int, sqrt(size(data, 3)))
    end
    if isnothing(ncol)
        ncol = ceil(Int, size(data, 3)/nrow)
    end

    fig = Figure(; figure_kw...)

    x_vals = [distance(path, i) for i=1:length(path)]
    positions, labels = ticks(path)

    for i = 1:size(data, 3)
        row = div(i-1, ncol) + 1
        col = (i-1) % ncol + 1

        ax = Axis(fig[row, col]; axis_kw...)
        ax.xlabel = string(label(path))
        ax.xticks = (positions, labels)

        if !isnothing(subtitles)
            ax.title = subtitles[i]
            ax.titlesize = subtitlesize
        end

        heatmap!(ax, x_vals, y, transpose(data[:, :, i]); colorrange=clims)
    end

    Colorbar(fig[nrow+1, 1:ncol], limits=clims, vertical=false)

    return fig
end

# plot! methods for multi-layer plotting

# plot! for Figure - delegates to the axis
function Makie.plot!(fig::Figure, obj; kwargs...)
    ax = get_axis(fig)
    return Makie.plot!(ax, obj; kwargs...)
end

# plot! for array of tuples/points - adds scatter (used by multi-layer plots)
function Makie.plot!(ax::Union{Axis, Axis3}, points::AbstractVector{<:Tuple}; seriestype=:scatter, kwargs...)
    if seriestype == :scatter
        length(first(points)) == 2 || error("Only 2D points supported for Axis")
        x_vals = [p[1] for p in points]
        y_vals = [p[2] for p in points]
        scatter!(ax, x_vals, y_vals; marker=:circle, strokewidth=0, kwargs...)
    elseif seriestype == :lines
        length(first(points)) == 2 || error("Only 2D points supported for Axis")
        x_vals = [p[1] for p in points]
        y_vals = [p[2] for p in points]
        lines!(ax, x_vals, y_vals; kwargs...)
    end
    return ax
end

# 11. Assignment plotting
function seriestype_for_data(data)
    t = Tuple(data)
    return seriestype(t...)
end

function seriestype(args...)
    return nothing
end

function seriestype(data::Data)
    return seriestype(Tuple(data)...)
end

function seriestype(x::AbstractVector{<:Number}, y::Union{AbstractVector{<:Number}, AbstractMatrix{<:Number}}, args...)
    return :path
end

function seriestype(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, z::Union{AbstractMatrix{<:Number}, AbstractArray{<:Number, 3}}, args...)
    return :heatmap
end

function Makie.plot(assignment::Assignment; figure_kw=(;), axis_kw=(;))
    data = assignment.data
    t = Tuple(data)

    st = seriestype(t...)

    if st == :path
        fig = Figure(; figure_kw...)
        ax = Axis(fig[1, 1]; title=str(assignment), titlesize=10,
                  xminorticksvisible=true, yminorticksvisible=true,
                  xminorticks=Makie.IntervalsBetween(10), yminorticks=Makie.IntervalsBetween(10),
                  xgridvisible=true, ygridvisible=true, xminorgridvisible=true, yminorgridvisible=true,
                  axis_kw...)

        if length(t) >= 2
            x_data = t[1]
            y_data = t[2]

            if isa(y_data, AbstractMatrix)
                for j in 1:size(y_data, 2)
                    lines!(ax, x_data, y_data[:, j])
                end
            else
                lines!(ax, x_data, y_data)
            end
        end
        return fig
    elseif st == :heatmap
        fig = Figure(; figure_kw...)
        ax = Axis(fig[1, 1]; title=str(assignment), titlesize=10, axis_kw...)

        if length(t) >= 3
            x_data = t[1]
            y_data = t[2]
            z_data = t[3]

            hm = heatmap!(ax, x_data, y_data, z_data; colorrange=extrema(z_data))
            Colorbar(fig[1, 2], hm)
        end
        return fig
    else
        # Default: try to plot as lines if it is vector data
        fig = Figure(; figure_kw...)
        ax = Axis(fig[1, 1]; title=str(assignment), titlesize=10, axis_kw...)

        for item in t
            if isa(item, AbstractVector{<:Number})
                lines!(ax, item)
            end
        end
        return fig
    end
end

end # module
