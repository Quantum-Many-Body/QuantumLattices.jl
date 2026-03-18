module QuantumLatticesMakieExt

import Makie: plot, plot!
using Makie: AbstractAxis, Axis, Axis3, Colorbar, DataAspect, Figure, axislegend, IntervalsBetween, limits!, lines!, scatter!, heatmap!, text!, wong_colors, xlims!, ylims!
using QuantumLattices: AbstractLattice, Assignment, Bond, BrillouinZone, FractionalReciprocalSpace, Neighbors, Point, ReciprocalCurve, ReciprocalPath, ReciprocalScatter, ReciprocalZone
using QuantumLattices: atol, bonds, dimension, distance, fractionals, getcontent, label, rank, scalartype, seriestype, str, ticks

# 1. Lattice plotting
function plot(lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; siteon=false, kwargs...)
    fig = Figure()
    kwargs = (title=string(getcontent(lattice, :name)), titlesize=10, kwargs...)
    ax = if dimension(lattice) <= 2
        Axis(fig[1, 1]; aspect=DataAspect(), kwargs...)
    else
        Axis3(fig[1, 1]; kwargs...)
    end
    plot!(ax, lattice, neighbors, filter; siteon, kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; siteon=false, kwargs...)
    dim = dimension(lattice)
    if isa(neighbors, Int) || 0 in keys(neighbors)
        coords = NTuple{dim, scalartype(lattice)}[]
        labels = String[]
        for i in eachindex(lattice)
            bond = Bond(Point(i, lattice[i], zero(lattice[i])))
            if filter(bond)
                push!(coords, ntuple(j -> lattice[i][j], dim))
                push!(labels, siteon ? string(i) : "")
            end
        end
        if length(coords) > 0
            coords = [[c[d] for c in coords] for d in 1:dim]
            scatter!(ax, coords...; marker=get(kwargs, :marker, :circle), strokewidth=get(kwargs, :strokewidth, 0))
            if siteon
                for (i, label) in enumerate(labels)
                    !isempty(label) && text!(ax, coords[1][i], coords[2][i]; text=label, fontsize=get(kwargs, :fontsize, 8))
                end
            end
        end
    end
    for bond in bonds(lattice, neighbors)
        length(bond)==2 && filter(bond) || continue
        coords = [[bond[1].rcoordinate[d], bond[2].rcoordinate[d]] for d in 1:dim]
        lines!(ax, coords...)
    end
    return ax
end

# 2. FractionalReciprocalSpace plotting (BrillouinZone, ReciprocalZone, ReciprocalScatter)
function plot(reciprocalspace::FractionalReciprocalSpace; fractional=false, autolims=true, kwargs...)
    N = rank(reciprocalspace)
    fig = Figure()
    kwargs = (titlesize=10, kwargs...)
    N > 0 && (kwargs = (xlabel=label(reciprocalspace, 1), kwargs...))
    N > 1 && (kwargs = (ylabel=label(reciprocalspace, 2), kwargs...))
    N > 2 && (kwargs = (zlabel=label(reciprocalspace, 3), kwargs...))
    if dimension(reciprocalspace) <= 2 || (fractional && N <= 2)
        ax = Axis(fig[1, 1]; aspect=DataAspect(), kwargs...)
    else
        ax = Axis3(fig[1, 1]; kwargs...)
    end
    plot!(ax, reciprocalspace; fractional=fractional, autolims=autolims, kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, reciprocalspace::FractionalReciprocalSpace; fractional=false, autolims=true, kwargs...)
    dim = fractional ? rank(reciprocalspace) : dimension(reciprocalspace)
    dim <= 2 ? (ax isa Axis || error("Axis required for $(dim)D reciprocal space")) : (ax isa Axis3 || error("Axis3 required for 3D reciprocal space"))
    coords = fractional ? fractionals(reciprocalspace) : collect(reciprocalspace)
    coords = [[c[d] for c in coords] for d in 1:dim]
    autolims || limits!(ax, extrema.(coords)...)
    scatter!(ax, coords...; marker=get(kwargs, :marker, :circle), strokewidth=get(kwargs, :strokewidth, 0))
    if dim == 2
        ax.xgridvisible = get(kwargs, :xgridvisible, true)
        ax.ygridvisible = get(kwargs, :ygridvisible, true)
        ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
        ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
    end
    for (name, default) in zip((:xlabel, :ylabel, :zlabel), ntuple(i->label(reciprocalspace, i), dim))
        setproperty!(ax, name, get(kwargs, name, default))
    end
    return ax
end

# 3. ReciprocalPath plotting
function plot(path::ReciprocalPath; kwargs...)
    dimension(path)==2 || error("Only 2D reciprocal path can be plotted.")
    fig = Figure()
    ax = Axis(fig[1, 1]; aspect=DataAspect(), titlesize=10, xlabel=label(path, 1), ylabel=label(path, 2), kwargs...)
    plot!(ax, path; kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, path::ReciprocalPath; kwargs...)
    ax isa Axis || error("Axis required for reciprocal path.")
    fontsize = get(kwargs, :fontsize, 8)
    scatter!(ax, [c[1] for c in path], [c[2] for c in path]; marker=get(kwargs, :marker, :circle), strokewidth=get(kwargs, :strokewidth, 0))
    for i = 1:length(path.contents)
        if length(path.contents[i]) > 0
            text!(ax, path.contents[i][1][1], path.contents[i][1][2]; text=string(path.labels[i].first), fontsize=fontsize)
            j = i%length(path.labels) + 1
            if length(path.contents[i])>1 && path.labels[i].second≠path.labels[j].first
                text!(ax, path.contents[i][end][1], path.contents[i][end][2]; text=string(path.labels[i].second), fontsize=fontsize)
            end
        end
    end
    return ax
end

# 4. ReciprocalCurve plotting
function plot(curve::ReciprocalCurve; kwargs...)
    dimension(curve)==2 || error("Only 2D reciprocal curve can be plotted.")
    fig = Figure()
    ax = Axis(fig[1, 1]; aspect=DataAspect(), titlesize=10, xlabel=label(curve, 1), ylabel=label(curve, 2), kwargs...)
    plot!(ax, curve; kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, curve::ReciprocalCurve; kwargs...)
    ax isa Axis || error("Axis required for reciprocal curve.")
    scatter!(ax, [c[1] for c in curve], [c[2] for c in curve]; marker=get(kwargs, :marker, :circle), strokewidth=get(kwargs, :strokewidth, 0))
    return ax
end

# 5. Heatmap for BrillouinZone/ReciprocalZone with data
function plot(reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractMatrix{<:Number}; subtitles=nothing, subtitlesize=8, clims=nothing, kwargs...)
    @assert length(reciprocalspace.reciprocals)==2 "plot error: only two dimensional reciprocal spaces are supported."
    fig = Figure()
    ax = Axis(fig[1, 1]; aspect=DataAspect(), titlesize=10, kwargs...)
    hm = plot!(ax, reciprocalspace, data; subplot=1, subtitles=subtitles, subtitlesize=subtitlesize, clims=clims, kwargs...)
    Colorbar(fig[1, 2], hm)
    return fig
end
function plot!(ax::AbstractAxis, reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractMatrix{<:Number}; subplot=1, subtitles=nothing, subtitlesize=8, clims=nothing, kwargs...)
    ax isa Axis || error("Axis required for heatmap")
    x = collect(range(reciprocalspace, 1))
    y = collect(range(reciprocalspace, 2))
    Δx, Δy = x[2]-x[1], y[2]-y[1]
    xlims!(ax, x[1]-Δx, x[end]+Δx)
    ylims!(ax, y[1]-Δy, y[end]+Δy)
    ax.xlabel = get(kwargs, :xlabel, label(reciprocalspace, 1))
    ax.ylabel = get(kwargs, :ylabel, label(reciprocalspace, 2))
    isnothing(subtitles) || (ax.title = subtitles[subplot]; ax.titlesize = subtitlesize)
    isnothing(clims) && (clims = extrema(data))
    return heatmap!(ax, x, y, transpose(data); colorrange=clims)
end

# 6. ReciprocalScatter with weights
function plot(reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}; fractional=true, weightmultiplier=1.0, weightcolors=nothing, weightlabels=nothing, kwargs...)
    fig = Figure()
    coordinates = fractional ? [tuple(c...) for c in fractionals(reciprocalscatter)] : [tuple(c...) for c in reciprocalscatter]
    length(first(coordinates)) == 2 || return fig
    ax = Axis(fig[1, 1]; aspect=DataAspect(), titlesize=10, kwargs...)
    plot!(ax, reciprocalscatter, weights; fractional=fractional, weightmultiplier=weightmultiplier, weightcolors=weightcolors, weightlabels=weightlabels, kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}; fractional=true, weightmultiplier=1.0, weightcolors=nothing, weightlabels=nothing, kwargs...)
    ax isa Axis || error("Axis required for reciprocal scatter")
    coordinates = fractional ? [tuple(c...) for c in fractionals(reciprocalscatter)] : [tuple(c...) for c in reciprocalscatter]
    length(first(coordinates)) == 2 || return ax
    x, y = [c[1] for c in coordinates], [c[2] for c in coordinates]
    xlims!(ax, extrema(x)...)
    ylims!(ax, extrema(y)...)
    show_legend = !isnothing(weightlabels)
    for (i, index) in enumerate(axes(weights, 2))
        color = isnothing(weightcolors) ? wong_colors()[mod1(i, 7)] : weightcolors[i]
        show_legend && scatter!(ax, [x[1]], [y[1]]; markersize=5, color=color, label=weightlabels[i])
        scatter!(ax, x, y; markersize=weights[:, index] .* weightmultiplier .+ atol, color=color, strokewidth=get(kwargs, :strokewidth, 0))
    end
    show_legend && axislegend(ax)
    return ax
end

# 7. Line plot for ReciprocalPath with data
function plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1];
        titlesize=10,
        xminorticksvisible=true,
        yminorticksvisible=true,
        xminorticks=IntervalsBetween(10),
        yminorticks=IntervalsBetween(10),
        xgridvisible=true,
        ygridvisible=true,
        xminorgridvisible=true,
        yminorgridvisible=true,
        kwargs...)
    plot!(ax, path, data; kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, path::ReciprocalPath, data::AbstractMatrix{<:Number}; kwargs...)
    ax isa Axis || error("Axis required for line plot")
    ax.xlabel = get(kwargs, :xlabel, string(label(path)))
    ax.xticks = get(kwargs, :xticks, ticks(path))
    xlims!(ax, 0, distance(path))
    for j in 1:size(data, 2)
        lines!(ax, [distance(path, i) for i=1:length(path)], data[:, j]; linewidth=get(kwargs, :linewidth, 1))
    end
    return ax
end

# 8. Scatter plot for ReciprocalPath with data and weights
function plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}; weightmultiplier=5.0, weightcolors=nothing, weightlabels=nothing, kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1];
        titlesize=10,
        xminorticksvisible=true,
        yminorticksvisible=true,
        xminorticks=IntervalsBetween(10),
        yminorticks=IntervalsBetween(10),
        kwargs...)
    plot!(ax, path, data, weights; weightmultiplier=weightmultiplier, weightcolors=weightcolors, weightlabels=weightlabels, kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}; weightmultiplier=5.0, weightcolors=nothing, weightlabels=nothing, kwargs...)
    ax isa Axis || error("Axis required for scatter plot")
    ax.xlabel = get(kwargs, :xlabel, string(label(path)))
    ax.xticks = get(kwargs, :xticks, ticks(path))
    xlims!(ax, 0, distance(path))
    x = [distance(path, i) for i=1:length(path)]
    show_legend = !isnothing(weightlabels)
    for i in 1:size(weights, 3)
        color = isnothing(weightcolors) ? wong_colors()[mod1(i, 7)] : weightcolors[i]
        show_legend && scatter!(ax, [0.0], [0.0]; markersize=5, color=color, label=weightlabels[i])
        for j in 1:size(data, 2)
            scatter!(ax, x, data[:, j]; markersize=weights[:, j, i] .* weightmultiplier .+ atol, color=color, alpha=get(kwargs, :alpha, 0.8), strokewidth=get(kwargs, :strokewidth, 0))
        end
    end
    show_legend && axislegend(ax)
    return ax
end

# 9. Heatmap for ReciprocalPath with y data
function plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number}; subtitles=nothing, subtitlesize=8, clims=nothing, kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1]; titlesize=10, kwargs...)
    hm = plot!(ax, path, y, data; subplot=1, subtitles=subtitles, subtitlesize=subtitlesize, clims=clims, kwargs...)
    Colorbar(fig[1, 2], hm)
    return fig
end
function plot!(ax::AbstractAxis, path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number}; subplot=1, subtitles=nothing, subtitlesize=8, clims=nothing, kwargs...)
    ax isa Axis || error("Axis required for heatmap")
    ax.xlabel = get(kwargs, :xlabel, string(label(path)))
    ax.xticks = get(kwargs, :xticks, ticks(path))
    xlims!(ax, 0, distance(path))
    ylims!(ax, minimum(y), maximum(y))
    isnothing(subtitles) || (ax.title = subtitles[subplot]; ax.titlesize = subtitlesize)
    isnothing(clims) && (clims = extrema(data))
    return heatmap!(ax, [distance(path, i) for i=1:length(path)], y, transpose(data); colorrange=clims)
end

# 10. Multiple heatmaps for 3D data
function plot(reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlesize=8, nrow=nothing, ncol=nothing, clims=nothing, kwargs...)
    isnothing(clims) && (clims = extrema(data))
    isnothing(nrow) && (nrow = round(Int, sqrt(size(data, 3))))
    isnothing(ncol) && (ncol = ceil(Int, size(data, 3)/nrow))
    fig = Figure()
    x = collect(range(reciprocalspace, 1))
    y = collect(range(reciprocalspace, 2))
    for i = 1:size(data, 3)
        ax = Axis(fig[div(i-1, ncol) + 1, (i-1) % ncol + 1]; aspect=DataAspect(), kwargs...)
        plot!(ax, reciprocalspace, data[:, :, i]; subplot=i, subtitles=subtitles, subtitlesize=subtitlesize, clims=clims, kwargs...)
    end
    Colorbar(fig[nrow+1, 1:ncol], limits=clims, vertical=false)
    return fig
end
function plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlesize=8, nrow=nothing, ncol=nothing, clims=nothing, kwargs...)
    isnothing(clims) && (clims = extrema(data))
    isnothing(nrow) && (nrow = round(Int, sqrt(size(data, 3))))
    isnothing(ncol) && (ncol = ceil(Int, size(data, 3)/nrow))
    fig = Figure()
    for i = 1:size(data, 3)
        ax = Axis(fig[div(i-1, ncol) + 1, (i-1) % ncol + 1]; kwargs...)
        plot!(ax, path, y, data[:, :, i]; subplot=i, subtitles=subtitles, subtitlesize=subtitlesize, clims=clims, kwargs...)
    end
    Colorbar(fig[nrow+1, 1:ncol], limits=clims, vertical=false)
    return fig
end

# 11. Assignment plotting
function plot(assignment::Assignment; kwargs...)
    fig = Figure()
    attr = seriestype(assignment.data)
    attr in (:path, :heatmap) || error("plot error: unsupported seriestype for Assignment plotting, only :path and :heatmap are supported.")
    kwargs = (title=str(assignment), titlesize=10, kwargs...)
    if attr == :path
        ax = Axis(fig[1, 1];
            xgridvisible=true,
            ygridvisible=true,
            xminorgridvisible=true,
            yminorgridvisible=true,
            xminorticksvisible=true,
            yminorticksvisible=true,
            xminorticks=IntervalsBetween(10),
            yminorticks=IntervalsBetween(10),
            kwargs...)
        data = Tuple(assignment.data)
        for j in 1:size(data[2], 2)
            lines!(ax, data[1], data[2][:, j]; linewidth=get(kwargs, :linewidth, 1))
        end
    else
        ax = Axis(fig[1, 1]; kwargs...)
        data = Tuple(assignment.data)
        hm = heatmap!(ax, data[1], data[2], transpose(data[3]); colorrange=get(kwargs, :colorrange, extrema(data[3])))
        Colorbar(fig[1, 2], hm)
    end
    return fig
end

end # module
