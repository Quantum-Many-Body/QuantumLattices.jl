module QuantumLatticesMakieExt

import Makie: plot, plot!
using Makie: AbstractAxis, Axis, Axis3, Colorbar, DataAspect, Figure, axislegend, IntervalsBetween, limits!, lines!, scatter!, heatmap!, text!, wong_colors, xlims!, ylims!
using QuantumLattices: AbstractLattice, Assignment, Bond, BrillouinZone, FractionalReciprocalSpace, Neighbors, Point, ReciprocalCurve, ReciprocalPath, ReciprocalScatter, ReciprocalZone
using QuantumLattices: atol, bonds, dimension, distance, fractionals, getcontent, label, rank, scalartype, seriestype, str, ticks

# 1. Lattice plotting
function plot(lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; siteon=false, figure_kw=(;), axis_kw=(;))
    fig = Figure(; figure_kw...)
    axis_kw = (; title=string(getcontent(lattice, :name)), titlesize=10, axis_kw...)
    ax = if dimension(lattice) <= 2
        Axis(fig[1, 1]; aspect=DataAspect(), axis_kw...)
    else
        Axis3(fig[1, 1]; axis_kw...)
    end
    plot!(ax, lattice, neighbors, filter; siteon)
    return fig
end
function plot!(ax::AbstractAxis, lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; siteon=false)
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
            scatter!(ax, coords...; marker=:circle, strokewidth=0)
            if siteon
                for (i, label) in enumerate(labels)
                    !isempty(label) && text!(ax, coords[1][i], coords[2][i]; text=label, fontsize=8)
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
function plot(reciprocalspace::FractionalReciprocalSpace; fractional=false, autolims=true, figure_kw=(;), axis_kw=(;))
    N = rank(reciprocalspace)
    fig = Figure(; figure_kw...)
    axis_kw = (titlesize=10, axis_kw...)
    N > 0 && (axis_kw = (xlabel=label(reciprocalspace, 1), axis_kw...))
    N > 1 && (axis_kw = (ylabel=label(reciprocalspace, 2), axis_kw...))
    N > 2 && (axis_kw = (zlabel=label(reciprocalspace, 3), axis_kw...))
    if dimension(reciprocalspace) <= 2 || (fractional && N <= 2)
        ax = Axis(fig[1, 1]; aspect=DataAspect(), axis_kw...)
    else
        ax = Axis3(fig[1, 1]; axis_kw...)
    end
    plot!(ax, reciprocalspace; fractional=fractional, autolims=autolims)
    return fig
end
function plot!(ax::AbstractAxis, reciprocalspace::FractionalReciprocalSpace; fractional=false, autolims=true)
    dim = fractional ? rank(reciprocalspace) : dimension(reciprocalspace)
    dim <= 2 ? (ax isa Axis || error("Axis required for $(dim)D reciprocal space")) : (ax isa Axis3 || error("Axis3 required for 3D reciprocal space"))
    coords = fractional ? fractionals(reciprocalspace) : collect(reciprocalspace)
    coords = [[c[d] for c in coords] for d in 1:dim]
    autolims || limits!(ax, extrema.(coords)...)
    scatter!(ax, coords...; marker=:circle, strokewidth=0)
    dim == 2 && (ax.xgridvisible = ax.ygridvisible = ax.xminorgridvisible = ax.yminorgridvisible = true)
    return ax
end

# 3. ReciprocalPath plotting
function plot(path::ReciprocalPath; figure_kw=(;), axis_kw=(;))
    dimension(path)==2 || error("Only 2D reciprocal path can be plotted.")
    fig = Figure(; figure_kw...)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), titlesize=10, xlabel=label(path, 1), ylabel=label(path, 2), axis_kw...)
    plot!(ax, path)
    return fig
end
function plot!(ax::AbstractAxis, path::ReciprocalPath)
    ax isa Axis || error("Axis required for reciprocal path.")
    scatter!(ax, [c[1] for c in path], [c[2] for c in path]; marker=:circle, strokewidth=0)
    for i = 1:length(path.contents)
        if length(path.contents[i]) > 0
            text!(ax, path.contents[i][1][1], path.contents[i][1][2]; text=string(path.labels[i].first), fontsize=8)
            j = i%length(path.labels) + 1
            if length(path.contents[i])>1 && path.labels[i].second≠path.labels[j].first
                text!(ax, path.contents[i][end][1], path.contents[i][end][2]; text=string(path.labels[i].second), fontsize=8)
            end
        end
    end
    return ax
end

# 4. ReciprocalCurve plotting
function plot(curve::ReciprocalCurve; figure_kw=(;), axis_kw=(;))
    dimension(curve)==2 || error("Only 2D reciprocal curve can be plotted.")
    fig = Figure(; figure_kw...)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), titlesize=10, xlabel=label(curve, 1), ylabel=label(curve, 2), axis_kw...)
    plot!(ax, curve)
    return fig
end
function plot!(ax::AbstractAxis, curve::ReciprocalCurve)
    ax isa Axis || error("Axis required for reciprocal curve.")
    scatter!(ax, [c[1] for c in curve], [c[2] for c in curve]; marker=:circle, strokewidth=0)
    return ax
end

# 5. Heatmap for BrillouinZone/ReciprocalZone with data
function plot(reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractMatrix{<:Number}; subtitles=nothing, subtitlesize=8, clims=nothing, figure_kw=(;), axis_kw=(;))
    @assert length(reciprocalspace.reciprocals)==2 "plot error: only two dimensional reciprocal spaces are supported."
    fig = Figure(; figure_kw...)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), titlesize=10, axis_kw...)
    hm = plot!(ax, reciprocalspace, data; subplot=1, subtitles=subtitles, subtitlesize=subtitlesize, clims=clims)
    Colorbar(fig[1, 2], hm)
    return fig
end
function plot!(ax::AbstractAxis, reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractMatrix{<:Number}; subplot=1, subtitles=nothing, subtitlesize=8, clims=nothing)
    ax isa Axis || error("Axis required for heatmap")
    x = collect(range(reciprocalspace, 1))
    y = collect(range(reciprocalspace, 2))
    Δx, Δy = x[2]-x[1], y[2]-y[1]
    ax.xlabel = label(reciprocalspace, 1)
    ax.ylabel = label(reciprocalspace, 2)
    xlims!(ax, x[1]-Δx, x[end]+Δx)
    ylims!(ax, y[1]-Δy, y[end]+Δy)
    isnothing(subtitles) || (ax.title = subtitles[subplot]; ax.titlesize = subtitlesize)
    isnothing(clims) && (clims = extrema(data))
    return heatmap!(ax, x, y, transpose(data); colorrange=clims)
end

# 6. ReciprocalScatter with weights
function plot(reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}; fractional=true, weightmultiplier=1.0, weightcolors=nothing, weightlabels=nothing, figure_kw=(;), axis_kw=(;))
    fig = Figure(; figure_kw...)
    coordinates = fractional ? [tuple(c...) for c in fractionals(reciprocalscatter)] : [tuple(c...) for c in reciprocalscatter]
    length(first(coordinates)) == 2 || return fig
    ax = Axis(fig[1, 1]; aspect=DataAspect(), titlesize=10, axis_kw...)
    plot!(ax, reciprocalscatter, weights; fractional=fractional, weightmultiplier=weightmultiplier, weightcolors=weightcolors, weightlabels=weightlabels)
    return fig
end
function plot!(ax::AbstractAxis, reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}; fractional=true, weightmultiplier=1.0, weightcolors=nothing, weightlabels=nothing)
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
        scatter!(ax, x, y; markersize=weights[:, index] .* weightmultiplier .+ atol, color=color, strokewidth=0)
    end
    show_legend && axislegend(ax)
    return ax
end

# 7. Line plot for ReciprocalPath with data
function plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}; figure_kw=(;), axis_kw=(;))
    fig = Figure(; figure_kw...)
    ax = Axis(fig[1, 1]; titlesize=10, xminorticksvisible=true, yminorticksvisible=true, xminorticks=IntervalsBetween(10), yminorticks=IntervalsBetween(10), xgridvisible=true, ygridvisible=true, xminorgridvisible=true, yminorgridvisible=true, axis_kw...)
    plot!(ax, path, data)
    return fig
end
function plot!(ax::AbstractAxis, path::ReciprocalPath, data::AbstractMatrix{<:Number})
    ax isa Axis || error("Axis required for line plot")
    ax.xlabel = string(label(path))
    ax.xticks = ticks(path)
    xlims!(ax, 0, distance(path))
    for j in 1:size(data, 2)
        lines!(ax, [distance(path, i) for i=1:length(path)], data[:, j])
    end
    return ax
end

# 8. Scatter plot for ReciprocalPath with data and weights
function plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}; weightmultiplier=5.0, weightcolors=nothing, weightlabels=nothing, figure_kw=(;), axis_kw=(;))
    fig = Figure(; figure_kw...)
    ax = Axis(fig[1, 1]; titlesize=10, xminorticksvisible=true, yminorticksvisible=true, xminorticks=IntervalsBetween(10), yminorticks=IntervalsBetween(10), axis_kw...)
    plot!(ax, path, data, weights; weightmultiplier=weightmultiplier, weightcolors=weightcolors, weightlabels=weightlabels)
    return fig
end
function plot!(ax::AbstractAxis, path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}; weightmultiplier=5.0, weightcolors=nothing, weightlabels=nothing)
    ax isa Axis || error("Axis required for scatter plot")
    ax.xlabel = string(label(path))
    ax.xticks = ticks(path)
    xlims!(ax, 0, distance(path))
    x = [distance(path, i) for i=1:length(path)]
    show_legend = !isnothing(weightlabels)
    for i in 1:size(weights, 3)
        color = isnothing(weightcolors) ? wong_colors()[mod1(i, 7)] : weightcolors[i]
        show_legend && scatter!(ax, [0.0], [0.0]; markersize=5, color=color, label=weightlabels[i])
        for j in 1:size(data, 2)
            scatter!(ax, x, data[:, j]; markersize=weights[:, j, i] .* weightmultiplier .+ atol, color=color, alpha=0.8, strokewidth=0)
        end
    end
    show_legend && axislegend(ax)
    return ax
end

# 9. Heatmap for ReciprocalPath with y data
function plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number}; subtitles=nothing, subtitlesize=8, clims=nothing, figure_kw=(;), axis_kw=(;))
    fig = Figure(; figure_kw...)
    ax = Axis(fig[1, 1]; titlesize=10, axis_kw...)
    hm = plot!(ax, path, y, data; subplot=1, subtitles=subtitles, subtitlesize=subtitlesize, clims=clims)
    Colorbar(fig[1, 2], hm)
    return fig
end
function plot!(ax::AbstractAxis, path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number}; subplot=1, subtitles=nothing, subtitlesize=8, clims=nothing)
    ax isa Axis || error("Axis required for heatmap")
    ax.xlabel = string(label(path))
    ax.xticks = ticks(path)
    xlims!(ax, 0, distance(path))
    ylims!(ax, minimum(y), maximum(y))
    isnothing(subtitles) || (ax.title = subtitles[subplot]; ax.titlesize = subtitlesize)
    isnothing(clims) && (clims = extrema(data))
    return heatmap!(ax, [distance(path, i) for i=1:length(path)], y, transpose(data); colorrange=clims)
end

# 10. Multiple heatmaps for 3D data
function plot(reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlesize=8, nrow=nothing, ncol=nothing, clims=nothing, figure_kw=(;), axis_kw=(;))
    isnothing(clims) && (clims = extrema(data))
    isnothing(nrow) && (nrow = round(Int, sqrt(size(data, 3))))
    isnothing(ncol) && (ncol = ceil(Int, size(data, 3)/nrow))
    fig = Figure(; figure_kw...)
    x = collect(range(reciprocalspace, 1))
    y = collect(range(reciprocalspace, 2))
    for i = 1:size(data, 3)
        ax = Axis(fig[div(i-1, ncol) + 1, (i-1) % ncol + 1]; aspect=DataAspect(), axis_kw...)
        plot!(ax, reciprocalspace, data[:, :, i]; subplot=i, subtitles=subtitles, subtitlesize=subtitlesize, clims=clims)
    end
    Colorbar(fig[nrow+1, 1:ncol], limits=clims, vertical=false)
    return fig
end
function plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlesize=8, nrow=nothing, ncol=nothing, clims=nothing, figure_kw=(;), axis_kw=(;))
    isnothing(clims) && (clims = extrema(data))
    isnothing(nrow) && (nrow = round(Int, sqrt(size(data, 3))))
    isnothing(ncol) && (ncol = ceil(Int, size(data, 3)/nrow))
    fig = Figure(; figure_kw...)
    for i = 1:size(data, 3)
        ax = Axis(fig[div(i-1, ncol) + 1, (i-1) % ncol + 1]; axis_kw...)
        plot!(ax, path, y, data[:, :, i]; subplot=i, subtitles=subtitles, subtitlesize=subtitlesize, clims=clims)
    end
    Colorbar(fig[nrow+1, 1:ncol], limits=clims, vertical=false)
    return fig
end

# 11. Assignment plotting
function plot(assignment::Assignment; figure_kw=(;), axis_kw=(;))
    fig = Figure(; figure_kw...)
    attr = seriestype(assignment.data)
    attr in (:path, :heatmap) || error("plot error: unsupported seriestype for Assignment plotting, only :path and :heatmap are supported.")
    axis_kw = (title=str(assignment), titlesize=10, axis_kw...)
    if attr == :path
        ax = Axis(fig[1, 1]; xgridvisible=true, ygridvisible=true, xminorgridvisible=true, yminorgridvisible=true, xminorticksvisible=true, yminorticksvisible=true, xminorticks=IntervalsBetween(10), yminorticks=IntervalsBetween(10), axis_kw...)
        data = Tuple(assignment.data)
        for j in 1:size(data[2], 2)
            lines!(ax, data[1], data[2][:, j])
        end
    else
        ax = Axis(fig[1, 1]; axis_kw...)
        data = Tuple(assignment.data)
        hm = heatmap!(ax, data[1], data[2], transpose(data[3]); colorrange=extrema(data[3]))
        Colorbar(fig[1, 2], hm)
    end
    return fig
end

end # module
