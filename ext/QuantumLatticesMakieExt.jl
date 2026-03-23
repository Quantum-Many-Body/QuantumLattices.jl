module QuantumLatticesMakieExt

using Makie: AbstractAxis, Axis, Axis3, Colorbar, DataAspect, Figure, Heatmap, IntervalsBetween, Lines, Scatter, Series, attribute_names, axislegend, heatmap!, limits!, lines!, scatter!, series!, text!, wong_colors, xlims!, ylims!
using QuantumLattices: AbstractLattice, Assignment, Bond, BrillouinZone, FractionalReciprocalSpace, Neighbors, Point, ReciprocalCurve, ReciprocalPath, ReciprocalScatter, ReciprocalZone
using QuantumLattices: atol, bonds, dimension, distance, fractionals, getcontent, label, rank, scalartype, str, ticks
import Makie: plot, plot!

# Helper function to filter kwargs for plot functions
_filter_kwargs(kwargs, plot_type::Type) = NamedTuple(k => v for (k, v) in kwargs if k in attribute_names(plot_type))

# 1. Lattice plotting
@inline plot(lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; kwargs...) = plot!(Figure(), lattice, neighbors, filter; kwargs...)
function plot!(fig::Figure, lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; kwargs...)
    ax = dimension(lattice) <= 2 ? Axis(fig[1, 1]) : Axis3(fig[1, 1])
    plot!(ax, lattice, neighbors, filter; kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; siteon::Bool=false, kwargs...)
    ax.title = get(kwargs, :title, string(getcontent(lattice, :name)))
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax isa Axis && (ax.aspect = get(kwargs, :aspect, DataAspect()))
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
            scatter!(ax, coords...; marker=:circle, strokewidth=0, _filter_kwargs(kwargs, Scatter)...)
            if siteon
                for (i, label) in enumerate(labels)
                    isempty(label) || text!(ax, coords[1][i], coords[2][i]; text=label, fontsize=get(kwargs, :fontsize, 8))
                end
            end
        end
    end
    for bond in bonds(lattice, neighbors)
        length(bond)==2 && filter(bond) || continue
        coords = [[bond[1].rcoordinate[d], bond[2].rcoordinate[d]] for d in 1:dim]
        lines!(ax, coords...; _filter_kwargs(kwargs, Lines)...)
    end
    return ax
end

# 2. FractionalReciprocalSpace plotting (BrillouinZone, ReciprocalZone, ReciprocalScatter)
@inline plot(reciprocalspace::FractionalReciprocalSpace; kwargs...) = plot!(Figure(), reciprocalspace; kwargs...)
function plot!(fig::Figure, reciprocalspace::FractionalReciprocalSpace; fractional::Bool=false, kwargs...)
    ax = (dimension(reciprocalspace) <= 2 || (fractional && rank(reciprocalspace) <= 2)) ? Axis(fig[1, 1]) : Axis3(fig[1, 1])
    plot!(ax, reciprocalspace; fractional, kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, reciprocalspace::FractionalReciprocalSpace; fractional::Bool=false, autolims::Bool=true, kwargs...)
    dim = fractional ? rank(reciprocalspace) : dimension(reciprocalspace)
    dim <= 2 ? (ax isa Axis || error("Axis required for $(dim)D reciprocal space")) : (ax isa Axis3 || error("Axis3 required for 3D reciprocal space"))
    ax.titlesize = get(kwargs, :titlesize, 10)
    if dim == 2
        ax.aspect = get(kwargs, :aspect, DataAspect())
        ax.xgridvisible = get(kwargs, :xgridvisible, true)
        ax.ygridvisible = get(kwargs, :ygridvisible, true)
        ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
        ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
        ax.xminorticksvisible = get(kwargs, :xminorticksvisible, true)
        ax.yminorticksvisible = get(kwargs, :yminorticksvisible, true)
        ax.xminorticks = get(kwargs, :xminorticks, IntervalsBetween(10))
        ax.yminorticks = get(kwargs, :yminorticks, IntervalsBetween(10))
    end
    for (name, default) in zip((:xlabel, :ylabel, :zlabel), ntuple(i->label(reciprocalspace, i), dim))
        setproperty!(ax, name, get(kwargs, name, default))
    end
    coords = fractional ? fractionals(reciprocalspace) : collect(reciprocalspace)
    coords = [[c[d] for c in coords] for d in 1:dim]
    autolims || limits!(ax, extrema.(coords)...)
    scatter!(ax, coords...; marker=:circle, strokewidth=0, _filter_kwargs(kwargs, Scatter)...)
    return ax
end

# 3. ReciprocalPath plotting
@inline plot(path::ReciprocalPath; kwargs...) = plot!(Figure(), path; kwargs...)
function plot!(fig::Figure, path::ReciprocalPath; kwargs...)
    plot!(Axis(fig[1, 1]), path; kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, path::ReciprocalPath; kwargs...)
    dimension(path) == 2 || error("Only 2D reciprocal path can be plotted.")
    ax.aspect = get(kwargs, :aspect, DataAspect())
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xlabel = get(kwargs, :xlabel, label(path, 1))
    ax.ylabel = get(kwargs, :ylabel, label(path, 2))
    ax.xgridvisible = get(kwargs, :xgridvisible, true)
    ax.ygridvisible = get(kwargs, :ygridvisible, true)
    ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
    ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
    ax.xminorticksvisible = get(kwargs, :xminorticksvisible, true)
    ax.yminorticksvisible = get(kwargs, :yminorticksvisible, true)
    ax.xminorticks = get(kwargs, :xminorticks, IntervalsBetween(10))
    ax.yminorticks = get(kwargs, :yminorticks, IntervalsBetween(10))
    scatter!(ax, [c[1] for c in path], [c[2] for c in path]; marker=:circle, strokewidth=0, _filter_kwargs(kwargs, Scatter)...)
    for i = 1:length(path.contents)
        if length(path.contents[i]) > 0
            text!(ax, path.contents[i][1][1], path.contents[i][1][2]; text=string(path.labels[i].first), fontsize=get(kwargs, :fontsize, 8))
            j = i%length(path.labels) + 1
            if length(path.contents[i])>1 && path.labels[i].second≠path.labels[j].first
                text!(ax, path.contents[i][end][1], path.contents[i][end][2]; text=string(path.labels[i].second), fontsize=get(kwargs, :fontsize, 8))
            end
        end
    end
    return ax
end

# 4. ReciprocalCurve plotting
@inline plot(curve::ReciprocalCurve; kwargs...) = plot!(Figure(), curve; kwargs...)
function plot!(fig::Figure, curve::ReciprocalCurve; kwargs...)
    plot!(Axis(fig[1, 1]), curve; kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, curve::ReciprocalCurve; kwargs...)
    dimension(curve) == 2 || error("Only 2D reciprocal curve can be plotted.")
    ax.aspect = get(kwargs, :aspect, DataAspect())
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xlabel = get(kwargs, :xlabel, label(curve, 1))
    ax.ylabel = get(kwargs, :ylabel, label(curve, 2))
    ax.xgridvisible = get(kwargs, :xgridvisible, true)
    ax.ygridvisible = get(kwargs, :ygridvisible, true)
    ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
    ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
    ax.xminorticksvisible = get(kwargs, :xminorticksvisible, true)
    ax.yminorticksvisible = get(kwargs, :yminorticksvisible, true)
    ax.xminorticks = get(kwargs, :xminorticks, IntervalsBetween(10))
    ax.yminorticks = get(kwargs, :yminorticks, IntervalsBetween(10))
    scatter!(ax, [c[1] for c in curve], [c[2] for c in curve]; marker=:circle, strokewidth=0, _filter_kwargs(kwargs, Scatter)...)
    return ax
end

# 5. Heatmap for BrillouinZone/ReciprocalZone with data
@inline plot(reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractMatrix{<:Number}; kwargs...) = plot!(Figure(), reciprocalspace, data; kwargs...)
function plot!(fig::Figure, reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractMatrix{<:Number}; kwargs...)
    ax = Axis(fig[1, 1])
    plot!(ax, reciprocalspace, data; kwargs...)
    Colorbar(fig[1, 2], ax.scene.plots[end])
    return fig
end
function plot!(ax::AbstractAxis, reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractMatrix{<:Number}; clims=nothing, kwargs...)
    length(reciprocalspace.reciprocals)==2 || error("Only two dimensional reciprocal spaces are supported.")
    ax.aspect = get(kwargs, :aspect, DataAspect())
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xlabel = get(kwargs, :xlabel, label(reciprocalspace, 1))
    ax.ylabel = get(kwargs, :ylabel, label(reciprocalspace, 2))
    x = collect(range(reciprocalspace, 1))
    y = collect(range(reciprocalspace, 2))
    Δx, Δy = x[2]-x[1], y[2]-y[1]
    xlims!(ax, get(kwargs, :xlims, (x[1]-Δx, x[end]+Δx))...)
    ylims!(ax, get(kwargs, :ylims, (y[1]-Δy, y[end]+Δy))...)
    isnothing(clims) && (clims = extrema(data))
    heatmap!(ax, x, y, transpose(data); colorrange=clims, _filter_kwargs(kwargs, Heatmap)...)
    return ax
end

# 6. ReciprocalScatter with weights
@inline plot(reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}; kwargs...) = plot!(Figure(), reciprocalscatter, weights; kwargs...)
function plot!(fig::Figure, reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}; kwargs...)
    plot!(Axis(fig[1, 1]), reciprocalscatter, weights; kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}; fractional::Bool=true, weightmultiplier::Real=5.0, weightcolors=nothing, weightlabels=nothing, kwargs...)
    ax isa Axis || error("Axis required for reciprocal scatter")
    ax.aspect = get(kwargs, :aspect, DataAspect())
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xgridvisible = get(kwargs, :xgridvisible, true)
    ax.ygridvisible = get(kwargs, :ygridvisible, true)
    ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
    ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
    ax.xminorticksvisible = get(kwargs, :xminorticksvisible, true)
    ax.yminorticksvisible = get(kwargs, :yminorticksvisible, true)
    ax.xminorticks = get(kwargs, :xminorticks, IntervalsBetween(10))
    ax.yminorticks = get(kwargs, :yminorticks, IntervalsBetween(10))
    coordinates = fractional ? [tuple(c...) for c in fractionals(reciprocalscatter)] : [tuple(c...) for c in reciprocalscatter]
    length(first(coordinates)) == 2 || error("ReciprocalScatter with weights requires 2D coordinates")
    x, y = [c[1] for c in coordinates], [c[2] for c in coordinates]
    xlims!(ax, get(kwargs, :xlims, extrema(x))...)
    ylims!(ax, get(kwargs, :ylims, extrema(y))...)
    for (i, index) in enumerate(axes(weights, 2))
        color = isnothing(weightcolors) ? wong_colors()[mod1(i, 7)] : weightcolors[i]
        isnothing(weightlabels) || scatter!(ax, [x[1]], [y[1]]; markersize=5, color=color, label=weightlabels[i])
        scatter!(ax, x, y; markersize=weights[:, index] .* weightmultiplier .+ atol, color=color, strokewidth=0, _filter_kwargs(kwargs, Scatter)...)
    end
    isnothing(weightlabels) || axislegend(ax)
    return ax
end

# 7. Line plot for ReciprocalPath with data
@inline plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}; kwargs...) = plot!(Figure(), path, data; kwargs...)
function plot!(fig::Figure, path::ReciprocalPath, data::AbstractMatrix{<:Number}; kwargs...)
    plot!(Axis(fig[1, 1]), path, data; kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, path::ReciprocalPath, data::AbstractMatrix{<:Number}; kwargs...)
    ax isa Axis || error("Axis required for line plot")
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xlabel = get(kwargs, :xlabel, string(label(path)))
    ax.xticks = get(kwargs, :xticks, ticks(path))
    ax.xminorticksvisible = get(kwargs, :xminorticksvisible, true)
    ax.yminorticksvisible = get(kwargs, :yminorticksvisible, true)
    ax.xminorticks = get(kwargs, :xminorticks, IntervalsBetween(10))
    ax.yminorticks = get(kwargs, :yminorticks, IntervalsBetween(10))
    ax.xgridvisible = get(kwargs, :xgridvisible, true)
    ax.ygridvisible = get(kwargs, :ygridvisible, true)
    ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
    ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
    xlims!(ax, get(kwargs, :xlims, (0, distance(path)))...)
    series!(ax, [distance(path, i) for i=1:length(path)], transpose(data); linewidth=1, _filter_kwargs(kwargs, Series)...)
    return ax
end

# 8. Scatter plot for ReciprocalPath with data and weights
@inline plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}; kwargs...) = plot!(Figure(), path, data, weights; kwargs...)
function plot!(fig::Figure, path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}; kwargs...)
    plot!(Axis(fig[1, 1]), path, data, weights; kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}; weightmultiplier::Real=5.0, weightcolors=nothing, weightlabels=nothing, kwargs...)
    ax isa Axis || error("Axis required for scatter plot")
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xlabel = get(kwargs, :xlabel, string(label(path)))
    ax.xticks = get(kwargs, :xticks, ticks(path))
    ax.xgridvisible = get(kwargs, :xgridvisible, true)
    ax.ygridvisible = get(kwargs, :ygridvisible, true)
    ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
    ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
    ax.xminorticksvisible = get(kwargs, :xminorticksvisible, true)
    ax.yminorticksvisible = get(kwargs, :yminorticksvisible, true)
    ax.xminorticks = get(kwargs, :xminorticks, IntervalsBetween(10))
    ax.yminorticks = get(kwargs, :yminorticks, IntervalsBetween(10))
    xlims!(ax, get(kwargs, :xlims, (0, distance(path)))...)
    x = [distance(path, i) for i=1:length(path)]
    show_legend = !isnothing(weightlabels)
    for i in 1:size(weights, 3)
        color = isnothing(weightcolors) ? wong_colors()[mod1(i, 7)] : weightcolors[i]
        show_legend && scatter!(ax, [0.0], [0.0]; markersize=5, color=color, label=weightlabels[i])
        for j in 1:size(data, 2)
            scatter!(ax, x, data[:, j]; markersize=weights[:, j, i] .* weightmultiplier .+ atol, color=color, alpha=0.8, strokewidth=0, _filter_kwargs(kwargs, Scatter)...)
        end
    end
    show_legend && axislegend(ax)
    return ax
end

# 9. Heatmap for ReciprocalPath with y data
@inline plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number}; kwargs...) = plot!(Figure(), path, y, data; kwargs...)
function plot!(fig::Figure, path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number}; kwargs...)
    ax = Axis(fig[1, 1])
    plot!(ax, path, y, data; kwargs...)
    Colorbar(fig[1, 2], ax.scene.plots[end])
    return fig
end
function plot!(ax::AbstractAxis, path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number}; clims=nothing, kwargs...)
    ax isa Axis || error("Axis required for heatmap")
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xlabel = get(kwargs, :xlabel, string(label(path)))
    ax.xticks = get(kwargs, :xticks, ticks(path))
    xlims!(ax, get(kwargs, :xlims, (0, distance(path)))...)
    ylims!(ax, get(kwargs, :ylims, extrema(y))...)
    isnothing(clims) && (clims = extrema(data))
    heatmap!(ax, [distance(path, i) for i=1:length(path)], y, transpose(data); colorrange=clims, _filter_kwargs(kwargs, Heatmap)...)
    return ax
end

# 10. Multiple heatmaps for 3D data (multi-panel only)
@inline plot(reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractArray{<:Number, 3}; kwargs...) = plot!(Figure(), reciprocalspace, data; kwargs...)
function plot!(fig::Figure, reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractArray{<:Number, 3}; subtitles=nothing, titlesize::Int=8, nrow=nothing, ncol=nothing, clims=nothing, kwargs...)
    isnothing(clims) && (clims = extrema(data))
    isnothing(nrow) && (nrow = round(Int, sqrt(size(data, 3))))
    isnothing(ncol) && (ncol = ceil(Int, size(data, 3)/nrow))
    for i = 1:size(data, 3)
        plot!(Axis(fig[div(i-1, ncol) + 1, (i-1) % ncol + 1]), reciprocalspace, data[:, :, i]; title=isnothing(subtitles) ? nothing : subtitles[i], titlesize=titlesize, clims=clims, kwargs...)
    end
    Colorbar(fig[nrow+1, 1:ncol], limits=clims, vertical=false)
    return fig
end
@inline plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractArray{<:Number, 3}; kwargs...) = plot!(Figure(), path, y, data; kwargs...)
function plot!(fig::Figure, path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractArray{<:Number, 3}; subtitles=nothing, titlesize::Int=8, nrow=nothing, ncol=nothing, clims=nothing, kwargs...)
    isnothing(clims) && (clims = extrema(data))
    isnothing(nrow) && (nrow = round(Int, sqrt(size(data, 3))))
    isnothing(ncol) && (ncol = ceil(Int, size(data, 3)/nrow))
    for i = 1:size(data, 3)
        plot!(Axis(fig[div(i-1, ncol) + 1, (i-1) % ncol + 1]), path, y, data[:, :, i]; title=isnothing(subtitles) ? nothing : subtitles[i], titlesize=titlesize, clims=clims, kwargs...)
    end
    Colorbar(fig[nrow+1, 1:ncol], limits=clims, vertical=false)
    return fig
end

# 11. path plotting for (x, y)
@inline plot(x::AbstractVector{<:Number}, y::AbstractMatrix{<:Number}; kwargs...) = plot!(Figure(), x, y; kwargs...)
function plot!(fig::Figure, x::AbstractVector{<:Number}, y::AbstractMatrix{<:Number}; kwargs...)
    plot!(Axis(fig[1, 1]), x, y; kwargs...)
    return fig
end
function plot!(ax::AbstractAxis, x::AbstractVector{<:Number}, y::AbstractMatrix{<:Number}; kwargs...)
    ax isa Axis || error("Axis required for line plot")
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xgridvisible = get(kwargs, :xgridvisible, true)
    ax.ygridvisible = get(kwargs, :ygridvisible, true)
    ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
    ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
    ax.xminorticksvisible = get(kwargs, :xminorticksvisible, true)
    ax.yminorticksvisible = get(kwargs, :yminorticksvisible, true)
    ax.xminorticks = get(kwargs, :xminorticks, IntervalsBetween(10))
    ax.yminorticks = get(kwargs, :yminorticks, IntervalsBetween(10))
    xlims!(ax, get(kwargs, :xlims, extrema(x))...)
    ylims!(ax, get(kwargs, :ylims, extrema(y))...)
    series!(ax, x, transpose(y); linewidth=1, _filter_kwargs(kwargs, Series)...)
    return ax
end

# 12. Assignment plotting
@inline plot(assignment::Assignment; kwargs...) = plot!(Figure(), assignment; kwargs...)
function plot!(fig::Figure, assignment::Assignment; kwargs...)
    ax = Axis(fig[1, 1])
    plot!(ax, assignment; kwargs...)
    ax.scene.plots[end] isa Heatmap && Colorbar(fig[1, 2], ax.scene.plots[end])
    return fig
end
function plot!(ax::AbstractAxis, assignment::Assignment; kwargs...)
    ax.title = get(kwargs, :title, str(assignment))
    data = Tuple(assignment.data)
    plot!(ax, data...; kwargs...)
    return ax
end

end # module
