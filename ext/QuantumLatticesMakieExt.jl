module QuantumLatticesMakieExt

import Makie
using QuantumLattices: AbstractLattice, Assignment, Bond, BrillouinZone, FractionalReciprocalSpace, Neighbors, Point, ReciprocalCurve, ReciprocalPath, ReciprocalScatter, ReciprocalZone
using QuantumLattices: atol, bonds, dimension, distance, fractionals, getcontent, label, rank, scalartype, str, ticks

# Helper function to filter kwargs for plot functions
_filter_kwargs(kwargs, plot_type::Type) = NamedTuple(k => v for (k, v) in kwargs if k in Makie.attribute_names(plot_type))

# 1. Lattice plotting
@inline Makie.plot(lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; kwargs...) = Makie.plot!(Makie.Figure(), lattice, neighbors, filter; kwargs...)
function Makie.plot!(fig::Makie.Figure, lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; kwargs...)
    ax = dimension(lattice) <= 2 ? Makie.Axis(fig[1, 1]) : Makie.Axis3(fig[1, 1])
    Makie.plot!(ax, lattice, neighbors, filter; kwargs...)
    return fig
end
function Makie.plot!(ax::Makie.AbstractAxis, lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; siteon::Bool=false, kwargs...)
    ax.title = get(kwargs, :title, string(getcontent(lattice, :name)))
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax isa Makie.Axis && (ax.aspect = get(kwargs, :aspect, Makie.DataAspect()))
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
            Makie.scatter!(ax, coords...; marker=:circle, strokewidth=0, _filter_kwargs(kwargs, Makie.Scatter)...)
            if siteon
                for (i, label) in enumerate(labels)
                    isempty(label) || Makie.text!(ax, coords[1][i], coords[2][i]; text=label, fontsize=get(kwargs, :fontsize, 8))
                end
            end
        end
    end
    return Makie.linesegments!(
        ax, [Makie.Point{dim}(point.rcoordinate...) for bond in bonds(lattice, neighbors) if length(bond)==2 && filter(bond) for point in bond];
        _filter_kwargs(kwargs, Makie.LineSegments)...
    )
end

# 2. FractionalReciprocalSpace plotting (BrillouinZone, ReciprocalZone, ReciprocalScatter)
@inline Makie.plot(reciprocalspace::FractionalReciprocalSpace; kwargs...) = Makie.plot!(Makie.Figure(), reciprocalspace; kwargs...)
function Makie.plot!(fig::Makie.Figure, reciprocalspace::FractionalReciprocalSpace; fractional::Bool=false, kwargs...)
    ax = (dimension(reciprocalspace) <= 2 || (fractional && rank(reciprocalspace) <= 2)) ? Makie.Axis(fig[1, 1]) : Makie.Axis3(fig[1, 1])
    Makie.plot!(ax, reciprocalspace; fractional, kwargs...)
    return fig
end
function Makie.plot!(ax::Makie.AbstractAxis, reciprocalspace::FractionalReciprocalSpace; fractional::Bool=false, autolims::Bool=true, kwargs...)
    dim = fractional ? rank(reciprocalspace) : dimension(reciprocalspace)
    dim <= 2 ? (ax isa Makie.Axis || error("Makie.Axis required for $(dim)D reciprocal space")) : (ax isa Makie.Axis3 || error("Makie.Axis3 required for 3D reciprocal space"))
    ax.titlesize = get(kwargs, :titlesize, 10)
    if dim == 2
        ax.aspect = get(kwargs, :aspect, Makie.DataAspect())
        ax.xgridvisible = get(kwargs, :xgridvisible, true)
        ax.ygridvisible = get(kwargs, :ygridvisible, true)
        ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
        ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
        ax.xminorticksvisible = get(kwargs, :xminorticksvisible, true)
        ax.yminorticksvisible = get(kwargs, :yminorticksvisible, true)
        ax.xminorticks = get(kwargs, :xminorticks, Makie.IntervalsBetween(10))
        ax.yminorticks = get(kwargs, :yminorticks, Makie.IntervalsBetween(10))
    end
    for (name, default) in zip((:xlabel, :ylabel, :zlabel), ntuple(i->label(reciprocalspace, i), dim))
        setproperty!(ax, name, get(kwargs, name, default))
    end
    coords = fractional ? fractionals(reciprocalspace) : collect(reciprocalspace)
    coords = [[c[d] for c in coords] for d in 1:dim]
    autolims || Makie.limits!(ax, extrema.(coords)...)
    return Makie.scatter!(ax, coords...; marker=:circle, strokewidth=0, _filter_kwargs(kwargs, Makie.Scatter)...)
end

# 3. ReciprocalPath plotting
@inline Makie.plot(path::ReciprocalPath; kwargs...) = Makie.plot!(Makie.Figure(), path; kwargs...)
function Makie.plot!(fig::Makie.Figure, path::ReciprocalPath; kwargs...)
    Makie.plot!(Makie.Axis(fig[1, 1]), path; kwargs...)
    return fig
end
function Makie.plot!(ax::Makie.AbstractAxis, path::ReciprocalPath; kwargs...)
    dimension(path) == 2 || error("Only 2D reciprocal path can be plotted.")
    ax.aspect = get(kwargs, :aspect, Makie.DataAspect())
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xlabel = get(kwargs, :xlabel, label(path, 1))
    ax.ylabel = get(kwargs, :ylabel, label(path, 2))
    ax.xgridvisible = get(kwargs, :xgridvisible, true)
    ax.ygridvisible = get(kwargs, :ygridvisible, true)
    ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
    ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
    ax.xminorticksvisible = get(kwargs, :xminorticksvisible, true)
    ax.yminorticksvisible = get(kwargs, :yminorticksvisible, true)
    ax.xminorticks = get(kwargs, :xminorticks, Makie.IntervalsBetween(10))
    ax.yminorticks = get(kwargs, :yminorticks, Makie.IntervalsBetween(10))
    for i = 1:length(path.contents)
        if length(path.contents[i]) > 0
            Makie.text!(ax, path.contents[i][1][1], path.contents[i][1][2]; text=string(path.labels[i].first), fontsize=get(kwargs, :fontsize, 8))
            j = i%length(path.labels) + 1
            if length(path.contents[i])>1 && path.labels[i].second≠path.labels[j].first
                Makie.text!(ax, path.contents[i][end][1], path.contents[i][end][2]; text=string(path.labels[i].second), fontsize=get(kwargs, :fontsize, 8))
            end
        end
    end
    return Makie.scatter!(ax, [c[1] for c in path], [c[2] for c in path]; marker=:circle, strokewidth=0, _filter_kwargs(kwargs, Makie.Scatter)...)
end

# 4. ReciprocalCurve plotting
@inline Makie.plot(curve::ReciprocalCurve; kwargs...) = Makie.plot!(Makie.Figure(), curve; kwargs...)
function Makie.plot!(fig::Makie.Figure, curve::ReciprocalCurve; kwargs...)
    Makie.plot!(Makie.Axis(fig[1, 1]), curve; kwargs...)
    return fig
end
function Makie.plot!(ax::Makie.AbstractAxis, curve::ReciprocalCurve; kwargs...)
    dimension(curve) == 2 || error("Only 2D reciprocal curve can be plotted.")
    ax.aspect = get(kwargs, :aspect, Makie.DataAspect())
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xlabel = get(kwargs, :xlabel, label(curve, 1))
    ax.ylabel = get(kwargs, :ylabel, label(curve, 2))
    ax.xgridvisible = get(kwargs, :xgridvisible, true)
    ax.ygridvisible = get(kwargs, :ygridvisible, true)
    ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
    ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
    ax.xminorticksvisible = get(kwargs, :xminorticksvisible, true)
    ax.yminorticksvisible = get(kwargs, :yminorticksvisible, true)
    ax.xminorticks = get(kwargs, :xminorticks, Makie.IntervalsBetween(10))
    ax.yminorticks = get(kwargs, :yminorticks, Makie.IntervalsBetween(10))
    return Makie.scatter!(ax, [c[1] for c in curve], [c[2] for c in curve]; marker=:circle, strokewidth=0, _filter_kwargs(kwargs, Makie.Scatter)...)
end

# 5. Heatmap for BrillouinZone/ReciprocalZone with data
@inline Makie.plot(reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractMatrix{<:Number}; kwargs...) = Makie.plot!(Makie.Figure(), reciprocalspace, data; kwargs...)
function Makie.plot!(fig::Makie.Figure, reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractMatrix{<:Number}; kwargs...)
    Makie.Colorbar(fig[1, 2], Makie.plot!(Makie.Axis(fig[1, 1]), reciprocalspace, data; kwargs...))
    return fig
end
function Makie.plot!(ax::Makie.AbstractAxis, reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractMatrix{<:Number}; clims=nothing, kwargs...)
    length(reciprocalspace.reciprocals)==2 || error("Only two dimensional reciprocal spaces are supported.")
    ax.aspect = get(kwargs, :aspect, Makie.DataAspect())
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xlabel = get(kwargs, :xlabel, label(reciprocalspace, 1))
    ax.ylabel = get(kwargs, :ylabel, label(reciprocalspace, 2))
    x = collect(range(reciprocalspace, 1))
    y = collect(range(reciprocalspace, 2))
    Δx, Δy = x[2]-x[1], y[2]-y[1]
    Makie.xlims!(ax, get(kwargs, :xlims, (x[1]-Δx, x[end]+Δx))...)
    Makie.ylims!(ax, get(kwargs, :ylims, (y[1]-Δy, y[end]+Δy))...)
    isnothing(clims) && (clims = extrema(data))
    return Makie.heatmap!(ax, x, y, transpose(data); colorrange=clims, _filter_kwargs(kwargs, Makie.Heatmap)...)
end

# 6. ReciprocalScatter with weights
@inline Makie.plot(reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}; kwargs...) = Makie.plot!(Makie.Figure(), reciprocalscatter, weights; kwargs...)
function Makie.plot!(fig::Makie.Figure, reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}; kwargs...)
    Makie.plot!(Makie.Axis(fig[1, 1]), reciprocalscatter, weights; kwargs...)
    return fig
end
function Makie.plot!(ax::Makie.AbstractAxis, reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}; fractional::Bool=true, weightmultiplier::Real=5.0, weightcolors=nothing, weightlabels=nothing, kwargs...)
    ax isa Makie.Axis || error("Makie.Axis required for reciprocal scatter")
    ax.aspect = get(kwargs, :aspect, Makie.DataAspect())
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xgridvisible = get(kwargs, :xgridvisible, true)
    ax.ygridvisible = get(kwargs, :ygridvisible, true)
    ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
    ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
    ax.xminorticksvisible = get(kwargs, :xminorticksvisible, true)
    ax.yminorticksvisible = get(kwargs, :yminorticksvisible, true)
    ax.xminorticks = get(kwargs, :xminorticks, Makie.IntervalsBetween(10))
    ax.yminorticks = get(kwargs, :yminorticks, Makie.IntervalsBetween(10))
    coordinates = fractional ? [tuple(c...) for c in fractionals(reciprocalscatter)] : [tuple(c...) for c in reciprocalscatter]
    length(first(coordinates)) == 2 || error("ReciprocalScatter with weights requires 2D coordinates")
    x, y = [c[1] for c in coordinates], [c[2] for c in coordinates]
    Makie.xlims!(ax, get(kwargs, :xlims, extrema(x))...)
    Makie.ylims!(ax, get(kwargs, :ylims, extrema(y))...)
    if !isnothing(weightlabels)
        for (i, label) in enumerate(weightlabels)
            Makie.scatter!(ax, [first(x)], [first(y)]; markersize=5, color=ith_color(weightcolors, i), label)
        end
        Makie.axislegend(ax)
    end
    num = size(weights, 2)
    return Makie.scatter!(
        ax, repeat(x, outer=num), repeat(y, outer=num);
        markersize=vec(weights) .* weightmultiplier .+ atol,
        color=reduce(vcat, [fill(ith_color(weightcolors, i), length(reciprocalscatter)) for i in 1:num]),
        strokewidth=0,
        _filter_kwargs(kwargs, Makie.Scatter)...
    )
end
@inline ith_color(weightcolors, i::Int) = weightcolors[i]
@inline ith_color(::Nothing, i::Int) = Makie.wong_colors()[mod1(i, 7)]

# 7. Line plot for ReciprocalPath with data
@inline Makie.plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}; kwargs...) = Makie.plot!(Makie.Figure(), path, data; kwargs...)
function Makie.plot!(fig::Makie.Figure, path::ReciprocalPath, data::AbstractMatrix{<:Number}; kwargs...)
    Makie.plot!(Makie.Axis(fig[1, 1]), path, data; kwargs...)
    return fig
end
function Makie.plot!(ax::Makie.AbstractAxis, path::ReciprocalPath, data::AbstractMatrix{<:Number}; kwargs...)
    ax isa Makie.Axis || error("Makie.Axis required for line plot")
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xlabel = get(kwargs, :xlabel, string(label(path)))
    ax.xticks = get(kwargs, :xticks, ticks(path))
    ax.xminorticksvisible = get(kwargs, :xminorticksvisible, true)
    ax.yminorticksvisible = get(kwargs, :yminorticksvisible, true)
    ax.xminorticks = get(kwargs, :xminorticks, Makie.IntervalsBetween(10))
    ax.yminorticks = get(kwargs, :yminorticks, Makie.IntervalsBetween(10))
    ax.xgridvisible = get(kwargs, :xgridvisible, true)
    ax.ygridvisible = get(kwargs, :ygridvisible, true)
    ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
    ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
    Makie.xlims!(ax, get(kwargs, :xlims, (0, distance(path)))...)
    return Makie.series!(ax, [distance(path, i) for i in eachindex(path)], transpose(data); linewidth=1, _filter_kwargs(kwargs, Makie.Series)...)
end

# 8. Scatter plot for ReciprocalPath with data and weights
@inline Makie.plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}; kwargs...) = Makie.plot!(Makie.Figure(), path, data, weights; kwargs...)
function Makie.plot!(fig::Makie.Figure, path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}; kwargs...)
    Makie.plot!(Makie.Axis(fig[1, 1]), path, data, weights; kwargs...)
    return fig
end
function Makie.plot!(ax::Makie.AbstractAxis, path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}; weightmultiplier::Real=5.0, weightcolors=nothing, weightlabels=nothing, kwargs...)
    ax isa Makie.Axis || error("Makie.Axis required for scatter plot")
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xlabel = get(kwargs, :xlabel, string(label(path)))
    ax.xticks = get(kwargs, :xticks, ticks(path))
    ax.xgridvisible = get(kwargs, :xgridvisible, true)
    ax.ygridvisible = get(kwargs, :ygridvisible, true)
    ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
    ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
    ax.xminorticksvisible = get(kwargs, :xminorticksvisible, true)
    ax.yminorticksvisible = get(kwargs, :yminorticksvisible, true)
    ax.xminorticks = get(kwargs, :xminorticks, Makie.IntervalsBetween(10))
    ax.yminorticks = get(kwargs, :yminorticks, Makie.IntervalsBetween(10))
    Makie.xlims!(ax, get(kwargs, :xlims, (0, distance(path)))...)
    x = repeat([distance(path, i) for i in eachindex(path)], outer=length(weights)÷length(path))
    y = repeat(vec(data), length(weights)÷length(data))
    if !isnothing(weightlabels)
        for (i, label) in enumerate(weightlabels)
            Makie.scatter!(ax, [first(x)], [first(y)]; markersize=5, color=ith_color(weightcolors, i), label=label)
        end
        Makie.axislegend(ax)
    end
    return Makie.scatter!(
        ax, x, y;
        markersize=vec(weights) .* weightmultiplier .+ atol,
        color=reduce(vcat, [fill(ith_color(weightcolors, i), length(data)) for i in 1:length(weights)÷length(data)]),
        alpha=0.8,
        strokewidth=0,
        _filter_kwargs(kwargs, Makie.Scatter)...
    )
end

# 9. Heatmap for ReciprocalPath with y data
@inline Makie.plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number}; kwargs...) = Makie.plot!(Makie.Figure(), path, y, data; kwargs...)
function Makie.plot!(fig::Makie.Figure, path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number}; kwargs...)
    Makie.Colorbar(fig[1, 2], Makie.plot!(Makie.Axis(fig[1, 1]), path, y, data; kwargs...))
    return fig
end
function Makie.plot!(ax::Makie.AbstractAxis, path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number}; clims=nothing, kwargs...)
    ax isa Makie.Axis || error("Makie.Axis required for heatmap")
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xlabel = get(kwargs, :xlabel, string(label(path)))
    ax.xticks = get(kwargs, :xticks, ticks(path))
    Makie.xlims!(ax, get(kwargs, :xlims, (0, distance(path)))...)
    Makie.ylims!(ax, get(kwargs, :ylims, extrema(y))...)
    isnothing(clims) && (clims = extrema(data))
    return Makie.heatmap!(ax, [distance(path, i) for i in eachindex(path)], y, transpose(data); colorrange=clims, _filter_kwargs(kwargs, Makie.Heatmap)...)
end

# 10. Multiple heatmaps for 3D data (multi-panel only)
@inline Makie.plot(reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractArray{<:Number, 3}; kwargs...) = Makie.plot!(Makie.Figure(), reciprocalspace, data; kwargs...)
function Makie.plot!(fig::Makie.Figure, reciprocalspace::Union{BrillouinZone, ReciprocalZone}, data::AbstractArray{<:Number, 3}; subtitles=nothing, titlesize::Int=8, nrow=nothing, ncol=nothing, clims=nothing, kwargs...)
    isnothing(clims) && (clims = extrema(data))
    isnothing(nrow) && (nrow = round(Int, sqrt(size(data, 3))))
    isnothing(ncol) && (ncol = ceil(Int, size(data, 3)/nrow))
    for i = 1:size(data, 3)
        Makie.plot!(Makie.Axis(fig[div(i-1, ncol)+1, (i-1)%ncol+1]), reciprocalspace, data[:, :, i]; title=isnothing(subtitles) ? nothing : subtitles[i], titlesize=titlesize, clims=clims, kwargs...)
    end
    Makie.Colorbar(fig[nrow+1, 1:ncol], limits=clims, vertical=false)
    return fig
end
@inline Makie.plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractArray{<:Number, 3}; kwargs...) = Makie.plot!(Makie.Figure(), path, y, data; kwargs...)
function Makie.plot!(fig::Makie.Figure, path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractArray{<:Number, 3}; subtitles=nothing, titlesize::Int=8, nrow=nothing, ncol=nothing, clims=nothing, kwargs...)
    isnothing(clims) && (clims = extrema(data))
    isnothing(nrow) && (nrow = round(Int, sqrt(size(data, 3))))
    isnothing(ncol) && (ncol = ceil(Int, size(data, 3)/nrow))
    for i = 1:size(data, 3)
        Makie.plot!(Makie.Axis(fig[div(i-1, ncol)+1, (i-1)%ncol+1]), path, y, data[:, :, i]; title=isnothing(subtitles) ? nothing : subtitles[i], titlesize=titlesize, clims=clims, kwargs...)
    end
    Makie.Colorbar(fig[nrow+1, 1:ncol], limits=clims, vertical=false)
    return fig
end

# 11. path plotting for (x, y)
@inline Makie.plot(x::AbstractVector{<:Number}, y::AbstractMatrix{<:Number}; kwargs...) = Makie.plot!(Makie.Figure(), x, y; kwargs...)
function Makie.plot!(fig::Makie.Figure, x::AbstractVector{<:Number}, y::AbstractMatrix{<:Number}; kwargs...)
    Makie.plot!(Makie.Axis(fig[1, 1]), x, y; kwargs...)
    return fig
end
function Makie.plot!(ax::Makie.AbstractAxis, x::AbstractVector{<:Number}, y::AbstractMatrix{<:Number}; kwargs...)
    ax isa Makie.Axis || error("Makie.Axis required for line plot")
    ax.titlesize = get(kwargs, :titlesize, 10)
    ax.xgridvisible = get(kwargs, :xgridvisible, true)
    ax.ygridvisible = get(kwargs, :ygridvisible, true)
    ax.xminorgridvisible = get(kwargs, :xminorgridvisible, true)
    ax.yminorgridvisible = get(kwargs, :yminorgridvisible, true)
    ax.xminorticksvisible = get(kwargs, :xminorticksvisible, true)
    ax.yminorticksvisible = get(kwargs, :yminorticksvisible, true)
    ax.xminorticks = get(kwargs, :xminorticks, Makie.IntervalsBetween(10))
    ax.yminorticks = get(kwargs, :yminorticks, Makie.IntervalsBetween(10))
    Makie.xlims!(ax, get(kwargs, :xlims, extrema(x))...)
    Makie.ylims!(ax, get(kwargs, :ylims, extrema(y))...)
    return Makie.series!(ax, x, transpose(y); linewidth=1, _filter_kwargs(kwargs, Makie.Series)...)
end

# 12. Assignment plotting
@inline Makie.plot(assignment::Assignment; kwargs...) = Makie.plot!(Makie.Figure(), assignment; kwargs...)
function Makie.plot!(fig::Makie.Figure, assignment::Assignment; kwargs...)
    plt = Makie.plot!(Makie.Axis(fig[1, 1]), assignment; kwargs...)
    isa(plt, Makie.Heatmap) && Makie.Colorbar(fig[1, 2], plt)
    return fig
end
function Makie.plot!(ax::Makie.AbstractAxis, assignment::Assignment; kwargs...)
    ax.title = get(kwargs, :title, str(assignment))
    return Makie.plot!(ax, Tuple(assignment.data)...; kwargs...)
end

end # module
