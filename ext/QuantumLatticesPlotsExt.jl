module QuantumLatticesPlotsExt

using RecipesBase: @layout, @recipe, @series
using QuantumLattices: AbstractLattice, Assignment, Bond, BrillouinZone, Data, FractionalReciprocalSpace, Neighbors, Point, ReciprocalCurve, ReciprocalPath, ReciprocalScatter, ReciprocalZone
using QuantumLattices: atol, bonds, dimension, distance, fractionals, getcontent, label, rank, scalartype, str, ticks
using QuantumLattices.Frameworks: seriestype

"""
    @recipe plot(lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; siteon=false)

Define recipe for visualization of a lattice.
"""
@recipe function plot(lattice::AbstractLattice, neighbors::Union{Int, Neighbors}, filter::Function=bond->true; siteon=false)
    title --> String(getcontent(lattice, :name))
    titlefontsize --> 10
    legend := false
    aspect_ratio := :equal
    framestyle --> :box
    if isa(neighbors, Int) || 0∈keys(neighbors)
        @series begin
            seriestype := :scatter
            markerstrokewidth --> 0
            coordinates = NTuple{dimension(lattice), scalartype(lattice)}[]
            sites = String[]
            for i in eachindex(lattice)
                bond = Bond(Point(i, lattice[i], zero(lattice[i])))
                if filter(bond)
                    push!(coordinates, Tuple(lattice[i]))
                    push!(sites, siteon ? string(i) : "")
                end
            end
            series_annotations := sites
            coordinates
        end
    end
    data = [Float64[] for _ in 1:dimension(lattice)]
    for bond in bonds(lattice, neighbors)
        if length(bond)>0 && filter(bond)
            for i = 1:dimension(lattice)
                for point in bond
                    push!(data[i], point.rcoordinate[i])
                end
                push!(data[i], NaN)
            end
        end
    end
    Tuple(data)
end

"""
    @recipe plot(reciprocalspace::FractionalReciprocalSpace; fractional=false, autolims=true)

Define recipe for visualization of a reciprocal space with fractional coordinates.

When `fractional` is `true`, fractional coordinates will be plotted. Otherwise Cartesian coordinates will be plotted.
"""
@recipe function plot(reciprocalspace::FractionalReciprocalSpace; fractional=false, autolims=true)
    seriestype --> :scatter
    markerstrokewidth --> 0
    titlefontsize --> 10
    aspect_ratio --> :equal
    legend --> false
    minorgrid --> true
    framestyle --> :box
    if fractional
        N = rank(reciprocalspace)
        coordinates = map(NTuple{N, scalartype(reciprocalspace)}, fractionals(reciprocalspace))
        N>0 && (xlabel-->label(reciprocalspace, 1); autolims || (xlims-->extrema(v->v[1], coordinates)))
        N>1 && (ylabel-->label(reciprocalspace, 2); autolims || (ylims-->extrema(v->v[2], coordinates)))
        N>2 && (zlabel-->label(reciprocalspace, 3); autolims || (zlims-->extrema(v->v[3], coordinates)))
        coordinates
    else
        dim = dimension(reciprocalspace)
        coordinates = map(NTuple{dim, scalartype(reciprocalspace)}, reciprocalspace)
        dim>0 && (xlabel-->label(reciprocalspace, 1); autolims || (xlims-->extrema(v->v[1], coordinates)))
        dim>1 && (ylabel-->label(reciprocalspace, 2); autolims || (ylims-->extrema(v->v[2], coordinates)))
        dim>2 && (zlabel-->label(reciprocalspace, 3); autolims || (zlims-->extrema(v->v[3], coordinates)))
        coordinates
    end
end

"""
    @recipe plot(path::ReciprocalPath)

Define the recipe for the visualization of a reciprocal path.
"""
@recipe function plot(path::ReciprocalPath)
    seriestype --> :scatter
    markerstrokewidth --> 0
    titlefontsize --> 10
    legend --> false
    aspect_ratio --> :equal
    minorgrid --> true
    framestyle --> :box
    dim = dimension(path)
    dim>0 && (xlabel --> label(path, 1))
    dim>1 && (ylabel --> label(path, 2))
    dim>2 && (zlabel --> label(path, 3))
    x, y, text = [], [], []
    for i = 1:length(path.contents)
        if length(path.contents[i])>0
            push!(x, path.contents[i][1][1])
            push!(y, path.contents[i][1][2])
            push!(text, string(path.labels[i].first))
            j = i%length(path.labels)+1
            if length(path.contents[i])>1 && path.labels[i].second≠path.labels[j].first
                push!(x, path.contents[i][end][1])
                push!(y, path.contents[i][end][2])
                push!(text, string(path.labels[i].second))
            end
        end
    end
    annotation := (x, y, text)
    coordinates = map(Tuple, path)
end

"""
    @recipe plot(curve::ReciprocalCurve)

Define the recipe for the visualization of a reciprocal curve.
"""
@recipe function plot(curve::ReciprocalCurve)
    seriestype --> :scatter
    markerstrokewidth --> 0
    titlefontsize --> 10
    legend --> false
    aspect_ratio --> :equal
    minorgrid --> true
    framestyle --> :box
    dim = dimension(curve)
    dim>0 && (xlabel --> label(curve, 1))
    dim>1 && (ylabel --> label(curve, 2))
    dim>2 && (zlabel --> label(curve, 3))
    map(Tuple, curve)
end

# plot utilities
"""
    @recipe plot(reciprocalspace::BrillouinZone, data::AbstractMatrix{<:Number})
    @recipe plot(reciprocalspace::ReciprocalZone, data::AbstractMatrix{<:Number})

Define the recipe for the heatmap visualization of data on a Brillouin/reciprocal zone.
"""
const heatmap = quote
    @assert length(reciprocalspace.reciprocals)==2 "plot error: only two dimensional reciprocal spaces are supported."
    seriestype --> :heatmap
    titlefontsize --> 10
    aspect_ratio --> :equal
    x, y = range(reciprocalspace, 1), range(reciprocalspace, 2)
    Δx, Δy= x[2]-x[1], y[2]-y[1]
    xlims --> (x[1]-Δx, x[end]+Δx)
    ylims --> (y[1]-Δy, y[end]+Δy)
    clims --> extrema(data)
    xlabel --> label(reciprocalspace, 1)
    ylabel --> label(reciprocalspace, 2)
    x, y, data
end
@eval @recipe plot(reciprocalspace::BrillouinZone, data::AbstractMatrix{<:Number}) = $heatmap
@eval @recipe plot(reciprocalspace::ReciprocalZone, data::AbstractMatrix{<:Number}) = $heatmap

"""
    @recipe plot(reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}; fractional=true, weightmultiplier=1.0, weightcolors=nothing, weightlabels=nothing)

Define the recipe for the scatter visualization of reciprocal points with a series of weights.
"""
@recipe function plot(reciprocalscatter::ReciprocalScatter, weights::AbstractMatrix{<:Number}; fractional=true, weightmultiplier=1.0, weightcolors=nothing, weightlabels=nothing)
    seriestype --> :scatter
    markerstrokewidth --> 0
    legend --> !isnothing(weightlabels)
    point = fractional ? [ntuple(i->0, rank(reciprocalscatter))] : [ntuple(i->0, dimension(reciprocalscatter))]
    fractional := fractional
    autolims := false
    for (i, index) in enumerate(axes(weights, 2))
        color = isnothing(weightcolors) ? i : weightcolors[i]
        @series begin
            markersize := atol
            markercolor := color
            label := isnothing(weightlabels) ? "" : weightlabels[i]
            point
        end
        @series begin
            markersize := weights[:, index]*weightmultiplier .+ atol
            markercolor := color
            label := false
            reciprocalscatter
        end
    end
end

"""
    @recipe plot(path::ReciprocalPath, data::AbstractMatrix{<:Number})

Define the recipe for the line visualization of data along a reciprocal path.
"""
const line = quote
    seriestype --> :path
    titlefontsize --> 10
    legend --> false
    xlims --> (0, distance(path))
    minorgrid --> true
    xminorticks --> 10
    yminorticks --> 10
    xticks --> ticks(path)
    xlabel --> string(label(path))
    [distance(path, i) for i=1:length(path)], data
end
@eval @recipe plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}) = $line

"""
    @recipe plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}; weightmultiplier=5.0, weightcolors=nothing, weightlabels=nothing)

Define the recipe for the scatter visualization of data along a reciprocal path with a series of weights.
"""
const scatter = quote
    seriestype --> :scatter
    markerstrokewidth := 0
    legend --> !isnothing(weightlabels)
    xlims --> (0, distance(path))
    alpha --> 0.8
    for i in 1:size(weights, 3)
        @series begin
            markersize := atol
            markercolor := isnothing(weightcolors) ? i : weightcolors[i]
            label := isnothing(weightlabels) ? i : weightlabels[i]
            [(0, 0)]
        end
    end
    for (i, index) in enumerate(axes(weights, 3))
        @series begin
            markersize := weights[:, :, index]*weightmultiplier .+ atol
            markercolor := isnothing(weightcolors) ? i : weightcolors[i]
            label := false
            path, data
        end
    end
end
@eval @recipe plot(path::ReciprocalPath, data::AbstractMatrix{<:Number}, weights::AbstractArray{<:Number, 3}; weightmultiplier=5.0, weightcolors=nothing, weightlabels=nothing) = $scatter

"""
    @recipe plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number})

Define the recipe for the heatmap visualization of data on x-y plain with x axis being a reciprocal path.
"""
@recipe function plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractMatrix{<:Number})
    seriestype --> :heatmap
    titlefontsize --> 10
    xticks --> ticks(path)
    xlabel --> string(label(path))
    xlims --> (0, distance(path))
    ylims --> (minimum(y), maximum(y))
    clims --> extrema(data)
    [distance(path, i) for i=1:length(path)], y, data
end

"""
    @recipe plot(reciprocalspace::BrillouinZone, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing)
    @recipe plot(reciprocalspace::ReciprocalZone, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing)
    @recipe plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing)

Define the recipe for the heatmap visualization of a series of data on
1) a Brillouin zone,
2) a reciprocal zone,
3) x-y plain with x axis being a reciprocal path.
"""
setup(expr::Expr) = quote
    seriestype := :heatmap
    isnothing(clims) && (clims = extrema(data))
    isnothing(nrow) && (nrow = round(Int, sqrt(size(data)[3])))
    isnothing(ncol) && (ncol = ceil(Int, size(data)[3]/nrow))
    layout := @layout [(nrow, ncol); b{0.05h}]
    colorbar := false
    for i = 1:size(data)[3]
        @series begin
            isnothing(subtitles) || begin
                title := subtitles[i]
                titlefontsize := subtitlefontsize
            end
            subplot := i
            clims := clims
            $expr
        end
    end
    subplot := nrow*ncol+1
    xlims := (minimum(clims), maximum(clims))
    xlabel := ""
    ylims := (0, 1)
    yticks := (0:1, ("", ""))
    ylabel := ""
    LinRange(clims..., 100), [0, 1], [LinRange(clims..., 100)'; LinRange(clims..., 100)']
end
@eval @recipe plot(reciprocalspace::BrillouinZone, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing) = $(setup(:(reciprocalspace, data[:, :, i])))
@eval @recipe plot(reciprocalspace::ReciprocalZone, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing) = $(setup(:(reciprocalspace, data[:, :, i])))
@eval @recipe plot(path::ReciprocalPath, y::AbstractVector{<:Number}, data::AbstractArray{<:Number, 3}; subtitles=nothing, subtitlefontsize=8, nrow=nothing, ncol=nothing, clims=nothing) = $(setup(:(path, y, data[:, :, i])))

"""
    @recipe plot(assignment::Assignment)

Define the recipe for the visualization of an assignment of an algorithm.
"""
@recipe function plot(assignment::Assignment)
    title --> str(assignment)
    titlefontsize --> 10
    attr = seriestype(assignment.data)
    isnothing(attr) || begin
        seriestype --> attr
        attr==:path && begin
            legend --> false
            minorgrid --> true
            xminorticks --> 10
            yminorticks --> 10
        end
    end
    Tuple(assignment.data)
end

end # module
