##################################################
# Visualization
##################################################

using RecipesBase
using Colors: RGB, HSV
using Statistics

elementcolor(::Hydrogen) = RGB(1., 1., 1.)
elementcolor(::Carbon) = RGB(0., 0., 0.)
elementcolor(::Nitrogen) = RGB(0., 0., 1.)
elementcolor(::Oxygen) = RGB(1., 0., 0.)
volume(::Hydrogen) = 5
volume(::Carbon) = 7
volume(::Nitrogen) = 8
volume(::Oxygen) = 9

@recipe function f(mol::Molecule; unitradius=15, local_axes=false)
    xs = [position(atom)[1] for atom in mol.atoms]
    ys = [position(atom)[2] for atom in mol.atoms]
    zs = [position(atom)[3] for atom in mol.atoms]
    max_distance = max(maximum(xs) - minimum(xs),
                       maximum(ys) - minimum(ys),
                       maximum(zs) - minimum(zs))
    max_distance += sqrt(max_distance)

    markercolor --> elementcolor.(mol.atoms)
    markersize --> unitradius .* ((3/4π) .* volume.(mol.atoms)) .^ 1/3
    @series begin
        segments = [Float64[] for _ in 1:3]
        linewidth --> 2
        linecolor --> :black
        seriestype := :path3d
        primary := false
        markershape := :none

        for bond in mol.bonds
            p_i = position(bond.atom1)
            p_j = position(bond.atom2)
            push!.(segments, p_i, p_j, NaN)
        end
        (segments[1], segments[2], segments[3])
    end

    if local_axes
        for (j, col) in zip(1:3, (:green, :blue, :red))
            @series begin
                segments = [Float64[] for _ in 1:3]
                linewidth --> 5
                linecolor --> col
                seriestype := :path3d
                primary := false
                markershape := :none

                for atom in mol.atoms
                    # for (p0, U) in zip(axes.origins, axes.axes)
                    p_i = position(atom)
                    p_j = p_i + atom.axes[:, j]
                    push!.(segments, p_i, p_j, NaN)
                end
                (segments[1], segments[2], segments[3])
            end
        end
    end

    xlims = mean(xs) - max_distance/2, mean(xs) + max_distance/2
    ylims = mean(ys) - max_distance/2, mean(ys) + max_distance/2
    zlims = mean(zs) - max_distance/2, mean(zs) + max_distance/2
    xticks = floor(xlims[1]):ceil(xlims[2])
    yticks = floor(ylims[1]):ceil(ylims[2])
    zticks = floor(zlims[1]):ceil(zlims[2])

    xlims --> xlims
    xticks --> (xticks, string.(xticks))
    ylims --> ylims
    yticks --> (yticks, string.(yticks))
    zlims --> zlims
    zticks --> (zticks, string.(zticks))

    (xs, ys, zs)
end

# Projection

@recipe function f(mol::Molecule, weights::Vector; unitradius=15, local_axes=false, pscale = 10)
    xs = [position(atom)[1] for atom in mol.atoms]
    ys = [position(atom)[2] for atom in mol.atoms]
    zs = [position(atom)[3] for atom in mol.atoms]
    max_distance = max(maximum(xs) - minimum(xs),
                       maximum(ys) - minimum(ys),
                       maximum(zs) - minimum(zs))
    max_distance += sqrt(max_distance)

    xlims = mean(xs) - max_distance/2, mean(xs) + max_distance/2
    ylims = mean(ys) - max_distance/2, mean(ys) + max_distance/2
    zlims = mean(zs) - max_distance/2, mean(zs) + max_distance/2
    xticks = floor(xlims[1]):ceil(xlims[2])
    yticks = floor(ylims[1]):ceil(ylims[2])
    zticks = floor(zlims[1]):ceil(zlims[2])

    xlims --> xlims
    xticks --> (xticks, string.(xticks))
    ylims --> ylims
    yticks --> (yticks, string.(yticks))
    zlims --> zlims
    zticks --> (zticks, string.(zticks))

    @series begin
        xs = [position(atom)[1] for atom in mol.atoms]
        ys = [position(atom)[2] for atom in mol.atoms]
        zs = [position(atom)[3] for atom in mol.atoms]
        markercolor --> elementcolor.(mol.atoms)
        markersize --> unitradius .* ((3/4π) .* volume.(mol.atoms)) .^ 1/3
        seriesalpha --> 0.1

        (xs, ys, zs)
    end

    @series begin
        segments = [Float64[] for _ in 1:3]
        linewidth --> 2
        linecolor --> :black
        seriestype := :path3d
        primary := false
        markershape := :none

        for bond in mol.bonds
            p_i = position(bond.atom1)
            p_j = position(bond.atom2)
            push!.(segments, p_i, p_j, NaN)
        end
        (segments[1], segments[2], segments[3])
    end

    if local_axes
        for (j, col) in zip(1:3, (:green, :blue, :red))
            @series begin
                segments = [Float64[] for _ in 1:3]
                linewidth --> 5
                linecolor --> col
                seriestype := :path3d
                primary := false
                markershape := :none

                for atom in mol.atoms
                    # for (p0, U) in zip(axes.origins, axes.axes)
                    p_i = position(atom)
                    p_j = p_i + atom.axes[:, j]
                    push!.(segments, p_i, p_j, NaN)
                end
                (segments[1], segments[2], segments[3])
            end
        end
    end

    @series begin
        xs, ys, zs = [Float64[] for _ in 1:3]
        markercolor := markercolor = []#HSV.(180/π .* angle.(weights), 1, 1)
        markersize := markersize = []#unitradius .* abs.(weights)

        i = 1
        for atom in mol.atoms
            x, y, z = position(atom)
            axes = unitradius .* atom.axes ./ pscale
            for orbital in atom.orbitals
                if orbital isa Orbital"*s"
                    push!(xs, x)
                    push!(ys, y)
                    push!(zs, z)
                    push!(markersize, unitradius*abs(weights[i]))
                    push!(markercolor, HSV(180/π*angle(weights[i]), 1, 1))
                end
                for (idx, O) in ((1, Orbital"*p_x"), (2, Orbital"*p_y"), (3, Orbital"*p_z"))
                    if orbital isa O
                        a = axes .* abs(weights[i])
                        push!(xs, x - a[1, idx], x + a[1, idx])
                        push!(ys, y - a[2, idx], y + a[2, idx])
                        push!(zs, z - a[3, idx], z + a[3, idx])
                        for _ in 1:2
                            push!(markersize, unitradius*abs(weights[i]))
                            push!(markercolor, HSV(180/π*angle(weights[i]), 1, 1))
                        end
                    end
                end
                i += 1
            end
        end

        (xs, ys, zs)
    end
end

# With arrows

@recipe function f(g::Function, mol::Molecule, arrows::Vector; unitradius=15)
    xs = [position(atom)[1] for atom in mol.atoms]
    ys = [position(atom)[2] for atom in mol.atoms]
    zs = [position(atom)[3] for atom in mol.atoms]
    max_distance = max(maximum(xs) - minimum(xs),
                       maximum(ys) - minimum(ys),
                       maximum(zs) - minimum(zs))
    max_distance += sqrt(max_distance)

    markercolor --> elementcolor.(mol.atoms)
    markersize --> unitradius .* ((3/4π) .* volume.(mol.atoms)) .^ 1/3
    @series begin
        segments = [Float64[] for _ in 1:3]
        linewidth --> 2
        linecolor --> :black
        seriestype := :path3d
        primary := false
        markershape := :none

        for bond in mol.bonds
            p_i = position(bond.atom1)
            p_j = position(bond.atom2)
            push!.(segments, p_i, p_j, NaN)
        end
        (segments[1], segments[2], segments[3])
    end

    for (j, col) in zip(1:3, (:green, :blue, :red))
        @series begin
            segments = [Float64[] for _ in 1:3]
            linewidth --> 5
            linecolor --> col
            seriestype := :path3d
            primary := false
            markershape := :none

            for (atom, arrow) in zip(filter(g, mol.atoms), arrows)
                # for (p0, U) in zip(axes.origins, axes.axes)
                p_i = position(atom)
                p_j = p_i + arrow[j]
                push!.(segments, p_i, p_j, NaN)
            end
            (segments[1], segments[2], segments[3])
        end
    end

    xlims = mean(xs) - max_distance/2, mean(xs) + max_distance/2
    ylims = mean(ys) - max_distance/2, mean(ys) + max_distance/2
    zlims = mean(zs) - max_distance/2, mean(zs) + max_distance/2
    xticks = floor(xlims[1]):ceil(xlims[2])
    yticks = floor(ylims[1]):ceil(ylims[2])
    zticks = floor(zlims[1]):ceil(zlims[2])

    xlims --> xlims
    xticks --> (xticks, string.(xticks))
    ylims --> ylims
    yticks --> (yticks, string.(yticks))
    zlims --> zlims
    zticks --> (zticks, string.(zticks))

    (xs, ys, zs)
end
