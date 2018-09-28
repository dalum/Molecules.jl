##################################################
# Visualization
##################################################

using Unitful: uconvert, NoUnits, unit, ustrip
using RecipesBase
using Colors: RGB
using Statistics
using Molecules; using Molecules: Hydrogen, Carbon, Nitrogen, Oxygen, Å

elementcolor(::Hydrogen) = RGB(1., 1., 1.)
elementcolor(::Carbon) = RGB(0., 0., 0.)
elementcolor(::Nitrogen) = RGB(0., 0., 1.)
elementcolor(::Oxygen) = RGB(1., 0., 0.)
radius(::Hydrogen) = 10
radius(::Carbon) = 15
radius(::Nitrogen) = 18
radius(::Oxygen) = 20

@recipe function f(mol::Molecule; relativesize=10)
    a = 1Å # characteristic length scale
    xs = [atom.position[1]  for atom in mol.atoms]
    ys = [atom.position[2]  for atom in mol.atoms]
    zs = [atom.position[3]  for atom in mol.atoms]
    max_distance = a + max(maximum(xs) - minimum(xs),
                           maximum(ys) - minimum(ys),
                           maximum(zs) - minimum(zs))

    markercolor --> elementcolor.(mol.atoms)
    markersize --> relativesize# * radius.(mol.atoms)
    @series begin
        segments = [Float64[] for _ in 1:3]
        linewidth --> 2
        linecolor --> :black
        seriestype := :path3d
        primary := false
        markershape := :none

        for bond in mol.bonds
            p_i = float.(uconvert.(NoUnits, bond.atom1.position ./ a))
            p_j = float.(uconvert.(NoUnits, bond.atom2.position ./ a))
            push!.(segments, p_i, p_j, NaN)
        end
        (segments[1], segments[2], segments[3])
    end

    xlims = uconvert.(NoUnits, (mean(xs) - max_distance/2, mean(xs) + max_distance/2) ./ a)
    ylims = uconvert.(NoUnits, (mean(ys) - max_distance/2, mean(ys) + max_distance/2) ./ a)
    zlims = uconvert.(NoUnits, (mean(zs) - max_distance/2, mean(zs) + max_distance/2) ./ a)
    xticks = floor(xlims[1]):(a.val):ceil(xlims[2])
    yticks = floor(ylims[1]):(a.val):ceil(ylims[2])
    zticks = floor(zlims[1]):(a.val):ceil(zlims[2])

    xlims --> xlims
    xticks --> (xticks, string.(xticks .* unit(a)))
    ylims --> ylims
    yticks --> (yticks, string.(yticks .* unit(a)))
    zlims --> zlims
    zticks --> (zticks, string.(zticks .* unit(a)))

    (ustrip.(xs), ustrip.(ys), ustrip.(zs))
end

# @recipe function f(axes::Axes)
#     a = 1Å # characteristic length scale
#     A = hcat(axes.origins...)
#     xs, ys, zs = A[1,:], A[2,:], A[3,:]
#     max_distance = a + max(maximum(xs) - minimum(xs),
#                            maximum(ys) - minimum(ys),
#                            maximum(zs) - minimum(zs))

#     xlims = uconvert.(NoUnits, (mean(xs) - max_distance/2, mean(xs) + max_distance/2) ./ a)
#     ylims = uconvert.(NoUnits, (mean(ys) - max_distance/2, mean(ys) + max_distance/2) ./ a)
#     zlims = uconvert.(NoUnits, (mean(zs) - max_distance/2, mean(zs) + max_distance/2) ./ a)
#     xticks = floor(xlims[1]):(a.val):ceil(xlims[2])
#     yticks = floor(ylims[1]):(a.val):ceil(ylims[2])
#     zticks = floor(zlims[1]):(a.val):ceil(zlims[2])

#     xlims --> xlims
#     xticks --> (xticks, string.(xticks .* unit(a)))
#     ylims --> ylims
#     yticks --> (yticks, string.(yticks .* unit(a)))
#     zlims --> zlims
#     zticks --> (zticks, string.(zticks .* unit(a)))

#     for (j, col) in zip(1:3, (:green, :blue, :red))
#         @series begin
#             segments = [Segments() for _ in 1:3]
#             linewidth --> 5
#             linecolor --> col
#             seriestype := :path3d
#             primary := false
#             markershape := :none

#             for (p0, U) in zip(axes.origins, axes.axes)
#                 p_i = float.(uconvert.(NoUnits, p0 / a))
#                 p_j = p_i + U[:, j]
#                 push!.(segments, p_i, p_j)
#             end
#             (segments[1].pts, segments[2].pts, segments[3].pts)
#         end
#     end
# end
