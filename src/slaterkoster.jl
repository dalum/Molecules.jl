"""
    slaterkoster(; kwargs...)

Create a function to be passed to `hamiltonian` or `overlap` for
calculating the matrix elements using the Slater-Koster scheme.

"""

slaterkoster(; kwargs...) = (args...) -> slaterkoster(args...; kwargs...)

function slaterkoster(r, a1; kwargs...)
    return kwargs[Symbol(first(string(last(a1))))]
end

function slaterkoster(r, a1, a2; kwargs...)
    o1, o2 = last(a1), last(a2)
    d = r ./ norm(r)
    s(x) = kwargs[Symbol(first(a1), first(a2), "_", first(string(last(a1))), first(string(last(a2))), x)]

    if (o1, o2) == (:s, :s)
        return s(:σ)
    end

    for (i, x) in [(1, :x), (2, :y), (3, :z)]
        if (o1, o2) == (:s, Symbol(:p, x))
            return d[i]*s(:σ)
        end
        if (o1, o2) == (Symbol(:p, x), :s)
            return d[i]*s(:σ)
        end

        for (j, y) in [(1, :x), (2, :y), (3, :z)]
            if (o1, o2) == (Symbol(:p, x), Symbol(:p, x))
                return d[i]^2*s(:σ) + (1 - d[i]^2)*s(:π)
            end
            if (o1, o2) == (Symbol(:p, x), Symbol(:p, y))
                return d[i]*d[j]*(s(:σ) - s(:π))
            end
        end
    end

    error("orbitals not supported: $a1, $a2")
end
