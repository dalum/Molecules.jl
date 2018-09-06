module Molecules

using LinearAlgebra

using Unitful: Unitful,
    @unit, @u_str,
    Length, Energy,
    eV, K
const kB = Unitful.k
@unit Å "Å" Angstrom 0.1u"nm" false

const localunits = Unitful.basefactors
function __init__()
    merge!(Unitful.basefactors, localunits)
end

using StaticArrays: SVector, SMatrix
using DataStructures: OrderedSet

export @Orbital_str, @orbital_str,
    Å, Atom, Bond, Molecule,
    makeaxes, slaterkoster

struct Orbital{N,L,M} end
macro Orbital_str(str::String)
    reg = r"^([0-9\*]+)([a-zA-Z]+)(?:_(?:\{(.+)\}|([a-zA-Z]+)))?$"
    m = match(reg, str)
    m === nothing && error("Invalid orbital")
    N = m[1] === nothing ? error("Invalid orbital") : m[1] == "*" ? nothing : parse(Int, m[1])
    L = m[2] === nothing ? error("Invalid orbital") : Symbol(m[2])
    M = m[3] === nothing ? m[4] === nothing ? Nothing : Symbol(m[4]) : Symbol(m[3])
    M === Symbol("*") && (M = nothing)

    N === nothing && M === nothing && return Orbital{<:Any,L,<:Any}
    N === nothing && return Orbital{<:Any,L,M}
    M === nothing && return Orbital{N,L,<:Any}
    return Orbital{N,L,M}
end

macro orbital_str(str::String)
    reg = r"^([0-9]+)([a-zA-Z]+)(?:_(?:\{(.+)\}|([a-zA-Z]+)))?$"
    m = match(reg, str)
    m === nothing && error("Invalid orbital")
    N = m[1] === nothing ? error("Invalid orbital") : parse(Int, m[1])
    L = m[2] === nothing ? error("Invalid orbital") : Symbol(m[2])
    M = m[3] === nothing ? m[4] === nothing ? Nothing : Symbol(m[4]) : Symbol(m[3])
    return Orbital{N,L,M}()
end

function Base.show(io::IO, ::Orbital{N,L,M}) where {N,L,M}
    if M === Nothing
        print(io, "orbital\"$N$L\"")
    elseif match(r"^[a-zA-Z]$", string(M)) !== nothing
        print(io, "orbital\"$N$(L)_$M\"")
    else
        print(io, "orbital\"$N$(L)_{$M}\"")
    end
end

struct Atom{SYMBOL,Z,N,T}
    position::SVector{3,Length{T}}
    axes::SMatrix{3,3,T}
    orbitals::OrderedSet{<:Orbital}

    function Atom{SYMBOL,Z,N}(x::Length{T}, y::Length{T}, z::Length{T}, axes::SMatrix{3,3,T,9}, orbitals::Orbital...) where {SYMBOL,Z,N,T}
        return new{SYMBOL,Z,N,T}(SVector{3,Length{T}}(x, y, z), axes, OrderedSet(orbitals))
    end
end # module

struct Bond{A1<:Atom,A2<:Atom}
    atoms::Tuple{A1,A2}
end
Bond(atom1, atom2) = Bond((atom1, atom2))

struct Molecule{T}
    position::SVector{3,Length{T}}
    axes::SMatrix{3,3,T}
    atoms::Vector{<:Atom}
    bonds::Vector{<:Bond}
end

##################################################
# Elements
##################################################

const Hydrogen = Atom{:H, 1, 0}
const Deuterium = Atom{:D, 1, 1}
const Tritium = Atom{:T, 1, 2}
const Carbon = Atom{:C, 6, 6}
const Carbon13 = Atom{:C, 6, 7}
const Nitrogen = Atom{:N, 7, 7}
const Oxygen = Atom{:O, 8, 8}

##################################################
# Coordinates
##################################################

function makeaxes(θ, ϕ)
    R1 = LinearAlgebra.Givens(3, 1, cos(ϕ), sin(ϕ))
    R2 = LinearAlgebra.Givens(2, 1, cos(θ), sin(θ))
    return SMatrix{3, 3}(R2*R1*Matrix(1.0I, 3, 3))
end

##################################################
# Hamiltonians
##################################################

"""
    hamiltonian(x)

Construct the matrix corresponding to the hamiltonian operator of `x`,
where `x` is an atom, bond or molecule.  For bonds, this will compute
the overlap integrals using the Slater-Koster scheme, for atoms, it
will return a diagonal matrix with the onsite energies of the
orbitals.

"""
function hamiltonian(atom::Atom{<:Any,<:Any,<:Any,T}) where {T}
    return [s == s′ ? onsiteenergy(atom, s) : zero(T)*eV for s in atom.orbitals, s′ in atom.orbitals]::Matrix{<:Energy{T}}
end

function hamiltonian(bond::Bond{<:Atom{<:Any,<:Any,<:Any,T},<:Atom{<:Any,<:Any,<:Any,T}}) where {T}
    atom1, atom2 = bond.atoms
    δr = atom2.position - atom1.position
    f = (val, δr, s1, s2) -> bondenergy(val, δr, atom1, atom2, s1, s2)
    return [slaterkoster(f, δr, s1, s2) for s1 in atom1.orbitals, s2 in atom2.orbitals]::Matrix{<:Energy{T}}
end


## Slater-Koster scheme

slaterkoster(f, δr::SVector{3}, s1::Orbital"*s", s2::Orbital"*s")::Energy = f(Val(:σ), norm(δr), s1, s2)

for (i, x) in [(1, :(:x)), (2, :(:y)), (3, :(:z))]
    @eval slaterkoster(f, δr::SVector{3}, s1::Orbital"*s", s2::Orbital{<:Any, :p, $x})::Energy = normalize(δr)[$i] * f(Val(:σ), norm(δr), s1, s2)

    @eval slaterkoster(f, δr::SVector{3}, s1::Orbital{<:Any,:p,$x}, s2::Orbital{<:Any,:s,Nothing})::Energy = normalize(δr)[$i] * f(Val(:σ), norm(δr), s1, s2)

    for (j, y) in [(1, :(:x)), (2, :(:y)), (3, :(:z))]
        if i == j
            @eval function slaterkoster(f, δr::SVector{3}, s1::Orbital{<:Any,:p,$x}, s2::Orbital{<:Any,:p,$x})::Energy
                d = normalize(δr)
                return d[$i]^2 * f(Val(:σ), norm(δr), s1, s2) + (1 - d[$i]^2) * f(Val(:π), norm(δr), s1, s2)
            end
        else
            @eval function slaterkoster(f, δr::SVector{3}, s1::Orbital{<:Any,:p,$x}, s2::Orbital{<:Any,:p,$y})::Energy
                d = normalize(δr)
                d[$i] * d[$j] * (f(Val(:σ), norm(δr), s1, s2) - f(Val(:π), norm(δr), s1, s2))
            end
        end
    end
end

include("energies.jl")

end
