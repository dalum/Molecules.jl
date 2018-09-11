module Molecules

using LinearAlgebra

using Unitful: Unitful,
    @unit, @u_str,
    Energy, Length, Temperature,
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
    hamiltonian, makeaxes, slaterkoster

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

mutable struct Atom{T,SYMBOL,Z,N}
    position::SVector{3,<:Length{T}}
    axes::SMatrix{3,3,T}
    orbitals::OrderedSet{<:Orbital}

    Atom{T,SYMBOL,Z,N}(axes::SMatrix{3,3,T,9}, x::Length{T}, y::Length{T}, z::Length{T}, orbitals::Orbital...) where {T,SYMBOL,Z,N} = new{T,SYMBOL,Z,N}(
        SVector{3}(x, y, z),
        axes,
        OrderedSet(orbitals))
end # module

Atom{<:Any,SYMBOL,Z,N}(x::Length{T}, y::Length{T}, z::Length{T}, orbitals::Orbital...) where {T,SYMBOL,Z,N} = Atom{T,SYMBOL,Z,N}(
    makeaxes(T, zero(T), zero(T)),
    x,
    y,
    z,
    orbitals...)

Atom{<:Any,SYMBOL,Z,N}(position::Vector{<:Length{T}}, orbitals::Orbital...) where {T,SYMBOL,Z,N} = Atom{T,SYMBOL,Z,N}(
    makeaxes(T, zero(T), zero(T)),
    position...,
    orbitals...)

struct Bond{A1<:Atom,A2<:Atom}
    atom1::A1
    atom2::A2

    function Bond{A1,A2}(atom1::A1, atom2::A2) where {A1<:Atom,A2<:Atom}
        if hash(atom1) > hash(atom2)
            return new{A1,A2}(atom1, atom2)
        else
            return new{A2,A1}(atom2, atom1)
        end
    end
end
Bond(atom1::A1, atom2::A2) where {A1, A2} = Bond{A1,A2}(atom1, atom2)

mutable struct Molecule{T}
    position::SVector{3,<:Length{T}}
    axes::SMatrix{3,3,T}
    atoms::OrderedSet{<:Atom{T}}
    bonds::Set{<:Bond}
end

Molecule(atoms::OrderedSet{<:Atom{T}}, bonds) where T = Molecule(
    SVector{3}(zero(T)*Å, zero(T)*Å, zero(T)*Å),
    makeaxes(T, zero(T), zero(T)),
    atoms,
    Set(bonds))

Molecule(atoms::AbstractVector{<:Atom{T}}, bonds) where T = Molecule(
    SVector{3}(zero(T)*Å, zero(T)*Å, zero(T)*Å),
    makeaxes(T, zero(T), zero(T)),
    OrderedSet(atoms),
    Set(bonds))

##################################################
# Chemical properties
##################################################

protonnumber(::Atom{T,SYMBOL,Z,N}) where {T,SYMBOL,Z,N} = Z
neutronnumber(::Atom{T,SYMBOL,Z,N}) where {T,SYMBOL,Z,N} = N
symbol(::Atom{T,SYMBOL,Z,N}) where {T,SYMBOL,Z,N} = SYMBOL

countorbitals(atom::Atom) = length(atom.orbitals)
countorbitals(mol::Molecule) = sum(map(countorbitals, mol.atoms))

function countelectrons(atom::Atom)
    n = protonnumber(atom)
    n <= 2 && return n
    n <= 10 && return n - 2 # Not very generic ...
end
countelectrons(mol::Molecule) = sum(countelectrons(atom) for atom in mol.atoms)

function chemicalpotential(mol::Molecule)
    H = hamiltonian(mol)
    n = Int(ceil(countelectrons(mol) / 2))
    energies = eigvals(H / 1eV) * 1eV
    return (energies[n] + energies[n + 1]) / 2, (energies[n + 1] - energies[n]) / 2
end

##################################################
# Elements
##################################################

const Hydrogen = Atom{<:Any, :H, 1, 0}
const Deuterium = Atom{<:Any, :D, 1, 1}
const Tritium = Atom{<:Any, :T, 1, 2}
const Carbon = Atom{<:Any, :C, 6, 6}
const Carbon13 = Atom{<:Any, :C, 6, 7}
const Nitrogen = Atom{<:Any, :N, 7, 7}
const Oxygen = Atom{<:Any, :O, 8, 8}

##################################################
# Coordinates
##################################################

function makeaxes(::Type{T}, θ::T, ϕ::T) where T
    R1 = LinearAlgebra.Givens(3, 1, cos(ϕ), sin(ϕ))
    R2 = LinearAlgebra.Givens(2, 1, cos(θ), sin(θ))
    return SMatrix{3, 3}(R2*R1*Matrix(one(T)*I, 3, 3))
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
function hamiltonian(atom::Atom)
    return [s == s′ ? onsiteenergy(atom, s) : 0.0eV for s in atom.orbitals, s′ in atom.orbitals]
end

function hamiltonian(atom1::Atom, atom2::Atom)
    δr = atom1.position - atom2.position
    f = (val, δr, s1, s2) -> bondenergy(val, δr, atom1, atom2, s1, s2)
    return [slaterkoster(f, δr, s1, s2) for s1 in atom1.orbitals, s2 in atom2.orbitals]
end

function hamiltonian(mol::Molecule)
    N = countorbitals(mol)
    H = fill(0.0eV, N, N)

    i = 1
    for atom1 in mol.atoms
        j = 1
        Ni = countorbitals(atom1)
        for atom2 in mol.atoms
            Nj = countorbitals(atom2)
            if atom1 === atom2
                H[i:i+Ni-1, j:j+Nj-1] = hamiltonian(atom1)
            elseif Bond(atom1, atom2) in mol.bonds
                H[i:i+Ni-1, j:j+Nj-1] = hamiltonian(atom1, atom2)
            end
            j += Nj
        end
        i += Ni
    end
    return H
end

function hamiltonian(mol1::Molecule, mol2::Molecule, bonds::Set{<:Bond})
    N1 = countorbitals(mol1)
    N2 = countorbitals(mol2)
    H = fill(0.0eV, N1, N2)

    i = 1
    for atom1 in mol1.atoms
        j = 1
        Ni = countorbitals(atom1)
        for atom2 in mol2.atoms
            Nj = countorbitals(atom2)
            if Bond(atom1, atom2) in bonds
                H[i:i+Ni-1, j:j+Nj-1] = hamiltonian(atom1, atom2)
            end
            j += Nj
        end
        i += Ni
    end
    return H
end

## Slater-Koster scheme

slaterkoster(f, δr::SVector{3}, s1::Orbital"*s", s2::Orbital"*s") = f(Val(:σ), norm(δr), s1, s2)

for (i, x) in [(1, :(:x)), (2, :(:y)), (3, :(:z))]
    @eval slaterkoster(f, δr::SVector{3}, s1::Orbital"*s", s2::Orbital{<:Any, :p, $x}) =
        normalize(δr)[$i] * f(Val(:σ), norm(δr), s1, s2)

    @eval slaterkoster(f, δr::SVector{3}, s1::Orbital{<:Any,:p,$x}, s2::Orbital{<:Any,:s,Nothing}) =
        normalize(δr)[$i] * f(Val(:σ), norm(δr), s1, s2)

    for (j, y) in [(1, :(:x)), (2, :(:y)), (3, :(:z))]
        if i == j
            @eval function slaterkoster(f, δr::SVector{3}, s1::Orbital{<:Any,:p,$x}, s2::Orbital{<:Any,:p,$x})
                d = normalize(δr)
                return d[$i]^2 * f(Val(:σ), norm(δr), s1, s2) + (1 - d[$i]^2) * f(Val(:π), norm(δr), s1, s2)
            end
        else
            @eval function slaterkoster(f, δr::SVector{3}, s1::Orbital{<:Any,:p,$x}, s2::Orbital{<:Any,:p,$y})
                d = normalize(δr)
                d[$i] * d[$j] * (f(Val(:σ), norm(δr), s1, s2) - f(Val(:π), norm(δr), s1, s2))
            end
        end
    end
end

include("energies.jl")

##################################################
# Angular momentum
##################################################

using Unitful: ħ

angularmomentum(vec::SVector{3}, args...) =
    vec[1]*angularmomentum(:x, args...) +
    vec[2]*angularmomentum(:y, args...) +
    vec[3]*angularmomentum(:z, args...)

angularmomentum(x::Symbol, args...) = angularmomentum(Val(x), args...)

angularmomentum(::Val, ::Orbital"*s", ::Orbital"*s") = 0.0ħ
angularmomentum(::Val, ::Orbital"*p_{*}", ::Orbital"*s") = 0.0ħ
angularmomentum(::Val, ::Orbital"*s", ::Orbital"*p_{*}") = 0.0ħ
angularmomentum(::Val, ::Orbital"*p_{*}", ::Orbital"*p_{*}") = 0.0ħ
angularmomentum(::Val{:x}, ::Orbital"*p_y", ::Orbital"*p_z") = -im*ħ
angularmomentum(::Val{:x}, ::Orbital"*p_z", ::Orbital"*p_y") = im*ħ
angularmomentum(::Val{:y}, ::Orbital"*p_z", ::Orbital"*p_x") = -im*ħ
angularmomentum(::Val{:y}, ::Orbital"*p_x", ::Orbital"*p_z") = im*ħ
angularmomentum(::Val{:z}, ::Orbital"*p_x", ::Orbital"*p_y") = -im*ħ
angularmomentum(::Val{:z}, ::Orbital"*p_y", ::Orbital"*p_x") = im*ħ

function angularmomentum(val::Val, atom::Atom)
    N = countorbitals(atom)
    L = fill(0.0im*ħ, N, N)
    for j in 1:N, i in 1:N
        L[i, j] = angularmomentum(val, atom.orbitals[i], atom.orbitals[j])
    end
    return L
end

function angularmomentum(val::Val, mol::Molecule)
    N = countorbitals(mol)
    L = fill(0.0im*ħ, N, N)

    i = 1
    for atom in mol.atoms
        Ni = countorbitals(atom)
        L[i:i+Ni-1, i:i+Ni-1] = angularmomentum(val, atom)
        i += Ni
    end

    return L
end

##################################################
# Distribution functions
##################################################

nB(E::Energy, T::Temperature = 300K) = 1 / (exp(E / (kB*T)) - 1)
nF(E::Energy, T::Temperature = 300K) = 1 / (exp(E / (kB*T)) + 1)

##################################################
# Leads
##################################################

# struct DensityOfStates{T}
#     f::T
# end

# (dos::DensityOfStates{<:AbstractMatrix})(E) = dos.f
# (dos::DensityOfStates{<:Function})(E) = dos.f(E)

"""
    Lead{N}

A lead with `N` channels.

"""
struct Lead{T}
    mol::Molecule{T}
    Γ::Energy{T}
end

"""
    selfenergy(l, mol)

Return the self-energy to `mol` from the lead, `l`.
"""

function selfenergy(l::Lead, mol::Molecule, bonds::Set{<:Bond})
    t = hamiltonian(l.mol, mol, bonds)
    H = hamiltonian(l.mol)
    return function (E::Energy)
        G = (E*I - H + I*im*l.Γ/2)^-1
        return t'G*t
    end
end

selfenergy(E::Energy, l::Lead, mol::Molecule, bonds::Set{<:Bond}) = selfenergy(l, mol, bonds)(E)

end #module
