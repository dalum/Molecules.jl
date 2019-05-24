module Molecules

using DataStructures: OrderedSet
using ForwardDiff
using LinearAlgebra
using SlaterKoster
using StaticArrays: SVector, SMatrix

import SlaterKoster: SlaterKosterTable, hamiltonian

export @Orbital_str, @orbital_str,
    chemicalpotential, hamiltonian, makeaxes, molecularorbitals, selfenergy, slaterkoster, vibrationalcoupling

const σ = Dict{Union{Symbol,Int},Matrix{Complex{Float64}}}()
const σ[0] = [1.0 0.0; 0.0 1.0]
const σ[:x] = [0.0 1.0; 1.0 0.0]
const σ[:y] = [0.0 -1.0im; 1.0im 0.0]
const σ[:z] = [1.0 0.0; 0.0 -1.0]
const ⊗ = kron

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

function shortform(::Orbital{N,L,M}) where {N,L,M}
    l = L
    m = M == Nothing ? "" : replace(string(M), r"[_^+\-{}]" => s"")
    return Symbol("$l$m")
end

mutable struct Atom{T1<:AbstractVector,T2<:AbstractMatrix,SYMBOL,Z,N}
    position::T1
    axes::T2
    orbitals::OrderedSet{Orbital}

    function Atom{T1,T2,SYMBOL,Z,N}(position::T1, axes::T2, orbitals::Orbital...) where {T1<:AbstractVector,T2<:AbstractMatrix,SYMBOL,Z,N}
        return new{T1,T2,SYMBOL,Z,N}(
            position,
            axes,
            OrderedSet{Orbital}(orbitals))
    end
end

function Atom{<:Any,<:Any,SYMBOL,Z,N}(position::T1, orbitals::Orbital...; axes=makeaxes(0., 0.)) where {T1,SYMBOL,Z,N}
    return Atom{T1,typeof(axes),SYMBOL,Z,N}(position, axes, orbitals...)
end

Base.getindex(atom::Atom, i::Integer) = collect(atom.orbitals)[i]

position(atom::Atom{<:AbstractVector}) = atom.position
position(atom::Atom{<:AbstractVector{<:ForwardDiff.Dual}}) = (x -> x.value).(atom.position)

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

mutable struct Molecule{T1<:AbstractVector,T2<:AbstractMatrix}
    position::T1
    axes::T2
    atoms::OrderedSet{Atom}
    bonds::Set{Bond}
end

Molecule(atoms, bonds) = Molecule(
    SVector{3}(0., 0., 0.),
    makeaxes(0., 0.),
    OrderedSet{Atom}(atoms),
    Set(bonds))

Molecule() = Molecule(
    SVector{3}(0., 0., 0.),
    makeaxes(0., 0.),
    OrderedSet{Atom}(),
    Set{Bond}())

Base.push!(mol::Molecule, atoms::Atom...) = push!(mol.atoms, atoms...)
Base.push!(mol::Molecule, bonds::Bond...) = push!(mol.bonds, bonds...)

Base.eachindex(mol::Molecule) = eachindex(collect(mol.atoms))
Base.iterate(mol::Molecule) = iterate(collect(mol.atoms))
Base.iterate(mol::Molecule, i) = iterate(collect(mol.atoms), i)
Base.getindex(mol::Molecule, i::Integer) = collect(mol.atoms)[i]
Base.lastindex(mol::Molecule) = lastindex(collect(mol.atoms))

function indices(mol::Molecule, atom::Atom)
    N0 = 1
    for atom′ in mol.atoms
        N1 = N0 + countorbitals(atom′)
        if atom === atom′
            return N0:(N1-1)
        end
        N0 = N1
    end
    return nothing
end

function indices(mol::Molecule, atom::Atom, orbital::Orbital)
    N0 = 1
    for atom′ in mol.atoms
        N1 = N0 + countorbitals(atom′)
        if atom === atom′
            for (i, orbital′) in enumerate(atom.orbitals)
                if orbital == orbital′
                    return (N0+i-1):(N0+i-1)
                end
            end
            return nothing
        end
        N0 = N1
    end
    return nothing
end

Base.Vector{T}(mol::Molecule) where {T} = normalize!(fill(oneunit(T), countorbitals(mol)))

function Base.Vector{T}(mol::Molecule, atom::Atom) where {T}
    u = fill(zero(T), countorbitals(mol))
    N = countorbitals(atom)
    u[indices(mol, atom)] .= 1 / sqrt(N)
    return u
end

function Base.Vector{T}(mol::Molecule, n::Integer) where {T}
    atom = mol[n]
    return Vector{T}(mol, atom)
end

function Base.Vector{T}(mol::Molecule, atom::Atom, orbital::Orbital) where {T}
    u = fill(zero(T), countorbitals(mol))
    u[indices(mol, atom, orbital)] .= 1
    return u
end

function Base.Vector{T}(mol::Molecule, atom::Atom, s::Integer) where {T}
    orbital = atom[s]
    return Vector{T}(mol, atom, orbital)
end

function Base.Vector{T}(mol::Molecule, n::Integer, s::Integer) where {T}
    atom = mol[n]
    orbital = atom[s]
    return Vector{T}(mol, atom, orbital)
end

##################################################
# Chemical properties
##################################################

protonnumber(::Atom{T1,T2,SYMBOL,Z,N}) where {T1,T2,SYMBOL,Z,N} = Z
neutronnumber(::Atom{T1,T2,SYMBOL,Z,N}) where {T1,T2,SYMBOL,Z,N} = N
symbol(::Atom{T1,T2,SYMBOL,Z,N}) where {T1,T2,SYMBOL,Z,N} = SYMBOL

countorbitals(atom::Atom)::Int = length(atom.orbitals)
countorbitals(mol::Molecule)::Int = sum(map(countorbitals, mol.atoms))

mass(skt::SlaterKosterTable, atom::Atom) = SlaterKoster.mass(skt, symbol(atom))

function electronvalence(atom::Atom)
    n = protonnumber(atom)
    n <= 2 && return n
    n <= 10 && return n - 2
    n <= 18 && return n - 10 # Not very generic ...
    n > 18 && error("atom numbers greater than 18 are currently not supported")
end
"""
    valencecount(mol)

Return the sum of electron valences of each atom in `mol`.

"""
valencecount(mol::Molecule, oxidation=0) = sum(electronvalence(atom) for atom in mol.atoms) - oxidation

function molecularorbitals(H)
    vals, vecs = eigen(H)
    s = sortperm(vals)
    return vals[s], vecs[:, s]
end

function chemicalpotential(mol::Molecule, levels, oxidation=0)
    n = valencecount(mol, oxidation)
    μ = real(levels[ceil(Int, n/2)] + levels[ceil(Int, (n + 1)/2)])/2
    Δ = real(levels[ceil(Int, (n + 1)/2)] - levels[ceil(Int, n/2)])/2
    return μ, Δ
end

include("hamiltonian.jl")

include("unitcell.jl")

##################################################
# Elements
##################################################

const Hydrogen = Atom{<:Any, <:Any, :H, 1, 0}
const Deuterium = Atom{<:Any, <:Any, :D, 1, 1}
const Tritium = Atom{<:Any, <:Any, :T, 1, 2}
const Carbon = Atom{<:Any, <:Any, :C, 6, 6}
const Carbon13 = Atom{<:Any, <:Any, :C, 6, 7}
const Nitrogen = Atom{<:Any, <:Any, :N, 7, 7}
const Oxygen = Atom{<:Any, <:Any, :O, 8, 8}

##################################################
# Coordinates
##################################################

function makeaxes(θ = 0.0, ϕ = 0.0, ρ = 0.0)
    R1 = LinearAlgebra.Givens(2, 1, cos(ρ), sin(ρ))
    R2 = LinearAlgebra.Givens(3, 1, cos(ϕ), sin(ϕ))
    R3 = LinearAlgebra.Givens(2, 1, cos(θ), sin(θ))
    return SMatrix{3, 3}(R3*(R2*(R1*Matrix(I, 3, 3))))
end

##################################################
# Basis change
##################################################

localbasis(atom::Atom) = localbasis(atom.axes, atom.orbitals...)

function localbasis(axes::AbstractMatrix{T}, ::Orbital"*s") where {T}
    return fill(one(T), 1, 1)
end

function localbasis(axes::AbstractMatrix{T}, ::Orbital"*s", ::Orbital"*p_x", ::Orbital"*p_y", ::Orbital"*p_z") where {T}
    #return one(Matrix{T}(undef, 4, 4))
    Z = zero(T) / oneunit(T)
    return [one(T) fill(Z, 1, 3)
            fill(Z, 3, 1) axes]
end

##################################################
# Position operator
##################################################

_position(::Val{:x}, atom::Atom) = atom.position[1]
_position(::Val{:y}, atom::Atom) = atom.position[2]
_position(::Val{:z}, atom::Atom) = atom.position[3]

position(x::Symbol, args...) = position(Val(x), args...)

position(::Val, s1::Orbital, s2::Orbital) = s1 == s2 ? 1.0 : 0.0

function position(val::Val, atom::Atom)
    N = countorbitals(atom)
    U = localbasis(atom)
    X = fill(0.0im, N, N)
    for (j, s2) in enumerate(atom.orbitals), (i, s1) in enumerate(atom.orbitals)
        X[i, j] = _position(val, atom) * position(val, s1, s2)
    end
    return U'X*U
end

function position(val::Val, mol::Molecule)
    N = countorbitals(mol)
    X = fill(0.0im, N, N)

    i = 1
    for atom in mol.atoms
        Ni = countorbitals(atom)
        X[i:i+Ni-1, i:i+Ni-1] = position(val, atom)
        i += Ni
    end

    return X
end

# function position(val::Val, mol::Molecule, n)
#     N = countorbitals(mol)
#     X = fill(0.0im, N, N)
#     inds = indices(mol, mol[n])
#     X[inds, inds] = position(val, mol[n])
#     return X
# end

##################################################
# Angular momentum
##################################################

angularmomentum(vec::SVector{3}, args...) =
    vec[1]*angularmomentum(:x, args...) +
    vec[2]*angularmomentum(:y, args...) +
    vec[3]*angularmomentum(:z, args...)

angularmomentum(x::Symbol, args...) = angularmomentum(Val(x), args...)

angularmomentum(::Val, ::Orbital"*s", ::Orbital"*s") = 0.0
angularmomentum(::Val, ::Orbital"*p_{*}", ::Orbital"*s") = 0.0
angularmomentum(::Val, ::Orbital"*s", ::Orbital"*p_{*}") = 0.0
angularmomentum(::Val, ::Orbital"*p_{*}", ::Orbital"*p_{*}") = 0.0
angularmomentum(::Val{:x}, ::Orbital"*p_y", ::Orbital"*p_z") = -im
angularmomentum(::Val{:x}, ::Orbital"*p_z", ::Orbital"*p_y") = im
angularmomentum(::Val{:y}, ::Orbital"*p_z", ::Orbital"*p_x") = -im
angularmomentum(::Val{:y}, ::Orbital"*p_x", ::Orbital"*p_z") = im
angularmomentum(::Val{:z}, ::Orbital"*p_x", ::Orbital"*p_y") = -im
angularmomentum(::Val{:z}, ::Orbital"*p_y", ::Orbital"*p_x") = im

function angularmomentum(val::Val, atom::Atom)
    N = countorbitals(atom)
    U = localbasis(atom)
    L = fill(0.0im, N, N)
    for (j, s2) in enumerate(atom.orbitals), (i, s1) in enumerate(atom.orbitals)
        L[i, j] = angularmomentum(val, s1, s2)
    end
    return U'L*U
end

function angularmomentum(val::Val, mol::Molecule)
    N = countorbitals(mol)
    L = fill(0.0im, N, N)

    i = 1
    for atom in mol.atoms
        Ni = countorbitals(atom)
        L[i:i+Ni-1, i:i+Ni-1] = angularmomentum(val, atom)
        i += Ni
    end

    return L
end

function angularmomentum(val::Val, mol::Molecule, n)
    N = countorbitals(mol)
    L = fill(0.0im, N, N)
    inds = indices(mol, mol[n])
    L[inds, inds] = angularmomentum(val, mol[n])
    return L
end

##################################################
# 
##################################################

# function inlocalbasis(x::Int, op, atom::Atom)
#     N = countorbitals(atom)
#     U = localbasis(atom)
#     op = U*op
#     L = fill(0.0im, N, N)
#     for (j, s2) in enumerate(atom.orbitals), (i, s1) in enumerate(atom.orbitals)
#         L[i, j] = angularmomentum(val, s1, s2)
#     end
#     return U'L*U
# end


# function inlocalbasis(x::Int, op, mol::Molecule)
    
# end

##################################################
# Distribution functions
##################################################

# In atomic units, room temperature is 300K / 315774K = 9.5e-4
nB(E, T = 9.5e-4) = 1/(exp(E/T) - 1)
nF(E, T = 9.5e-4) = 1/(exp(E/T) + 1)

##################################################
# Vibrations
##################################################

"""
    vibrationalcoupling(mol, atom, ω0, dr)

Calculate the numerical derivative of the `hamiltonian(mol)` w.r.t. a
displacement along `dr` of the position of `atom`.  `ω0` is the
characteristic frequency of the atomic vibration.

"""
function vibrationalcoupling(mol::Molecule, atom::Atom, ω0, dr)
    m = Molecules.mass(atom)
    H0 = hamiltonian(mol)
    atom.position += dr
    H′ = hamiltonian(mol)
    atom.position -= dr
    return (H′ - H0) / norm(dr) * √(1 / (2 * m * ω0))
end

##################################################
# Leads
##################################################

"""
    Lead{N}

A lead with `N` channels.

"""
struct Lead{T}
    mol::Molecule
    Γ::T
end

"""
    selfenergy(l, mol)

Return the self-energy to `mol` from the lead, `l`.
"""

function selfenergy(l::Lead, H::AbstractMatrix, t::AbstractMatrix)
    return function (E)
        G = inv(E*I - H + I*im*l.Γ/2)
        return t'G*t
    end
end

function selfenergy(T, skt::SlaterKosterTable, l::Lead, mol::Molecule, bonds::Set{<:Bond}; f=identity)
    t = f.(hamiltonian(T, skt, l.mol, mol, bonds))
    H = f.(hamiltonian(T, skt, l.mol))
    return selfenergy(l, H, t)
end

##################################################
# Plots recipes
##################################################

include("recipe.jl")
include("blender.jl")

end #module
