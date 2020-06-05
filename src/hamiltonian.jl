##################################################
# Hamiltonians
##################################################

"""
    hamiltonian([T::Type], f, args...; kwargs...)

Construct the matrix corresponding to the hamiltonian operator of an
atom, bond or molecule.  `f` is a function that accepts the atom types
and distances between the atoms to compute the hamiltonian matrix
element.

"""
hamiltonian(f, args...; kwargs...) = hamiltonian(Complex{Float64}, f, args...; kwargs...)

function hamiltonian(T::Type, f, atom::Atom; electric_field=(x, y, z) -> 0.0, kwargs...)
    N = countorbitals(atom)
    V = electric_field(position(atom)...)
    #f = 0
    H = fill(zero(T), N, N)
    U = localbasis(atom)
    for (i, s1) in enumerate(atom.orbitals)
        H[i,i] += f(chemical_symbol(atom) => shortform(s1))
    end
    return U'H*U - I*V
end

function hamiltonian(T::Type, f, atom1::Atom, atom2::Atom; fixed_distance=nothing, kwargs...)
    N1 = countorbitals(atom1)
    N2 = countorbitals(atom2)
    U1 = localbasis(atom1)
    U2 = localbasis(atom2)
    H = fill(zero(T), N1, N2)

    if fixed_distance === nothing
        r = atom1.position - atom2.position
    else
        r = fixed_distance .* normalize(atom1.position - atom2.position)
    end

    for (j, s2) in enumerate(atom2.orbitals), (i, s1) in enumerate(atom1.orbitals)
        H[i,j] += f(r, chemical_symbol(atom1) => shortform(s1), chemical_symbol(atom2) => shortform(s2))
    end
    return U1'H*U2
end

function hamiltonian(T::Type, f, mol::Molecule; kwargs...)
    N = countorbitals(mol)
    H = fill(zero(T), N, N)

    i = 1
    for atom1 in mol.atoms
        j = 1
        Ni = countorbitals(atom1)
        for atom2 in mol.atoms
            Nj = countorbitals(atom2)
            if atom1 === atom2
                H[i:i+Ni-1, j:j+Nj-1] += hamiltonian(T, f, atom1; kwargs...)
            elseif Bond(atom1, atom2) in mol.bonds
                H[i:i+Ni-1, j:j+Nj-1] += hamiltonian(T, f, atom1, atom2; kwargs...)
            end
            j += Nj
        end
        i += Ni
    end
    return H
end

function hamiltonian(T::Type, f, mol1::Molecule, mol2::Molecule, bonds::Set{<:Bond}; kwargs...)
    N1 = countorbitals(mol1)
    N2 = countorbitals(mol2)
    H = fill(zero(T), N1, N2)

    i = 1
    for atom1 in mol1.atoms
        j = 1
        Ni = countorbitals(atom1)
        for atom2 in mol2.atoms
            Nj = countorbitals(atom2)
            if Bond(atom1, atom2) in bonds
                H[i:i+Ni-1, j:j+Nj-1] += hamiltonian(T, f, atom1, atom2; kwargs...)
            end
            j += Nj
        end
        i += Ni
    end
    return H
end
