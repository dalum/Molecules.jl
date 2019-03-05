struct UnitCell
    unit::Molecule
    nextunit::Molecule
    bonds::Set{Bond}
end

function hamiltonian(T::Type, skt::SlaterKosterTable, uc::UnitCell, ξ)
    mol = uc.unit
    mol′ = uc.nextunit
    N = countorbitals(mol)
    H = fill(zero(T), N, N)

    i = 1
    for (atom1, atom1′) in zip(mol.atoms, mol′.atoms)
        Ni = countorbitals(atom1)

        j = 1
        for (atom2, atom2′) in zip(mol.atoms, mol′.atoms)
            Nj = countorbitals(atom2)

            if atom1 === atom2
                H[i:i+Ni-1, j:j+Nj-1] += hamiltonian(T, skt, atom1)
            elseif Bond(atom1, atom2) in mol.bonds
                H[i:i+Ni-1, j:j+Nj-1] += hamiltonian(T, skt, atom1, atom2)
            end

            if Bond(atom1, atom2′) in uc.bonds
                H[i:i+Ni-1, j:j+Nj-1] += hamiltonian(T, skt, atom1, atom2′)*exp(-im*ξ)
            end
            if Bond(atom1′, atom2) in uc.bonds
                H[i:i+Ni-1, j:j+Nj-1] += hamiltonian(T, skt, atom1′, atom2)*exp(im*ξ)
            end

            j += Nj
        end

        i += Ni
    end
    return H
end
