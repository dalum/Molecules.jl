struct UnitCell
    unit::Molecule
    nextunit::Molecule
    bonds::Set{Bond}
end

function hamiltonian(T::Type, f, uc::UnitCell)
    mol = uc.unit
    mol′ = uc.nextunit
    N = countorbitals(mol)
    H0 = fill(zero(T), N, N)
    t1 = fill(zero(T), N, N)

    i = 1
    for (atom1, atom1′) in zip(mol.atoms, mol′.atoms)
        Ni = countorbitals(atom1)

        j = 1
        for (atom2, atom2′) in zip(mol.atoms, mol′.atoms)
            Nj = countorbitals(atom2)

            if atom1 === atom2
                H0[i:i+Ni-1, j:j+Nj-1] += hamiltonian(T, f, atom1)
            elseif Bond(atom1, atom2) in mol.bonds
                H0[i:i+Ni-1, j:j+Nj-1] += hamiltonian(T, f, atom1, atom2)
            end

            if Bond(atom1, atom2′) in uc.bonds
                t1[i:i+Ni-1, j:j+Nj-1] += hamiltonian(T, f, atom1, atom2′)
            end

            j += Nj
        end

        i += Ni
    end
    return H0, t1
end

function overlap(T::Type, f, uc::UnitCell)
    mol = uc.unit
    mol′ = uc.nextunit
    N = countorbitals(mol)
    S0 = fill(zero(T), N, N)
    S1 = fill(zero(T), N, N)

    i = 1
    for (atom1, atom1′) in zip(mol.atoms, mol′.atoms)
        Ni = countorbitals(atom1)

        j = 1
        for (atom2, atom2′) in zip(mol.atoms, mol′.atoms)
            Nj = countorbitals(atom2)

            if atom1 === atom2
                S0[i:i+Ni-1, j:j+Nj-1] += overlap(T, f, atom1)
            elseif Bond(atom1, atom2) in mol.bonds
                S0[i:i+Ni-1, j:j+Nj-1] += overlap(T, f, atom1, atom2)
            end

            if Bond(atom1, atom2′) in uc.bonds
                S1[i:i+Ni-1, j:j+Nj-1] += overlap(T, f, atom1, atom2′)
            end

            j += Nj
        end

        i += Ni
    end
    return S0, S1
end
