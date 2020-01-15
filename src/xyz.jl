function write_xyz(mol::Molecule, filename)
    lines = String[]
    push!(lines, string(length(mol.atoms)))
    push!(lines, "Properties=species:S:1:pos:R:3")
    for atom in mol.atoms
        push!(lines, join([chemical_symbol(atom), position(atom)...], "\t"))
    end
    Base.write(filename, join(lines, "\n"))
end

function read_xyz(filename)
    lines = readlines(filename)
    N = parse(Int, lines[1])
    properties = split(lines[2])
    atoms = Atom[]
    for (idx, line) in enumerate(lines[3:end])
        args = split(line)
        sym = Symbol(args[1])
        pos = map(x -> parse(Float64, x), args[2:4])
        push!(atoms, chemical_symbol_to_atom_type(sym)(pos))
    end
    @assert N == length(atoms)
    return N, atoms
end

function read_xyz!(mol::Molecule, filename)
    N, atoms = read_xyz(filename)
    @assert length(mol) == N
    for (atom1, atom2) in zip(mol.atoms, atoms)
        @assert chemical_symbol(atom1) == chemical_symbol(atom2)
        atom1.position = atom2.position
    end
    return mol
end
