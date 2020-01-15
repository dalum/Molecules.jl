radius(::Hydrogen) = 0.2
radius(::Carbon) = 0.3
radius(::Nitrogen) = 8
radius(::Oxygen) = 9

function writetoblenderscript(mol::Molecule, filename)
    lines = String[]
    push!(lines, "import bpy", "import numpy")

    for element in ("Hydrogen", "Carbon", "Nitrogen", "Oxygen")
        push!(lines, """$(lowercase(element))_material = bpy.data.materials.new(name='$(element)Material')""")
    end

    for atom in mol.atoms
        x, y, z = position(atom)
        push!(lines, """bpy.ops.mesh.primitive_uv_sphere_add(radius=$(radius(atom)), location=($x, $y, $z))""")
    end

    for bond in mol.bonds
        x, y, z = (position(bond.atom2) .+ position(bond.atom1)) ./ 2
        δ = position(bond.atom2) .- position(bond.atom1)
        d = norm(δ)
        #z += d / 2
        θ = π/2 - atan(δ[3], sqrt(δ[1]^2 + δ[2]^2))
        ϕ = atan(δ[2], δ[1])

        push!(lines, """bpy.ops.mesh.primitive_cylinder_add(depth=$(norm(δ)), radius=0.1, location=($x, $y, $z), rotation=(0.0, $θ, $ϕ))""")
    end

    script = join(lines, "\n")
    write(filename, script)
end
