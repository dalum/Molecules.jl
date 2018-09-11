module TestMolecules

using Test
using LinearAlgebra
using Molecules
using Unitful: eV, ħ
using StaticArrays: SVector
using Molecules: Bond, Carbon, Hydrogen, Lead, Molecule, Orbital,
    angularmomentum, chemicalpotential, countorbitals, hamiltonian, onsiteenergy, selfenergy

##################################################
# Setup
##################################################

atom1 = Hydrogen(.0Å, .0Å, .0Å, orbital"1s")
atom2 = Carbon(1.0Å, .0Å, .0Å, orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z")
atom3 = Hydrogen(1.0Å, .0Å, 1.0Å, orbital"1s")
atom4 = Carbon(1.0Å, .0Å, 1.0Å, orbital"2s", orbital"2p_x", orbital"2p_y", orbital"2p_z")

bond12 = Bond(atom1, atom2)
bond23 = Bond(atom2, atom3)
mol1 = Molecule([atom1, atom2, atom3], [bond12, bond23])

h1 = hamiltonian(atom1)
h2 = hamiltonian(atom2)
h3 = hamiltonian(atom3)
t12 = hamiltonian(atom1, atom2)
t23 = hamiltonian(atom2, atom3)
t13 = hamiltonian(atom1, atom3)

H = [h1 t12 0t13
     t12' h2 t23
     0t13' t23' h3]


Lx = [0.0ħ  0.0ħ  0.0ħ  0.0ħ
      0.0ħ  0.0ħ  0.0ħ  0.0ħ
      0.0ħ  0.0ħ  0.0ħ -im*ħ
      0.0ħ  0.0ħ  im*ħ  0.0ħ]

Ly = [0.0ħ  0.0ħ  0.0ħ  0.0ħ
      0.0ħ  0.0ħ  0.0ħ  im*ħ
      0.0ħ  0.0ħ  0.0ħ  0.0ħ
      0.0ħ -im*ħ  0.0ħ  0.0ħ]

Lz = [0.0ħ  0.0ħ  0.0ħ  0.0ħ
      0.0ħ  0.0ħ -im*ħ  0.0ħ
      0.0ħ  im*ħ  0.0ħ  0.0ħ
      0.0ħ  0.0ħ  0.0ħ  0.0ħ]

vec = SVector(randn(3)...)

z12 = fill(0.0im*ħ, countorbitals(atom1), countorbitals(atom2))
z13 = fill(0.0im*ħ, countorbitals(atom1), countorbitals(atom3))
z23 = fill(0.0im*ħ, countorbitals(atom2), countorbitals(atom3))

atom6 = Hydrogen(-1.0Å, 0.0Å, 0.0Å, orbital"1s")
mol2 = Molecule([atom6], Bond[])
lead1 = Lead(mol2, 100.0eV)

μ1, Δ1 = chemicalpotential(mol1)
Σ1 = selfenergy(μ1 + Δ1, lead1, mol1, Set([Bond(atom1, atom6)]))
H1 = hamiltonian(mol1)

##################################################
# Tests
##################################################

@testset "Orbitals" begin
    @test orbital"1s" isa Orbital{1,:s,Nothing}
    @test orbital"1s" isa Orbital"*s_{*}"
    @test orbital"2p_x" ≡ orbital"2p_{x}"
    @test_throws LoadError eval(@macroexpand orbital"p_x")
    @test_throws LoadError eval(@macroexpand orbital"2p_x^2")
    @test Orbital"1s" ≡ Orbital{1,:s,Nothing}
    @test Orbital"2p_y" ≡ Orbital{2,:p,:y}
end

@testset "Atoms" begin
    @test countorbitals(atom1) == 1
    @test countorbitals(atom2) == 4
    @test hamiltonian(atom1) == fill(onsiteenergy(atom1, orbital"1s"), (1, 1))
    @test hamiltonian(atom2) == [onsiteenergy(atom2, orbital"2s") 0.0eV 0.0eV 0.0eV
                                 0.0eV onsiteenergy(atom2, orbital"2p_x") 0.0eV 0.0eV
                                 0.0eV 0.0eV onsiteenergy(atom2, orbital"2p_y") 0.0eV
                                 0.0eV 0.0eV 0.0eV onsiteenergy(atom2, orbital"2p_z")]
end

@testset "Molecules" begin
    @test countorbitals(mol1) == 6
    @test hamiltonian(mol1) == H
end

@testset "Energies" begin
    # Onsite energies are negative (electrons are bound)

    @test all(x -> x ≤ 0.0eV, hamiltonian(atom1))

    t12 = hamiltonian(atom2, atom4)

    # Check that energies match the expected bonding/anti--bonding
    # behaviour

    # s - *
    @test t12[1,1] < 0.0eV
    @test t12[1,2] ≈ 0.0eV
    @test t12[1,3] ≈ 0.0eV
    # s - p_z
    @test t12[1,4] > 0.0eV
    @test t12[4,1] < 0.0eV

    # p_x - *
    @test t12[2,1] ≈ 0.0eV
    @test t12[2,2] < 0.0eV
    @test t12[2,3] ≈ 0.0eV
    @test t12[2,4] ≈ 0.0eV
end

@testset "Angular momentum" begin
    @test angularmomentum(:z, orbital"1p_x", orbital"1p_y") == -im*ħ
    @test angularmomentum(:z, orbital"1p_y", orbital"1p_x") == im*ħ
    @test angularmomentum(:x, atom1) == fill(0.0ħ, 1, 1)
    @test angularmomentum(:x, atom2) == Lx
    @test angularmomentum(:y, atom2) == Ly
    @test angularmomentum(:z, atom2) == Lz
    @test angularmomentum(vec, atom2) == vec[1]*Lx + vec[2]*Ly + vec[3]*Lz
    @test angularmomentum(vec, mol1) == [angularmomentum(vec, atom1) z12 z13
                                        z12' angularmomentum(vec, atom2) z23
                                        z13' z23' angularmomentum(vec, atom3)]
end

@testset "Leads" begin
    @test size(H1) == size(Σ1)
    @test Σ1[1,1] ≠ 0.0eV
    Σ1′ = copy(Σ1)
    Σ1′[1, 1] = 0.0eV
    @test all(x -> x == 0.0eV, Σ1′)
end

end #module
