#
# Bond energies taken from
# http://www.wiredchemist.com/chemistry/data/bond_energies_lengths.html
#

# Constants

const a0 = 4π * Unitful.ϵ0 * Unitful.ħ^2 / (Unitful.me * Unitful.q^2) |> Å

const HHssσ = -432 * Unitful.kJ / Unitful.mol / Unitful.Na |> eV

const HCssσ = HHssσ * exp(-(3 * 109u"pm" / 2 - 2 * 74u"pm") / a0)
const HCspσ = (- √3 * 411 * Unitful.kJ / Unitful.mol / Unitful.Na - HCssσ) / √2 |> eV

const CCssσ = HHssσ * exp(-(134u"pm" - 2 * 74u"pm") / a0)
const CCspσ = HCspσ * exp(-(134u"pm" - 3/2 * 109u"pm") / a0)
const CCppσ = (3 * (347 * Unitful.kJ / Unitful.mol / Unitful.Na |> eV) - CCssσ) / 2
const CCppπ = -(602 - 347) * Unitful.kJ / Unitful.mol / Unitful.Na |> eV

##################################################
# Onsite energies
##################################################

@inline onsiteenergy(::Hydrogen, ::Orbital"1s") = -13.59843449eV
@inline onsiteenergy(::Carbon, ::Orbital"2s") = -15.9625933eV
@inline onsiteenergy(::Carbon, ::Orbital"2p_{*}") = -11.2602880eV

##################################################
# Bond energies
##################################################

τ = 0//1

@inline bondenergy(::Val{:σ}, δr, ::Hydrogen, ::Hydrogen, ::Orbital"1s", ::Orbital"1s") = 2*HHssσ * exp(-(δr - 1Å)/a0*τ)

@inline bondenergy(::Val{:σ}, δr, ::Hydrogen, ::Carbon, ::Orbital"1s", ::Orbital"2s") = 2*HCssσ * exp(-8(δr - 1Å)/a0*τ)
@inline bondenergy(::Val{:σ}, δr, ::Carbon, ::Hydrogen, ::Orbital"2s", ::Orbital"1s") = 2*HCssσ * exp(-8(δr - 1Å)/a0*τ)
@inline bondenergy(::Val{:σ}, δr, ::Hydrogen, ::Carbon, ::Orbital"1s", ::Orbital"2p_{*}") = HCspσ * exp(-8(δr - 1Å)/a0*τ)
@inline bondenergy(::Val{:σ}, δr, ::Carbon, ::Hydrogen, ::Orbital"2p_{*}", ::Orbital"1s") = -HCspσ * exp(-8(δr - 1Å)/a0*τ)

@inline bondenergy(::Val{:σ}, δr, ::Carbon, ::Carbon, ::Orbital"2s", ::Orbital"2s") = 2*CCssσ * exp(-12(δr - 1Å)/a0*τ)
@inline bondenergy(::Val{:σ}, δr, ::Carbon, ::Carbon, ::Orbital"2s", ::Orbital"2p_{*}") = 2*CCspσ * exp(-12(δr - 1Å)/a0*τ)
@inline bondenergy(::Val{:σ}, δr, ::Carbon, ::Carbon, ::Orbital"2p_{*}", ::Orbital"2s") = -CCspσ * exp(-12(δr - 1Å)/a0*τ)
@inline bondenergy(::Val{:σ}, δr, ::Carbon, ::Carbon, ::Orbital"2p_{*}", ::Orbital"2p_{*}") = CCppσ * exp(-12(δr - 1Å)/a0*τ)
@inline bondenergy(::Val{:π}, δr, ::Carbon, ::Carbon, ::Orbital"2p_{*}", ::Orbital"2p_{*}") = CCppπ * exp(-12(δr - 1Å)/a0*τ)

##################################################
# Slater-Koster scheme
##################################################

function slaterkoster(f, δr::SVector{3}, atom1::Atom, atom2::Atom, s1::Orbital"*s", s2::Orbital"*s")
    return f(Val(:σ), norm(δr), atom1, atom2, s1, s2)
end

for (i, x) in [(1, :(:x)), (2, :(:y)), (3, :(:z))]
    @eval @inline function slaterkoster(f, δr::SVector{3}, atom1::Atom, atom2::Atom, s1::Orbital"*s", s2::Orbital{<:Any, :p, $x})
        d = normalize(δr)
        return d[$i] * f(Val(:σ), norm(δr), atom1, atom2, s1, s2)
    end

    @eval @inline function slaterkoster(f, δr::SVector{3}, atom1::Atom, atom2::Atom, s1::Orbital{<:Any,:p,$x}, s2::Orbital{<:Any,:s,Nothing})
        d = normalize(δr)
        return d[$i] * f(Val(:σ), norm(δr), atom1, atom2, s1, s2)
    end

    for (j, y) in [(1, :(:x)), (2, :(:y)), (3, :(:z))]
        if i == j
            @eval @inline function slaterkoster(f, δr::SVector{3}, atom1::Atom, atom2::Atom, s1::Orbital{<:Any,:p,$x}, s2::Orbital{<:Any,:p,$x})
                d = normalize(δr)
                return d[$i]^2 * f(Val(:σ), norm(δr), atom1, atom2, s1, s2) + (1 - d[$i]^2) * f(Val(:π), norm(δr), atom1, atom2, s1, s2)
            end
        else
            @eval @inline function slaterkoster(f, δr::SVector{3}, atom1::Atom, atom2::Atom, s1::Orbital{<:Any,:p,$x}, s2::Orbital{<:Any,:p,$y})
                d = normalize(δr)
                return d[$i] * d[$j] * (f(Val(:σ), norm(δr), atom1, atom2, s1, s2) - f(Val(:π), norm(δr), atom1, atom2, s1, s2))
            end
        end
    end
end
