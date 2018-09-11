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

onsiteenergy(::Hydrogen, ::Orbital"1s") = -13.59843449eV
onsiteenergy(::Carbon, ::Orbital"2s") = -15.9625933eV
onsiteenergy(::Carbon, ::Orbital"2p_{*}") = -11.2602880eV

##################################################
# Bond energies
##################################################

bondenergy(::Val{:σ}, δr, ::Hydrogen, ::Hydrogen, ::Orbital"1s", ::Orbital"1s") = HHssσ

bondenergy(::Val{:σ}, δr, ::Hydrogen, ::Carbon, ::Orbital"1s", ::Orbital"2s") = HCssσ
bondenergy(::Val{:σ}, δr, ::Carbon, ::Hydrogen, ::Orbital"2s", ::Orbital"1s") = HCssσ
bondenergy(::Val{:σ}, δr, ::Hydrogen, ::Carbon, ::Orbital"1s", ::Orbital"2p_{*}") = HCspσ
bondenergy(::Val{:σ}, δr, ::Carbon, ::Hydrogen, ::Orbital"2p_{*}", ::Orbital"1s") = -HCspσ

bondenergy(::Val{:σ}, δr, ::Carbon, ::Carbon, ::Orbital"2s", ::Orbital"2s") = CCssσ
bondenergy(::Val{:σ}, δr, ::Carbon, ::Carbon, ::Orbital"2s", ::Orbital"2p_{*}") = CCspσ
bondenergy(::Val{:σ}, δr, ::Carbon, ::Carbon, ::Orbital"2p_{*}", ::Orbital"2s") = -CCspσ
bondenergy(::Val{:σ}, δr, ::Carbon, ::Carbon, ::Orbital"2p_{*}", ::Orbital"2p_{*}") = CCppσ
bondenergy(::Val{:π}, δr, ::Carbon, ::Carbon, ::Orbital"2p_{*}", ::Orbital"2p_{*}") = CCppπ
