/**
 *  @file MaskellSolidSolnPhase.cpp Implementation file for an ideal solid
 *      solution model with incompressible thermodynamics (see @ref
 *      thermoprops and @link Cantera::MaskellSolidSolnPhase
 *      MaskellSolidSolnPhase@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/MaskellSolidSolnPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"

#include <cassert>

namespace Cantera
{

MaskellSolidSolnPhase::MaskellSolidSolnPhase()
{
    warn_deprecated("class MaskellSolidSolnPhase", "To be removed after Cantera 3.0");
}

void MaskellSolidSolnPhase::getActivityConcentrations(CanteraDouble* c) const
{
    getActivityCoefficients(c);
    for (size_t sp = 0; sp < m_kk; ++sp) {
        c[sp] *= moleFraction(sp);
    }
}

// Molar Thermodynamic Properties of the Solution

CanteraDouble MaskellSolidSolnPhase::enthalpy_mole() const
{
    const CanteraDouble h0 = RT() * mean_X(m_h0_RT);
    const CanteraDouble r = moleFraction(product_species_index);
    const CanteraDouble fmval = fm(r);
    return h0 + r * fmval * h_mixing;
}

CanteraDouble xlogx(CanteraDouble x)
{
    return x * std::log(x);
}

CanteraDouble MaskellSolidSolnPhase::entropy_mole() const
{
    const CanteraDouble s0 = GasConstant * mean_X(m_s0_R);
    const CanteraDouble r = moleFraction(product_species_index);
    const CanteraDouble fmval = fm(r);
    const CanteraDouble rfm = r * fmval;
    return s0 + GasConstant * (xlogx(1-rfm) - xlogx(rfm) - xlogx(1-r-rfm) - xlogx((1-fmval)*r) - xlogx(1-r) - xlogx(r));
}

// Mechanical Equation of State

void MaskellSolidSolnPhase::calcDensity()
{
    const vector<CanteraDouble>& vbar = getStandardVolumes();

    vector<CanteraDouble> moleFracs(m_kk);
    Phase::getMoleFractions(&moleFracs[0]);
    CanteraDouble vtotal = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        vtotal += vbar[i] * moleFracs[i];
    }
    Phase::assignDensity(meanMolecularWeight() / vtotal);
}

void MaskellSolidSolnPhase::setPressure(CanteraDouble p)
{
    m_Pcurrent = p;
}

// Chemical Potentials and Activities

void MaskellSolidSolnPhase::getActivityCoefficients(CanteraDouble* ac) const
{
    static const int cacheId = m_cache.getId();
    CachedArray cached = m_cache.getArray(cacheId);
    if (!cached.validate(temperature(), pressure(), stateMFNumber())) {
        cached.value.resize(2);

        const CanteraDouble r = moleFraction(product_species_index);
        const CanteraDouble pval = p(r);
        const CanteraDouble rfm = r * fm(r);
        const CanteraDouble A = (std::pow(1 - rfm, pval) * std::pow(rfm, pval) * std::pow(r - rfm, 1 - pval)) /
                             (std::pow(1 - r - rfm, 1 + pval) * (1 - r));
        const CanteraDouble B = pval * h_mixing / RT();
        cached.value[product_species_index] = A * std::exp(B);
        cached.value[reactant_species_index] = 1 / (A * r * (1-r) ) * std::exp(-B);
    }
    std::copy(cached.value.begin(), cached.value.end(), ac);
}

void MaskellSolidSolnPhase::getChemPotentials(CanteraDouble* mu) const
{
    const CanteraDouble r = moleFraction(product_species_index);
    const CanteraDouble pval = p(r);
    const CanteraDouble rfm = r * fm(r);
    const CanteraDouble DgbarDr = pval * h_mixing +
                               RT() *
                               std::log( (std::pow(1 - rfm, pval) * std::pow(rfm, pval) * std::pow(r - rfm, 1 - pval) * r) /
                               (std::pow(1 - r - rfm, 1 + pval) * (1 - r)) );
    mu[product_species_index] = RT() * m_g0_RT[product_species_index] + DgbarDr;
    mu[reactant_species_index] = RT() * m_g0_RT[reactant_species_index] - DgbarDr;
}

void MaskellSolidSolnPhase::getChemPotentials_RT(CanteraDouble* mu) const
{
    warn_deprecated("MaskellSolidSolnPhase::getChemPotentials_RT",
                    "To be removed after Cantera 3.0. Use getChemPotentials instead.");
    getChemPotentials(mu);
    for (size_t sp=0; sp < m_kk; ++sp) {
        mu[sp] *= 1.0 / RT();
    }
}

// Partial Molar Properties

void MaskellSolidSolnPhase::getPartialMolarEnthalpies(CanteraDouble* hbar) const
{
    throw NotImplementedError("MaskellSolidSolnPhase::getPartialMolarEnthalpies");
}

void MaskellSolidSolnPhase::getPartialMolarEntropies(CanteraDouble* sbar) const
{
    throw NotImplementedError("MaskellSolidSolnPhase::getPartialMolarEntropies");
}

void MaskellSolidSolnPhase::getPartialMolarCp(CanteraDouble* cpbar) const
{
    throw NotImplementedError("MaskellSolidSolnPhase::getPartialMolarCp");
}

void MaskellSolidSolnPhase::getPartialMolarVolumes(CanteraDouble* vbar) const
{
    getStandardVolumes(vbar);
}

void MaskellSolidSolnPhase::getPureGibbs(CanteraDouble* gpure) const
{
    for (size_t sp=0; sp < m_kk; ++sp) {
        gpure[sp] = RT() * m_g0_RT[sp];
    }
}

void MaskellSolidSolnPhase::getStandardChemPotentials(CanteraDouble* mu) const
{
    // What is the difference between this and getPureGibbs? IdealSolidSolnPhase
    // gives the same for both
    getPureGibbs(mu);
}

// Utility Functions

void MaskellSolidSolnPhase::initThermo()
{
    if (!m_input.empty()) {
        set_h_mix(m_input.convert("excess-enthalpy", "J/kmol"));
        setProductSpecies(m_input["product-species"].asString());
    }
    VPStandardStateTP::initThermo();
}

void MaskellSolidSolnPhase::getParameters(AnyMap& phaseNode) const
{
    VPStandardStateTP::getParameters(phaseNode);
    phaseNode["excess-enthalpy"].setQuantity(h_mixing, "J/kmol");
    phaseNode["product-species"] = speciesName(product_species_index);
}

void MaskellSolidSolnPhase::setProductSpecies(const string& name)
{
    product_species_index = static_cast<int>(speciesIndex(name));
    if (product_species_index == -1) {
        throw CanteraError("MaskellSolidSolnPhase::setProductSpecies",
                           "Species '{}' not found", name);
    }
    reactant_species_index = (product_species_index == 0) ? 1 : 0;
}

CanteraDouble MaskellSolidSolnPhase::s() const
{
    return 1 + std::exp(h_mixing / RT());
}

CanteraDouble MaskellSolidSolnPhase::fm(const CanteraDouble r) const
{
    return (1 - std::sqrt(1 - 4*r*(1-r)/s())) / (2*r);
}

CanteraDouble MaskellSolidSolnPhase::p(const CanteraDouble r) const
{
    const CanteraDouble sval = s();
    return (1 - 2*r) / std::sqrt(sval*sval - 4 * sval * r + 4 * sval * r * r);
}

} // end namespace Cantera
