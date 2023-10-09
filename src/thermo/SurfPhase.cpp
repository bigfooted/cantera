/**
 *  @file SurfPhase.cpp
 *  Definitions for a simple thermodynamic model of a surface phase
 *  derived from ThermoPhase, assuming an ideal solution model
 *  (see @ref thermoprops and class
 *  @link Cantera::SurfPhase SurfPhase@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/EdgePhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

SurfPhase::SurfPhase(const string& infile, const string& id_)
{
    setNDim(2);
    initThermoFile(infile, id_);
}

CanteraDouble SurfPhase::enthalpy_mole() const
{
    if (m_n0 <= 0.0) {
        return 0.0;
    }
    _updateThermo();
    return mean_X(m_h0);
}

CanteraDouble SurfPhase::intEnergy_mole() const
{
    return enthalpy_mole();
}

CanteraDouble SurfPhase::entropy_mole() const
{
    _updateThermo();
    CanteraDouble s = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        s += moleFraction(k) * (m_s0[k] -
            GasConstant * log(std::max(concentration(k) * size(k)/m_n0, SmallNumber)));
    }
    return s;
}

CanteraDouble SurfPhase::cp_mole() const
{
    _updateThermo();
    return mean_X(m_cp0);
}

CanteraDouble SurfPhase::cv_mole() const
{
    return cp_mole();
}

void SurfPhase::getPartialMolarEnthalpies(CanteraDouble* hbar) const
{
    getEnthalpy_RT(hbar);
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT();
    }
}

void SurfPhase::getPartialMolarEntropies(CanteraDouble* sbar) const
{
    getEntropy_R(sbar);
    getActivityConcentrations(m_work.data());
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] -= log(std::max(m_work[k], SmallNumber)) - logStandardConc(k);
        sbar[k] *= GasConstant;
    }
}

void SurfPhase::getPartialMolarCp(CanteraDouble* cpbar) const
{
    getCp_R(cpbar);
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
}

// HKM 9/1/11  The partial molar volumes returned here are really partial molar areas.
//             Partial molar volumes for this phase should actually be equal to zero.
void SurfPhase::getPartialMolarVolumes(CanteraDouble* vbar) const
{
    getStandardVolumes(vbar);
}

void SurfPhase::getStandardChemPotentials(CanteraDouble* mu0) const
{
    _updateThermo();
    copy(m_mu0.begin(), m_mu0.end(), mu0);
}

void SurfPhase::getChemPotentials(CanteraDouble* mu) const
{
    _updateThermo();
    copy(m_mu0.begin(), m_mu0.end(), mu);
    getActivityConcentrations(m_work.data());
    for (size_t k = 0; k < m_kk; k++) {
        mu[k] += RT() * (log(std::max(m_work[k], SmallNumber)) - logStandardConc(k));
    }
}

void SurfPhase::getActivityConcentrations(CanteraDouble* c) const
{
    getConcentrations(c);
}

CanteraDouble SurfPhase::standardConcentration(size_t k) const
{
    return m_n0/size(k);
}

CanteraDouble SurfPhase::logStandardConc(size_t k) const
{
    return m_logn0 - m_logsize[k];
}

void SurfPhase::getPureGibbs(CanteraDouble* g) const
{
    _updateThermo();
    copy(m_mu0.begin(), m_mu0.end(), g);
}

void SurfPhase::getGibbs_RT(CanteraDouble* grt) const
{
    _updateThermo();
    scale(m_mu0.begin(), m_mu0.end(), grt, 1.0/RT());
}

void SurfPhase::getEnthalpy_RT(CanteraDouble* hrt) const
{
    _updateThermo();
    scale(m_h0.begin(), m_h0.end(), hrt, 1.0/RT());
}

void SurfPhase::getEntropy_R(CanteraDouble* sr) const
{
    _updateThermo();
    scale(m_s0.begin(), m_s0.end(), sr, 1.0/GasConstant);
}

void SurfPhase::getCp_R(CanteraDouble* cpr) const
{
    _updateThermo();
    scale(m_cp0.begin(), m_cp0.end(), cpr, 1.0/GasConstant);
}

void SurfPhase::getStandardVolumes(CanteraDouble* vol) const
{
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = 0.0;
    }
}

void SurfPhase::getGibbs_RT_ref(CanteraDouble* grt) const
{
    getGibbs_RT(grt);
}

void SurfPhase::getEnthalpy_RT_ref(CanteraDouble* hrt) const
{
    getEnthalpy_RT(hrt);
}

void SurfPhase::getEntropy_R_ref(CanteraDouble* sr) const
{
    getEntropy_R(sr);
}

void SurfPhase::getCp_R_ref(CanteraDouble* cprt) const
{
    getCp_R(cprt);
}

bool SurfPhase::addSpecies(shared_ptr<Species> spec)
{
    bool added = ThermoPhase::addSpecies(spec);
    if (added) {
        m_h0.push_back(0.0);
        m_s0.push_back(0.0);
        m_cp0.push_back(0.0);
        m_mu0.push_back(0.0);
        m_work.push_back(0.0);
        m_speciesSize.push_back(spec->size);
        m_logsize.push_back(log(spec->size));
        if (m_kk == 1) {
            vector<CanteraDouble> cov{1.0};
            setCoverages(cov.data());
        }
    }
    return added;
}

void SurfPhase::setMolarDensity(const CanteraDouble vm) {
    warn_deprecated("SurfPhase::setMolarDensity", "To be removed after Cantera 3.0");
    if (vm != 0.0) {
        throw CanteraError("SurfPhase::setMolarDensity",
                            "The volume of an interface is zero");
    }
}

void SurfPhase::setSiteDensity(CanteraDouble n0)
{
    if (n0 <= 0.0) {
        throw CanteraError("SurfPhase::setSiteDensity",
                           "Site density must be positive. Got {}", n0);
    }
    m_n0 = n0;
    assignDensity(n0 * meanMolecularWeight());
    m_logn0 = log(m_n0);
}

void SurfPhase::setCoverages(const CanteraDouble* theta)
{
    CanteraDouble sum = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        sum += theta[k] / size(k);
    }
    if (sum <= 0.0) {
        throw CanteraError("SurfPhase::setCoverages",
                           "Sum of Coverage fractions is zero or negative");
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_work[k] = theta[k] / (sum * size(k));
    }
    setMoleFractions(m_work.data());
}

void SurfPhase::setCoveragesNoNorm(const CanteraDouble* theta)
{
    CanteraDouble sum = 0.0;
    CanteraDouble sum2 = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        sum += theta[k] / size(k);
        sum2 += theta[k];
    }
    if (sum <= 0.0) {
        throw CanteraError("SurfPhase::setCoverages",
                           "Sum of Coverage fractions is zero or negative");
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_work[k] = theta[k] * sum2 / (sum * size(k));
    }
    setMoleFractions_NoNorm(m_work.data());
}

void SurfPhase::getCoverages(CanteraDouble* theta) const
{
    CanteraDouble sum_X = 0.0;
    CanteraDouble sum_X_s = 0.0;
    getMoleFractions(theta);
    for (size_t k = 0; k < m_kk; k++) {
        sum_X += theta[k];
        sum_X_s += theta[k] * size(k);
    }
    for (size_t k = 0; k < m_kk; k++) {
        theta[k] *= size(k) * sum_X / sum_X_s;
    }
}

void SurfPhase::setCoveragesByName(const string& cov)
{
    setCoveragesByName(parseCompString(cov, speciesNames()));
}

void SurfPhase::setCoveragesByName(const Composition& cov)
{
    vector<CanteraDouble> cv(m_kk, 0.0);
    bool ifound = false;
    for (size_t k = 0; k < m_kk; k++) {
        CanteraDouble c = getValue(cov, speciesName(k), 0.0);
        if (c > 0.0) {
            ifound = true;
            cv[k] = c;
        }
    }
    if (!ifound) {
        throw CanteraError("SurfPhase::setCoveragesByName",
                           "Input coverages are all zero or negative");
    }
    setCoverages(cv.data());
}

void SurfPhase::setState(const AnyMap& state) {
    if (state.hasKey("coverages")) {
        if (state["coverages"].is<string>()) {
            setCoveragesByName(state["coverages"].asString());
        } else {
            setCoveragesByName(state["coverages"].asMap<CanteraDouble>());
        }
    }
    ThermoPhase::setState(state);
}

void SurfPhase::compositionChanged()
{
    ThermoPhase::compositionChanged();
    assignDensity(m_n0 * meanMolecularWeight());
}

void SurfPhase::_updateThermo(bool force) const
{
    CanteraDouble tnow = temperature();
    if (m_tlast != tnow || force) {
        m_spthermo.update(tnow, m_cp0.data(), m_h0.data(), m_s0.data());
        m_tlast = tnow;
        for (size_t k = 0; k < m_kk; k++) {
            m_h0[k] *= GasConstant * tnow;
            m_s0[k] *= GasConstant;
            m_cp0[k] *= GasConstant;
            m_mu0[k] = m_h0[k] - tnow*m_s0[k];
        }
        m_tlast = tnow;
    }
}

void SurfPhase::initThermo()
{
    if (m_input.hasKey("site-density")) {
        // Units are kmol/m^2 for surface phases or kmol/m for edge phases
        setSiteDensity(m_input.convert("site-density",
            Units(1.0, 0, -static_cast<CanteraDouble>(m_ndim), 0, 0, 0, 1)));
    }
}

void SurfPhase::getParameters(AnyMap& phaseNode) const
{
    ThermoPhase::getParameters(phaseNode);
    phaseNode["site-density"].setQuantity(
        m_n0, Units(1.0, 0, -static_cast<CanteraDouble>(m_ndim), 0, 0, 0, 1));
}

EdgePhase::EdgePhase(const string& infile, const string& id_)
{
    setNDim(1);
    initThermoFile(infile, id_);
}

}
