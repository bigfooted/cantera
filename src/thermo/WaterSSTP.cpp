/**
 *  @file WaterSSTP.cpp
 * Definitions for a ThermoPhase class consisting of pure water (see @ref thermoprops
 * and class @link Cantera::WaterSSTP WaterSSTP@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/WaterSSTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{
WaterSSTP::WaterSSTP(const string& inputFile, const string& id)
{
    initThermoFile(inputFile, id);
}

string WaterSSTP::phaseOfMatter() const {
    const vector<string> phases = {
        "gas", "liquid", "supercritical", "unstable-liquid", "unstable-gas"
    };
    return phases[m_sub.phaseState()];
}

void WaterSSTP::initThermo()
{
    SingleSpeciesTP::initThermo();

    // Calculate the molecular weight. Note while there may be a very good
    // calculated weight in the steam table class, using this weight may lead to
    // codes exhibiting mass loss issues. We need to grab the elemental atomic
    // weights used in the Element class and calculate a consistent H2O
    // molecular weight based on that.
    size_t nH = elementIndex("H");
    if (nH == npos) {
        throw CanteraError("WaterSSTP::initThermo",
                           "H not an element");
    }
    CanteraDouble mw_H = atomicWeight(nH);
    size_t nO = elementIndex("O");
    if (nO == npos) {
        throw CanteraError("WaterSSTP::initThermo",
                           "O not an element");
    }
    CanteraDouble mw_O = atomicWeight(nO);
    m_mw = 2.0 * mw_H + mw_O;
    setMolecularWeight(0,m_mw);

    // Set the baseline
    CanteraDouble T = 298.15;
    Phase::setDensity(7.0E-8);
    Phase::setTemperature(T);

    CanteraDouble presLow = 1.0E-2;
    CanteraDouble oneBar = 1.0E5;
    CanteraDouble dd = m_sub.density(T, presLow, WATER_GAS, 7.0E-8);
    setDensity(dd);
    setTemperature(T);
    SW_Offset = 0.0;
    CanteraDouble s = entropy_mole();
    s -= GasConstant * log(oneBar/presLow);
    if (s != 188.835E3) {
        SW_Offset = 188.835E3 - s;
    }
    s = entropy_mole();
    s -= GasConstant * log(oneBar/presLow);

    CanteraDouble h = enthalpy_mole();
    if (h != -241.826E6) {
        EW_Offset = -241.826E6 - h;
    }
    h = enthalpy_mole();

    // Set the initial state of the system to 298.15 K and 1 bar.
    setTemperature(298.15);
    CanteraDouble rho0 = m_sub.density(298.15, OneAtm, WATER_LIQUID);
    setDensity(rho0);

    m_waterProps = make_unique<WaterProps>(&m_sub);

    // Set the flag to say we are ready to calculate stuff
    m_ready = true;
}

void WaterSSTP::getEnthalpy_RT(CanteraDouble* hrt) const
{
    *hrt = (m_sub.enthalpy_mass() * m_mw + EW_Offset) / RT();
}

void WaterSSTP::getIntEnergy_RT(CanteraDouble* ubar) const
{
    *ubar = (m_sub.intEnergy_mass() * m_mw + EW_Offset)/ RT();
}

void WaterSSTP::getEntropy_R(CanteraDouble* sr) const
{
    sr[0] = (m_sub.entropy_mass() * m_mw + SW_Offset) / GasConstant;
}

void WaterSSTP::getGibbs_RT(CanteraDouble* grt) const
{
    *grt = (m_sub.gibbs_mass() * m_mw + EW_Offset) / RT()
           - SW_Offset / GasConstant;
    if (!m_ready) {
        throw CanteraError("waterSSTP::getGibbs_RT", "Phase not ready");
    }
}

void WaterSSTP::getStandardChemPotentials(CanteraDouble* gss) const
{
    *gss = (m_sub.gibbs_mass() * m_mw + EW_Offset - SW_Offset*temperature());
    if (!m_ready) {
        throw CanteraError("waterSSTP::getStandardChemPotentials",
                           "Phase not ready");
    }
}

void WaterSSTP::getCp_R(CanteraDouble* cpr) const
{
    cpr[0] = m_sub.cp_mass() * m_mw / GasConstant;
}

CanteraDouble WaterSSTP::cv_mole() const
{
    return m_sub.cv_mass() * m_mw;
}

void WaterSSTP::getEnthalpy_RT_ref(CanteraDouble* hrt) const
{
    CanteraDouble p = pressure();
    CanteraDouble T = temperature();
    CanteraDouble dens = density();
    int waterState = WATER_GAS;
    CanteraDouble rc = m_sub.Rhocrit();
    if (dens > rc) {
        waterState = WATER_LIQUID;
    }
    CanteraDouble dd = m_sub.density(T, OneAtm, waterState, dens);
    if (dd <= 0.0) {
        throw CanteraError("WaterSSTP::getEnthalpy_RT_ref", "error");
    }
    CanteraDouble h = m_sub.enthalpy_mass() * m_mw;
    *hrt = (h + EW_Offset) / RT();
    dd = m_sub.density(T, p, waterState, dens);
}

void WaterSSTP::getGibbs_RT_ref(CanteraDouble* grt) const
{
    CanteraDouble p = pressure();
    CanteraDouble T = temperature();
    CanteraDouble dens = density();
    int waterState = WATER_GAS;
    CanteraDouble rc = m_sub.Rhocrit();
    if (dens > rc) {
        waterState = WATER_LIQUID;
    }
    CanteraDouble dd = m_sub.density(T, OneAtm, waterState, dens);
    if (dd <= 0.0) {
        throw CanteraError("WaterSSTP::getGibbs_RT_ref", "error");
    }
    m_sub.setState_TD(T, dd);
    CanteraDouble g = m_sub.gibbs_mass() * m_mw;
    *grt = (g + EW_Offset - SW_Offset*T)/ RT();
    dd = m_sub.density(T, p, waterState, dens);
}

void WaterSSTP::getGibbs_ref(CanteraDouble* g) const
{
    getGibbs_RT_ref(g);
    for (size_t k = 0; k < m_kk; k++) {
        g[k] *= RT();
    }
}

void WaterSSTP::getEntropy_R_ref(CanteraDouble* sr) const
{
    CanteraDouble p = pressure();
    CanteraDouble T = temperature();
    CanteraDouble dens = density();
    int waterState = WATER_GAS;
    CanteraDouble rc = m_sub.Rhocrit();
    if (dens > rc) {
        waterState = WATER_LIQUID;
    }
    CanteraDouble dd = m_sub.density(T, OneAtm, waterState, dens);

    if (dd <= 0.0) {
        throw CanteraError("WaterSSTP::getEntropy_R_ref", "error");
    }
    m_sub.setState_TD(T, dd);

    CanteraDouble s = m_sub.entropy_mass() * m_mw;
    *sr = (s + SW_Offset)/ GasConstant;
    dd = m_sub.density(T, p, waterState, dens);
}

void WaterSSTP::getCp_R_ref(CanteraDouble* cpr) const
{
    CanteraDouble p = pressure();
    CanteraDouble T = temperature();
    CanteraDouble dens = density();
    int waterState = WATER_GAS;
    CanteraDouble rc = m_sub.Rhocrit();
    if (dens > rc) {
        waterState = WATER_LIQUID;
    }
    CanteraDouble dd = m_sub.density(T, OneAtm, waterState, dens);
    m_sub.setState_TD(T, dd);
    if (dd <= 0.0) {
        throw CanteraError("WaterSSTP::getCp_R_ref", "error");
    }
    CanteraDouble cp = m_sub.cp_mass() * m_mw;
    *cpr = cp / GasConstant;
    dd = m_sub.density(T, p, waterState, dens);
}

void WaterSSTP::getStandardVolumes_ref(CanteraDouble* vol) const
{
    CanteraDouble p = pressure();
    CanteraDouble T = temperature();
    CanteraDouble dens = density();
    int waterState = WATER_GAS;
    CanteraDouble rc = m_sub.Rhocrit();
    if (dens > rc) {
        waterState = WATER_LIQUID;
    }
    CanteraDouble dd = m_sub.density(T, OneAtm, waterState, dens);
    if (dd <= 0.0) {
        throw CanteraError("WaterSSTP::getStandardVolumes_ref", "error");
    }
    *vol = meanMolecularWeight() /dd;
    dd = m_sub.density(T, p, waterState, dens);
}

CanteraDouble WaterSSTP::pressure() const
{
    return m_sub.pressure();
}

void WaterSSTP::setPressure(CanteraDouble p)
{
    CanteraDouble T = temperature();
    CanteraDouble dens = density();
    CanteraDouble pp = m_sub.psat(T);
    int waterState = WATER_SUPERCRIT;
    if (T < m_sub.Tcrit()) {
        if (p >= pp) {
            waterState = WATER_LIQUID;
            dens = 1000.;
        } else if (!m_allowGasPhase) {
            throw CanteraError("WaterSSTP::setPressure",
                "Model assumes liquid phase; pressure p = {} lies below\n"
                "the saturation pressure (P_sat = {}).", p, pp);
        }
    }

    CanteraDouble dd = m_sub.density(T, p, waterState, dens);
    if (dd <= 0.0) {
        throw CanteraError("WaterSSTP::setPressure", "Error");
    }
    setDensity(dd);
}

CanteraDouble WaterSSTP::isothermalCompressibility() const
{
    return m_sub.isothermalCompressibility();
}

CanteraDouble WaterSSTP::thermalExpansionCoeff() const
{
    return m_sub.coeffThermExp();
}

CanteraDouble WaterSSTP::dthermalExpansionCoeffdT() const
{
    CanteraDouble pres = pressure();
    CanteraDouble dens_save = density();
    CanteraDouble T = temperature();
    CanteraDouble tt = T - 0.04;
    CanteraDouble dd = m_sub.density(tt, pres, WATER_LIQUID, dens_save);
    if (dd < 0.0) {
        throw CanteraError("WaterSSTP::dthermalExpansionCoeffdT",
            "Unable to solve for the density at T = {}, P = {}", tt, pres);
    }
    CanteraDouble vald = m_sub.coeffThermExp();
    m_sub.setState_TD(T, dens_save);
    CanteraDouble val2 = m_sub.coeffThermExp();
    return (val2 - vald) / 0.04;
}

CanteraDouble WaterSSTP::critTemperature() const
{
    return m_sub.Tcrit();
}

CanteraDouble WaterSSTP::critPressure() const
{
    return m_sub.Pcrit();
}

CanteraDouble WaterSSTP::critDensity() const
{
    return m_sub.Rhocrit();
}

void WaterSSTP::setTemperature(const CanteraDouble temp)
{
    if (temp < 273.16) {
        throw CanteraError("WaterSSTP::setTemperature",
            "Model assumes liquid phase; temperature T = {} lies below\n"
            "the triple point temperature (T_triple = 273.16).", temp);
    }
    Phase::setTemperature(temp);
    m_sub.setState_TD(temp, density());
}

void WaterSSTP::setDensity(const CanteraDouble dens)
{
    Phase::setDensity(dens);
    m_sub.setState_TD(temperature(), dens);
}

CanteraDouble WaterSSTP::satPressure(CanteraDouble t) {
    CanteraDouble tsave = temperature();
    CanteraDouble dsave = density();
    CanteraDouble pp = m_sub.psat(t);
    m_sub.setState_TD(tsave, dsave);
    return pp;
}

CanteraDouble WaterSSTP::vaporFraction() const
{
    if (temperature() >= m_sub.Tcrit()) {
        CanteraDouble dens = density();
        if (dens >= m_sub.Rhocrit()) {
            return 0.0;
        }
        return 1.0;
    }
    // If below tcrit we always return 0 from this class
    return 0.0;
}

}
