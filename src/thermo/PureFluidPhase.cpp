/**
 * @file PureFluidPhase.cpp Definitions for a ThermoPhase object for a pure
 *     fluid phase consisting of gas, liquid, mixed-gas-liquid and supercritical
 *     fluid (see @ref thermoprops and class @link Cantera::PureFluidPhase
 *     PureFluidPhase@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PureFluidPhase.h"

#include "cantera/tpx/Sub.h"
#include "cantera/tpx/utils.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"

namespace Cantera
{

void PureFluidPhase::initThermo()
{
    if (m_input.hasKey("pure-fluid-name")) {
        setSubstance(m_input["pure-fluid-name"].asString());
    }

    m_sub.reset(tpx::newSubstance(m_tpx_name));

    m_mw = m_sub->MolWt();
    setMolecularWeight(0,m_mw);

    CanteraDouble cp0_R, h0_RT, s0_R, p;
    CanteraDouble T0 = 298.15;
    if (T0 < m_sub->Tcrit()) {
        m_sub->Set(tpx::PropertyPair::TX, T0, 1.0);
        p = 0.01*m_sub->P();
    } else {
        p = 0.001*m_sub->Pcrit();
    }
    p = 0.001 * p;
    m_sub->Set(tpx::PropertyPair::TP, T0, p);

    m_spthermo.update_single(0, T0, &cp0_R, &h0_RT, &s0_R);
    CanteraDouble s_R = s0_R - log(p/refPressure());
    m_sub->setStdState(h0_RT*GasConstant*298.15/m_mw,
                       s_R*GasConstant/m_mw, T0, p);
    debuglog("PureFluidPhase::initThermo: initialized phase "
             +name()+"\n", m_verbose);
}

void PureFluidPhase::getParameters(AnyMap& phaseNode) const
{
    ThermoPhase::getParameters(phaseNode);
    phaseNode["pure-fluid-name"] = m_sub->name();
}

vector<string> PureFluidPhase::fullStates() const
{
    return {"TD", "UV", "DP", "HP", "SP", "SV",
            "ST", "TV", "PV", "UP", "VH", "TH", "SH", "TPQ"};
}

vector<string> PureFluidPhase::partialStates() const
{
    return {"TP", "TQ", "PQ"};
}

string PureFluidPhase::phaseOfMatter() const
{
    if (temperature() >= critTemperature() || pressure() >= critPressure()) {
        return "supercritical";
    } else if (m_sub->TwoPhase() == 1) {
        return "liquid-gas-mix";
    } else if (pressure() < m_sub->Ps()) {
        return "gas";
    } else {
        return "liquid";
    }
}

CanteraDouble PureFluidPhase::minTemp(size_t k) const
{
    return m_sub->Tmin();
}

CanteraDouble PureFluidPhase::maxTemp(size_t k) const
{
    return m_sub->Tmax();
}

CanteraDouble PureFluidPhase::enthalpy_mole() const
{
    return m_sub->h() * m_mw;
}

CanteraDouble PureFluidPhase::intEnergy_mole() const
{
    return m_sub->u() * m_mw;
}

CanteraDouble PureFluidPhase::entropy_mole() const
{
    return m_sub->s() * m_mw;
}

CanteraDouble PureFluidPhase::gibbs_mole() const
{
    return m_sub->g() * m_mw;
}

CanteraDouble PureFluidPhase::cp_mole() const
{
    return m_sub->cp() * m_mw;
}

CanteraDouble PureFluidPhase::cv_mole() const
{
    return m_sub->cv() * m_mw;
}

CanteraDouble PureFluidPhase::pressure() const
{
    return m_sub->P();
}

void PureFluidPhase::setPressure(CanteraDouble p)
{
    Set(tpx::PropertyPair::TP, temperature(), p);
    ThermoPhase::setDensity(1.0/m_sub->v());
}

void PureFluidPhase::setTemperature(CanteraDouble T)
{
    ThermoPhase::setTemperature(T);
    Set(tpx::PropertyPair::TV, T, m_sub->v());
}

void PureFluidPhase::setDensity(CanteraDouble rho)
{
    ThermoPhase::setDensity(rho);
    Set(tpx::PropertyPair::TV, m_sub->Temp(), 1.0/rho);
}

void PureFluidPhase::Set(tpx::PropertyPair::type n, CanteraDouble x, CanteraDouble y) const
{
    m_sub->Set(n, x, y);
}

CanteraDouble PureFluidPhase::isothermalCompressibility() const
{
    return m_sub->isothermalCompressibility();
}

CanteraDouble PureFluidPhase::thermalExpansionCoeff() const
{
    return m_sub->thermalExpansionCoeff();
}

tpx::Substance& PureFluidPhase::TPX_Substance()
{
    return *m_sub;
}

void PureFluidPhase::getPartialMolarEnthalpies(CanteraDouble* hbar) const
{
    hbar[0] = enthalpy_mole();
}

void PureFluidPhase::getPartialMolarEntropies(CanteraDouble* sbar) const
{
    sbar[0] = entropy_mole();
}

void PureFluidPhase::getPartialMolarIntEnergies(CanteraDouble* ubar) const
{
    ubar[0] = intEnergy_mole();
}

void PureFluidPhase::getPartialMolarCp(CanteraDouble* cpbar) const
{
    cpbar[0] = cp_mole();
}

void PureFluidPhase::getPartialMolarVolumes(CanteraDouble* vbar) const
{
    vbar[0] = 1.0 / molarDensity();
}

Units PureFluidPhase::standardConcentrationUnits() const
{
    return Units(1.0);
}

void PureFluidPhase::getActivityConcentrations(CanteraDouble* c) const
{
    c[0] = 1.0;
}

CanteraDouble PureFluidPhase::standardConcentration(size_t k) const
{
    return 1.0;
}

void PureFluidPhase::getActivities(CanteraDouble* a) const
{
    a[0] = 1.0;
}

void PureFluidPhase::getStandardChemPotentials(CanteraDouble* mu) const
{
    mu[0] = gibbs_mole();
}

void PureFluidPhase::getEnthalpy_RT(CanteraDouble* hrt) const
{
    hrt[0] = enthalpy_mole() / RT();
}

void PureFluidPhase::getEntropy_R(CanteraDouble* sr) const
{
    sr[0] = entropy_mole() / GasConstant;
}

void PureFluidPhase::getGibbs_RT(CanteraDouble* grt) const
{
    grt[0] = gibbs_mole() / RT();
}

void PureFluidPhase::getEnthalpy_RT_ref(CanteraDouble* hrt) const
{
    CanteraDouble rhoSave = density();
    CanteraDouble t = temperature();
    CanteraDouble plow = 1.0E-8;
    Set(tpx::PropertyPair::TP, t, plow);
    getEnthalpy_RT(hrt);
    Set(tpx::PropertyPair::TV, t, 1 / rhoSave);

}

void PureFluidPhase::getGibbs_RT_ref(CanteraDouble* grt) const
{
    CanteraDouble rhoSave = density();
    CanteraDouble t = temperature();
    CanteraDouble pref = refPressure();
    CanteraDouble plow = 1.0E-8;
    Set(tpx::PropertyPair::TP, t, plow);
    getGibbs_RT(grt);
    grt[0] += log(pref/plow);
    Set(tpx::PropertyPair::TV, t, 1 / rhoSave);
}

void PureFluidPhase::getGibbs_ref(CanteraDouble* g) const
{
    getGibbs_RT_ref(g);
    g[0] *= RT();
}

void PureFluidPhase::getEntropy_R_ref(CanteraDouble* er) const
{
    CanteraDouble rhoSave = density();
    CanteraDouble t = temperature();
    CanteraDouble pref = refPressure();
    CanteraDouble plow = 1.0E-8;
    Set(tpx::PropertyPair::TP, t, plow);
    getEntropy_R(er);
    er[0] -= log(pref/plow);
    Set(tpx::PropertyPair::TV, t, 1 / rhoSave);
}

CanteraDouble PureFluidPhase::critTemperature() const
{
    return m_sub->Tcrit();
}

CanteraDouble PureFluidPhase::critPressure() const
{
    return m_sub->Pcrit();
}

CanteraDouble PureFluidPhase::critDensity() const
{
    return 1.0/m_sub->Vcrit();
}

CanteraDouble PureFluidPhase::satTemperature(CanteraDouble p) const
{
    return m_sub->Tsat(p);
}

/* The next several functions set the state. They run the Substance::Set
 * function, followed by setting the state of the ThermoPhase object
 * to the newly computed temperature and density of the Substance.
 */

void PureFluidPhase::setState_HP(CanteraDouble h, CanteraDouble p, CanteraDouble tol)
{
    Set(tpx::PropertyPair::HP, h, p);
    setState_TD(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_UV(CanteraDouble u, CanteraDouble v, CanteraDouble tol)
{
    Set(tpx::PropertyPair::UV, u, v);
    setState_TD(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_SV(CanteraDouble s, CanteraDouble v, CanteraDouble tol)
{
    Set(tpx::PropertyPair::SV, s, v);
    setState_TD(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_SP(CanteraDouble s, CanteraDouble p, CanteraDouble tol)
{
    Set(tpx::PropertyPair::SP, s, p);
    setState_TD(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_ST(CanteraDouble s, CanteraDouble t, CanteraDouble tol)
{
    Set(tpx::PropertyPair::ST, s, t);
    setState_TD(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_TV(CanteraDouble t, CanteraDouble v, CanteraDouble tol)
{
    Set(tpx::PropertyPair::TV, t, v);
    setState_TD(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_PV(CanteraDouble p, CanteraDouble v, CanteraDouble tol)
{
    Set(tpx::PropertyPair::PV, p, v);
    setState_TD(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_UP(CanteraDouble u, CanteraDouble p, CanteraDouble tol)
{
    Set(tpx::PropertyPair::UP, u, p);
    setState_TD(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_VH(CanteraDouble v, CanteraDouble h, CanteraDouble tol)
{
    Set(tpx::PropertyPair::VH, v, h);
    setState_TD(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_TH(CanteraDouble t, CanteraDouble h, CanteraDouble tol)
{
    Set(tpx::PropertyPair::TH, t, h);
    setState_TD(m_sub->Temp(), 1.0/m_sub->v());
}

void PureFluidPhase::setState_SH(CanteraDouble s, CanteraDouble h, CanteraDouble tol)
{
    Set(tpx::PropertyPair::SH, s, h);
    setState_TD(m_sub->Temp(), 1.0/m_sub->v());
}

CanteraDouble PureFluidPhase::satPressure(CanteraDouble t)
{
    Set(tpx::PropertyPair::TV, t, m_sub->v());
    return m_sub->Ps();
}

CanteraDouble PureFluidPhase::vaporFraction() const
{
    return m_sub->x();
}

void PureFluidPhase::setState_Tsat(CanteraDouble t, CanteraDouble x)
{
    Set(tpx::PropertyPair::TX, t, x);
    ThermoPhase::setTemperature(t);
    ThermoPhase::setDensity(1.0/m_sub->v());
}

void PureFluidPhase::setState_Psat(CanteraDouble p, CanteraDouble x)
{
    Set(tpx::PropertyPair::PX, p, x);
    ThermoPhase::setTemperature(m_sub->Temp());
    ThermoPhase::setDensity(1.0/m_sub->v());
}

string PureFluidPhase::report(bool show_thermo, CanteraDouble threshold) const
{
    fmt::memory_buffer b;
    // This is the width of the first column of names in the report.
    int name_width = 18;

    string blank_leader = fmt::format("{:{}}", "", name_width);

    string one_property = fmt::format("{{:>{}}}   {{:<.5g}} {{}}\n", name_width);

    string two_prop_header = "{}   {:^15}   {:^15}\n";
    string kg_kmol_header = fmt::format(
        two_prop_header, blank_leader, "1 kg", "1 kmol"
    );
    string two_prop_sep = fmt::format(
        "{}   {:-^15}   {:-^15}\n", blank_leader, "", ""
    );
    string two_property = fmt::format(
        "{{:>{}}}   {{:15.5g}}   {{:15.5g}}  {{}}\n", name_width
    );

    if (name() != "") {
        fmt_append(b, "\n  {}:\n", name());
    }
    fmt_append(b, "\n");
    fmt_append(b, one_property, "temperature", temperature(), "K");
    fmt_append(b, one_property, "pressure", pressure(), "Pa");
    fmt_append(b, one_property, "density", density(), "kg/m^3");
    fmt_append(b, one_property,
               "mean mol. weight", meanMolecularWeight(), "kg/kmol");
    fmt_append(b, "{:>{}}   {:<.5g}\n",
               "vapor fraction", name_width, vaporFraction());
    fmt_append(b, "{:>{}}   {}\n",
               "phase of matter", name_width, phaseOfMatter());

    if (show_thermo) {
        fmt_append(b, "\n");
        fmt_append(b, kg_kmol_header);
        fmt_append(b, two_prop_sep);
        fmt_append(b, two_property,
                   "enthalpy", enthalpy_mass(), enthalpy_mole(), "J");
        fmt_append(b, two_property,
                   "internal energy", intEnergy_mass(), intEnergy_mole(), "J");
        fmt_append(b, two_property,
                   "entropy", entropy_mass(), entropy_mole(), "J/K");
        fmt_append(b, two_property,
                   "Gibbs function", gibbs_mass(), gibbs_mole(), "J");
        fmt_append(b, two_property,
                   "heat capacity c_p", cp_mass(), cp_mole(), "J/K");
        fmt_append(b, two_property,
                   "heat capacity c_v", cv_mass(), cv_mole(), "J/K");
    }

    return to_string(b);
}

}
