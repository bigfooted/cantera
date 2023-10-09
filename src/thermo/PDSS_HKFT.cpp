/**
 *  @file PDSS_HKFT.cpp
 *    Definitions for the class PDSS_HKFT (pressure dependent standard state)
 *    which handles calculations for a single species in a phase using the
 *    HKFT standard state
 *    (see @ref pdssthermo and class @link Cantera::PDSS_HKFT PDSS_HKFT@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PDSS_HKFT.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/VPStandardStateTP.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"

namespace Cantera
{
// Set the default to error exit if there is an input file inconsistency
int PDSS_HKFT::s_InputInconsistencyErrorExit = 1;

PDSS_HKFT::PDSS_HKFT()
{
    m_pres = OneAtm;
    m_units.setDefaults({"cm", "gmol", "cal", "bar"});
}

CanteraDouble PDSS_HKFT::enthalpy_mole() const
{
    // Ok we may change this evaluation method in the future.
    CanteraDouble h = gibbs_mole() + m_temp * entropy_mole();
    return h;
}

CanteraDouble PDSS_HKFT::enthalpy_mole2() const
{
    warn_deprecated("PDSS_HKFT::enthalpy_mole2", "To be removed after Cantera 3.0");
    CanteraDouble enthTRPR = m_Mu0_tr_pr + 298.15 * m_units.convertTo(m_Entrop_tr_pr, "J/kmol");
    return deltaH() + enthTRPR;
}

CanteraDouble PDSS_HKFT::intEnergy_mole() const
{
    return enthalpy_RT() - molarVolume() * m_pres;
}

CanteraDouble PDSS_HKFT::entropy_mole() const
{
    return m_units.convertTo(m_Entrop_tr_pr, "J/kmol") + deltaS();
}

CanteraDouble PDSS_HKFT::gibbs_mole() const
{
    return m_Mu0_tr_pr + deltaG();
}

CanteraDouble PDSS_HKFT::cp_mole() const
{
    CanteraDouble pbar = m_pres * 1.0E-5;
    CanteraDouble c1term = m_c1;
    CanteraDouble c2term = m_c2 / (m_temp - 228.) / (m_temp - 228.);
    CanteraDouble a3term = -m_a3 / (m_temp - 228.) / (m_temp - 228.) / (m_temp - 228.) * 2.0 * m_temp * (pbar - m_presR_bar);
    CanteraDouble a4term = -m_a4 / (m_temp - 228.) / (m_temp - 228.) / (m_temp - 228.) * 2.0 * m_temp
                        * log((2600. + pbar)/(2600. + m_presR_bar));

    CanteraDouble omega_j;
    CanteraDouble domega_jdT;
    CanteraDouble d2omega_jdT2;
    if (m_charge_j == 0.0) {
        omega_j = m_omega_pr_tr;
        domega_jdT = 0.0;
        d2omega_jdT2 = 0.0;
    } else {
        CanteraDouble nu = 166027;
        CanteraDouble r_e_j_pr_tr = m_charge_j * m_charge_j / (m_omega_pr_tr/nu + m_charge_j/3.082);
        CanteraDouble gval = gstar(m_temp, m_pres, 0);
        CanteraDouble dgvaldT = gstar(m_temp, m_pres, 1);
        CanteraDouble d2gvaldT2 = gstar(m_temp, m_pres, 2);

        CanteraDouble r_e_j = r_e_j_pr_tr + fabs(m_charge_j) * gval;
        CanteraDouble dr_e_jdT = fabs(m_charge_j) * dgvaldT;
        CanteraDouble d2r_e_jdT2 = fabs(m_charge_j) * d2gvaldT2;
        CanteraDouble r_e_j2 = r_e_j * r_e_j;

        CanteraDouble charge2 = m_charge_j * m_charge_j;
        CanteraDouble r_e_H = 3.082 + gval;
        CanteraDouble r_e_H2 = r_e_H * r_e_H;
        omega_j = nu * (charge2 / r_e_j - m_charge_j / r_e_H);
        domega_jdT = nu * (-(charge2 / r_e_j2 * dr_e_jdT)
                            +(m_charge_j / r_e_H2 * dgvaldT));
        d2omega_jdT2 = nu * (2.0*charge2*dr_e_jdT*dr_e_jdT/(r_e_j2*r_e_j) - charge2*d2r_e_jdT2/r_e_j2
                             -2.0*m_charge_j*dgvaldT*dgvaldT/(r_e_H2*r_e_H) + m_charge_j*d2gvaldT2 /r_e_H2);
    }

    CanteraDouble relepsilon = m_waterProps->relEpsilon(m_temp, m_pres, 0);
    CanteraDouble drelepsilondT = m_waterProps->relEpsilon(m_temp, m_pres, 1);
    CanteraDouble Y = drelepsilondT / (relepsilon * relepsilon);
    CanteraDouble d2relepsilondT2 = m_waterProps->relEpsilon(m_temp, m_pres, 2);

    CanteraDouble X = d2relepsilondT2 / (relepsilon* relepsilon) - 2.0 * relepsilon * Y * Y;
    CanteraDouble Z = -1.0 / relepsilon;

    CanteraDouble yterm = 2.0 * m_temp * Y * domega_jdT;
    CanteraDouble xterm = omega_j * m_temp * X;
    CanteraDouble otterm = m_temp * d2omega_jdT2 * (Z + 1.0);
    CanteraDouble rterm = - m_domega_jdT_prtr * (m_Z_pr_tr + 1.0);

    CanteraDouble Cp_calgmol = c1term + c2term + a3term + a4term + yterm + xterm + otterm + rterm;

    return m_units.convertTo(Cp_calgmol, "J/kmol");
}

CanteraDouble PDSS_HKFT::molarVolume() const
{
    // Initially do all calculations in (cal/gmol/Pa)

    CanteraDouble a1term = m_a1 * 1.0E-5;
    CanteraDouble a2term = m_a2 / (2600.E5 + m_pres);
    CanteraDouble a3term = m_a3 * 1.0E-5/ (m_temp - 228.);
    CanteraDouble a4term = m_a4 / (m_temp - 228.) / (2600.E5 + m_pres);

    CanteraDouble omega_j;
    CanteraDouble domega_jdP;
    if (m_charge_j == 0.0) {
        omega_j = m_omega_pr_tr;
        domega_jdP = 0.0;
    } else {
        CanteraDouble nu = 166027.;
        CanteraDouble charge2 = m_charge_j * m_charge_j;
        CanteraDouble r_e_j_pr_tr = charge2 / (m_omega_pr_tr/nu + m_charge_j/3.082);

        CanteraDouble gval = gstar(m_temp, m_pres, 0);
        CanteraDouble dgvaldP = gstar(m_temp, m_pres, 3);

        CanteraDouble r_e_j = r_e_j_pr_tr + fabs(m_charge_j) * gval;
        CanteraDouble r_e_H = 3.082 + gval;

        omega_j = nu * (charge2 / r_e_j - m_charge_j / r_e_H);
        CanteraDouble dr_e_jdP = fabs(m_charge_j) * dgvaldP;
        domega_jdP = - nu * (charge2 / (r_e_j * r_e_j) * dr_e_jdP)
                     + nu * m_charge_j / (r_e_H * r_e_H) * dgvaldP;
    }

    CanteraDouble drelepsilondP = m_waterProps->relEpsilon(m_temp, m_pres, 3);
    CanteraDouble relepsilon = m_waterProps->relEpsilon(m_temp, m_pres, 0);
    CanteraDouble Q = drelepsilondP / (relepsilon * relepsilon);
    CanteraDouble Z = -1.0 / relepsilon;
    CanteraDouble wterm = - domega_jdP * (Z + 1.0);
    CanteraDouble qterm = - omega_j * Q;
    CanteraDouble molVol_calgmolPascal = a1term + a2term + a3term + a4term + wterm + qterm;

    // Convert to m**3 / kmol from (cal/gmol/Pa)
    return m_units.convertTo(molVol_calgmolPascal, "J/kmol");
}

CanteraDouble PDSS_HKFT::density() const
{
    return m_mw / molarVolume();
}

CanteraDouble PDSS_HKFT::gibbs_RT_ref() const
{
    CanteraDouble m_psave = m_pres;
    m_pres = m_waterSS->pref_safe(m_temp);
    CanteraDouble ee = gibbs_RT();
    m_pres = m_psave;
    return ee;
}

CanteraDouble PDSS_HKFT::enthalpy_RT_ref() const
{
    CanteraDouble m_psave = m_pres;
    m_pres = m_waterSS->pref_safe(m_temp);
    CanteraDouble hh = enthalpy_RT();
    m_pres = m_psave;
    return hh;
}

CanteraDouble PDSS_HKFT::entropy_R_ref() const
{
    CanteraDouble m_psave = m_pres;
    m_pres = m_waterSS->pref_safe(m_temp);
    CanteraDouble ee = entropy_R();
    m_pres = m_psave;
    return ee;
}

CanteraDouble PDSS_HKFT::cp_R_ref() const
{
    CanteraDouble m_psave = m_pres;
    m_pres = m_waterSS->pref_safe(m_temp);
    CanteraDouble ee = cp_R();
    m_pres = m_psave;
    return ee;
}

CanteraDouble PDSS_HKFT::molarVolume_ref() const
{
    CanteraDouble m_psave = m_pres;
    m_pres = m_waterSS->pref_safe(m_temp);
    CanteraDouble ee = molarVolume();
    m_pres = m_psave;
    return ee;
}

void PDSS_HKFT::setState_TP(CanteraDouble temp, CanteraDouble pres)
{
    setTemperature(temp);
    setPressure(pres);
}

void PDSS_HKFT::initThermo()
{
    PDSS::initThermo();
    if (m_input.hasKey("h0")) {
        m_deltaH_formation_tr_pr = m_input.convert("h0", "cal/gmol");
    }
    if (m_input.hasKey("g0")) {
        m_deltaG_formation_tr_pr = m_input.convert("g0", "cal/gmol");
    }
    if (m_input.hasKey("s0")) {
        m_Entrop_tr_pr = m_input.convert("s0", "cal/gmol/K");
    }
    auto& units = m_input.units();
    if (m_input.hasKey("a")) {
        auto& a = m_input["a"].asVector<AnyValue>(4);
        m_a1 = units.convert(a[0], "cal/gmol/bar");
        m_a2 = units.convert(a[1], "cal/gmol");
        m_a3 = units.convert(a[2], "cal*K/gmol/bar");
        m_a4 = units.convert(a[3], "cal*K/gmol");
    }
    if (m_input.hasKey("c")) {
        auto& c = m_input["c"].asVector<AnyValue>(2);
        m_c1 = units.convert(c[0], "cal/gmol/K");
        m_c2 = units.convert(c[1], "cal*K/gmol");
    }
    if (m_input.hasKey("omega")) {
        m_omega_pr_tr = m_input.convert("omega", "cal/gmol");
    }

    // Ok, if we are missing one, then we construct its value from the other two.
    // This code has been internally verified.
    m_charge_j = m_tp->charge(m_spindex);
    if (std::isnan(m_deltaH_formation_tr_pr)) {
        convertDGFormation();
        CanteraDouble Hcalc = m_Mu0_tr_pr + 298.15 * m_units.convertTo(m_Entrop_tr_pr, "J/kmol");
        m_deltaH_formation_tr_pr = m_units.convertFrom(Hcalc, "J/kmol");
    } else if (std::isnan(m_deltaG_formation_tr_pr)) {
        CanteraDouble DHjmol = m_units.convertTo(m_deltaH_formation_tr_pr, "J/kmol");
        m_Mu0_tr_pr = DHjmol - 298.15 * m_units.convertTo(m_Entrop_tr_pr, "J/kmol");
        m_deltaG_formation_tr_pr = m_units.convertFrom(m_Mu0_tr_pr, "J/kmol");
        CanteraDouble tmp = m_Mu0_tr_pr;
        convertDGFormation();
        CanteraDouble totalSum = m_Mu0_tr_pr - tmp;
        m_Mu0_tr_pr = tmp;
        m_deltaG_formation_tr_pr = m_units.convertFrom(m_Mu0_tr_pr - totalSum, "J/kmol");
    } else if (std::isnan(m_Entrop_tr_pr)) {
        convertDGFormation();
        CanteraDouble DHjmol = m_units.convertTo(m_deltaH_formation_tr_pr, "J/kmol");
        m_Entrop_tr_pr = m_units.convertFrom((DHjmol - m_Mu0_tr_pr) / 298.15, "J/kmol");
    }

    m_waterSS = &dynamic_cast<PDSS_Water&>(*m_tp->providePDSS(0));

    // Section to initialize m_Z_pr_tr and m_Y_pr_tr
    m_temp = 273.15 + 25.;
    m_pres = OneAtm;
    CanteraDouble relepsilon = m_waterProps->relEpsilon(m_temp, m_pres, 0);
    m_waterSS->setState_TP(m_temp, m_pres);
    m_densWaterSS = m_waterSS->density();
    m_Z_pr_tr = -1.0 / relepsilon;
    CanteraDouble drelepsilondT = m_waterProps->relEpsilon(m_temp, m_pres, 1);
    m_Y_pr_tr = drelepsilondT / (relepsilon * relepsilon);
    m_waterProps = make_unique<WaterProps>(m_waterSS);
    m_presR_bar = OneAtm / 1.0E5;
    m_presR_bar = 1.0;
    convertDGFormation();

    // Ok, we have mu. Let's check it against the input value
    // of DH_F to see that we have some internal consistency
    CanteraDouble Hcalc = m_Mu0_tr_pr + 298.15 * m_units.convertTo(m_Entrop_tr_pr, "J/kmol");
    CanteraDouble DHjmol = m_units.convertTo(m_deltaH_formation_tr_pr, "J/kmol");

    // If the discrepancy is greater than 100 cal gmol-1, print
    // an error and exit.
    if (fabs(Hcalc -DHjmol) > m_units.convertTo(100, "J/kmol")) {
        string sname = m_tp->speciesName(m_spindex);
        if (s_InputInconsistencyErrorExit) {
            throw CanteraError("PDSS_HKFT::initThermo", "For {}, DHjmol is"
                " not consistent with G and S: {} vs {} J/kmol",
                sname, Hcalc, m_units.convertTo(m_deltaH_formation_tr_pr, "J/kmol"));
        } else {
            warn_user("PDSS_HKFT::initThermo",
                "DHjmol for {} is not consistent with G and S: calculated {} "
                "vs input {} J/kmol; continuing with consistent DHjmol = {}",
                sname, Hcalc, m_units.convertTo(m_deltaH_formation_tr_pr, "J/kmol"),
                Hcalc);
            m_deltaH_formation_tr_pr = m_units.convertFrom(Hcalc, "J/kmol");
        }
    }
    CanteraDouble nu = 166027;
    CanteraDouble r_e_j_pr_tr;
    if (m_charge_j != 0.0) {
        r_e_j_pr_tr = m_charge_j * m_charge_j / (m_omega_pr_tr/nu + m_charge_j/3.082);
    } else {
        r_e_j_pr_tr = 0.0;
    }

    if (m_charge_j == 0.0) {
        m_domega_jdT_prtr = 0.0;
    } else {
        CanteraDouble gval = gstar(m_temp, m_pres, 0);
        CanteraDouble dgvaldT = gstar(m_temp, m_pres, 1);
        CanteraDouble r_e_j = r_e_j_pr_tr + fabs(m_charge_j) * gval;
        CanteraDouble dr_e_jdT = fabs(m_charge_j) * dgvaldT;
        m_domega_jdT_prtr = - nu * (m_charge_j * m_charge_j / (r_e_j * r_e_j) * dr_e_jdT)
                             + nu * m_charge_j / (3.082 + gval) / (3.082 + gval) * dgvaldT;
    }
}

void PDSS_HKFT::setDeltaH0(CanteraDouble dh0) {
    m_deltaH_formation_tr_pr = m_units.convertFrom(dh0, "J/kmol");
}

void PDSS_HKFT::setDeltaG0(CanteraDouble dg0) {
    m_deltaG_formation_tr_pr = m_units.convertFrom(dg0, "J/kmol");
}

void PDSS_HKFT::setS0(CanteraDouble s0) {
    m_Entrop_tr_pr = m_units.convertFrom(s0, "J/kmol/K");
}

void PDSS_HKFT::set_a(CanteraDouble* a) {
    m_a1 = m_units.convertFrom(a[0], "J/kmol/Pa");
    m_a2 = m_units.convertFrom(a[1], "J/kmol");
    m_a3 = m_units.convertFrom(a[2], "J*K/kmol/Pa");
    m_a4 = m_units.convertFrom(a[3], "J*K/kmol");
}

void PDSS_HKFT::set_c(CanteraDouble* c) {
    m_c1 = m_units.convertFrom(c[0], "J/kmol/K");
    m_c2 = m_units.convertFrom(c[1], "J*K/kmol");
}

void PDSS_HKFT::setOmega(CanteraDouble omega) {
    m_omega_pr_tr = m_units.convertFrom(omega, "J/kmol");
}

void PDSS_HKFT::getParameters(AnyMap& eosNode) const
{
    PDSS::getParameters(eosNode);
    eosNode["model"] = "HKFT";
    eosNode["h0"].setQuantity(m_deltaH_formation_tr_pr, "cal/gmol");
    eosNode["g0"].setQuantity(m_deltaG_formation_tr_pr, "cal/gmol");
    eosNode["s0"].setQuantity(m_Entrop_tr_pr, "cal/gmol/K");

    vector<AnyValue> a(4), c(2);
    a[0].setQuantity(m_a1, "cal/gmol/bar");
    a[1].setQuantity(m_a2, "cal/gmol");
    a[2].setQuantity(m_a3, "cal*K/gmol/bar");
    a[3].setQuantity(m_a4, "cal*K/gmol");
    eosNode["a"] = std::move(a);

    c[0].setQuantity(m_c1, "cal/gmol/K");
    c[1].setQuantity(m_c2, "cal*K/gmol");
    eosNode["c"] = std::move(c);

    eosNode["omega"].setQuantity(m_omega_pr_tr, "cal/gmol");
}

CanteraDouble PDSS_HKFT::deltaH() const
{
    warn_deprecated("PDSS_HKFT::deltaH", "To be removed after Cantera 3.0");
    CanteraDouble pbar = m_pres * 1.0E-5;
    CanteraDouble c1term = m_c1 * (m_temp - 298.15);
    CanteraDouble a1term = m_a1 * (pbar - m_presR_bar);
    CanteraDouble a2term = m_a2 * log((2600. + pbar)/(2600. + m_presR_bar));
    CanteraDouble c2term = -m_c2 * (1.0/(m_temp - 228.) - 1.0/(298.15 - 228.));
    CanteraDouble a3tmp = (2.0 * m_temp - 228.)/ (m_temp - 228.) /(m_temp - 228.);
    CanteraDouble a3term = m_a3 * a3tmp * (pbar - m_presR_bar);
    CanteraDouble a4term = m_a4 * a3tmp * log((2600. + pbar)/(2600. + m_presR_bar));
    CanteraDouble omega_j;
    CanteraDouble domega_jdT;
    if (m_charge_j == 0.0) {
        omega_j = m_omega_pr_tr;
        domega_jdT = 0.0;
    } else {
        CanteraDouble nu = 166027;
        CanteraDouble r_e_j_pr_tr = m_charge_j * m_charge_j / (m_omega_pr_tr/nu + m_charge_j/3.082);
        CanteraDouble gval = gstar(m_temp, m_pres, 0);
        CanteraDouble r_e_j = r_e_j_pr_tr + fabs(m_charge_j) * gval;
        CanteraDouble dgvaldT = gstar(m_temp, m_pres, 1);
        CanteraDouble dr_e_jdT = fabs(m_charge_j) * dgvaldT;
        omega_j = nu * (m_charge_j * m_charge_j / r_e_j - m_charge_j / (3.082 + gval));
        domega_jdT = - nu * (m_charge_j * m_charge_j / (r_e_j * r_e_j) * dr_e_jdT)
                     + nu * m_charge_j / (3.082 + gval) / (3.082 + gval) * dgvaldT;
    }

    CanteraDouble relepsilon = m_waterProps->relEpsilon(m_temp, m_pres, 0);
    CanteraDouble drelepsilondT = m_waterProps->relEpsilon(m_temp, m_pres, 1);

    CanteraDouble Y = drelepsilondT / (relepsilon * relepsilon);
    CanteraDouble Z = -1.0 / relepsilon;

    CanteraDouble yterm = m_temp * omega_j * Y;
    CanteraDouble yrterm = - 298.15 * m_omega_pr_tr * m_Y_pr_tr;

    CanteraDouble wterm = - omega_j * (Z + 1.0);
    CanteraDouble wrterm = + m_omega_pr_tr * (m_Z_pr_tr + 1.0);

    CanteraDouble otterm = m_temp * domega_jdT * (Z + 1.0);
    CanteraDouble otrterm = - m_temp * m_domega_jdT_prtr * (m_Z_pr_tr + 1.0);

    CanteraDouble deltaH_calgmol = c1term + a1term + a2term + c2term + a3term + a4term
                                + yterm + yrterm + wterm + wrterm + otterm + otrterm;

    return m_units.convertTo(deltaH_calgmol, "J/kmol");
}

CanteraDouble PDSS_HKFT::deltaG() const
{
    CanteraDouble pbar = m_pres * 1.0E-5;
    CanteraDouble sterm = - m_Entrop_tr_pr * (m_temp - 298.15);
    CanteraDouble c1term = -m_c1 * (m_temp * log(m_temp/298.15) - (m_temp - 298.15));
    CanteraDouble a1term = m_a1 * (pbar - m_presR_bar);
    CanteraDouble a2term = m_a2 * log((2600. + pbar)/(2600. + m_presR_bar));
    CanteraDouble c2term = -m_c2 * ((1.0/(m_temp - 228.) - 1.0/(298.15 - 228.)) * (228. - m_temp)/228.
                                 - m_temp / (228.*228.) * log((298.15*(m_temp-228.)) / (m_temp*(298.15-228.))));
    CanteraDouble a3term = m_a3 / (m_temp - 228.) * (pbar - m_presR_bar);
    CanteraDouble a4term = m_a4 / (m_temp - 228.) * log((2600. + pbar)/(2600. + m_presR_bar));

    CanteraDouble omega_j;
    if (m_charge_j == 0.0) {
        omega_j = m_omega_pr_tr;
    } else {
        CanteraDouble nu = 166027;
        CanteraDouble r_e_j_pr_tr = m_charge_j * m_charge_j / (m_omega_pr_tr/nu + m_charge_j/3.082);
        CanteraDouble gval = gstar(m_temp, m_pres, 0);
        CanteraDouble r_e_j = r_e_j_pr_tr + fabs(m_charge_j) * gval;
        omega_j = nu * (m_charge_j * m_charge_j / r_e_j - m_charge_j / (3.082 + gval));
    }

    CanteraDouble relepsilon = m_waterProps->relEpsilon(m_temp, m_pres, 0);
    CanteraDouble Z = -1.0 / relepsilon;
    CanteraDouble wterm = - omega_j * (Z + 1.0);
    CanteraDouble wrterm = m_omega_pr_tr * (m_Z_pr_tr + 1.0);
    CanteraDouble yterm = m_omega_pr_tr * m_Y_pr_tr * (m_temp - 298.15);
    CanteraDouble deltaG_calgmol = sterm + c1term + a1term + a2term + c2term + a3term + a4term + wterm + wrterm + yterm;

    return m_units.convertTo(deltaG_calgmol, "J/kmol");
}

CanteraDouble PDSS_HKFT::deltaS() const
{
    CanteraDouble pbar = m_pres * 1.0E-5;

    CanteraDouble c1term = m_c1 * log(m_temp/298.15);
    CanteraDouble c2term = -m_c2 / 228. * ((1.0/(m_temp - 228.) - 1.0/(298.15 - 228.))
                                        + 1.0 / 228. * log((298.15*(m_temp-228.)) / (m_temp*(298.15-228.))));
    CanteraDouble a3term = m_a3 / (m_temp - 228.) / (m_temp - 228.) * (pbar - m_presR_bar);
    CanteraDouble a4term = m_a4 / (m_temp - 228.) / (m_temp - 228.) * log((2600. + pbar)/(2600. + m_presR_bar));

    CanteraDouble omega_j;
    CanteraDouble domega_jdT;
    if (m_charge_j == 0.0) {
        omega_j = m_omega_pr_tr;
        domega_jdT = 0.0;
    } else {
        CanteraDouble nu = 166027;
        CanteraDouble r_e_j_pr_tr = m_charge_j * m_charge_j / (m_omega_pr_tr/nu + m_charge_j/3.082);
        CanteraDouble gval = gstar(m_temp, m_pres, 0);
        CanteraDouble dgvaldT = gstar(m_temp, m_pres, 1);
        CanteraDouble r_e_j = r_e_j_pr_tr + fabs(m_charge_j) * gval;
        CanteraDouble dr_e_jdT = fabs(m_charge_j) * dgvaldT;
        omega_j = nu * (m_charge_j * m_charge_j / r_e_j - m_charge_j / (3.082 + gval));
        domega_jdT = - nu * (m_charge_j * m_charge_j / (r_e_j * r_e_j) * dr_e_jdT)
                     + nu * m_charge_j / (3.082 + gval) / (3.082 + gval) * dgvaldT;
    }

    CanteraDouble relepsilon = m_waterProps->relEpsilon(m_temp, m_pres, 0);
    CanteraDouble drelepsilondT = m_waterProps->relEpsilon(m_temp, m_pres, 1);
    CanteraDouble Y = drelepsilondT / (relepsilon * relepsilon);
    CanteraDouble Z = -1.0 / relepsilon;
    CanteraDouble wterm = omega_j * Y;
    CanteraDouble wrterm = - m_omega_pr_tr * m_Y_pr_tr;
    CanteraDouble otterm = domega_jdT * (Z + 1.0);
    CanteraDouble otrterm = - m_domega_jdT_prtr * (m_Z_pr_tr + 1.0);
    CanteraDouble deltaS_calgmol = c1term + c2term + a3term + a4term + wterm + wrterm + otterm + otrterm;

    return m_units.convertTo(deltaS_calgmol, "J/kmol/K");
}

CanteraDouble PDSS_HKFT::ag(const CanteraDouble temp, const int ifunc) const
{
    static CanteraDouble ag_coeff[3] = { -2.037662, 5.747000E-3, -6.557892E-6};
    if (ifunc == 0) {
        return ag_coeff[0] + ag_coeff[1] * temp + ag_coeff[2] * temp * temp;
    } else if (ifunc == 1) {
        return ag_coeff[1] + ag_coeff[2] * 2.0 * temp;
    }
    if (ifunc != 2) {
        return 0.0;
    }
    return ag_coeff[2] * 2.0;
}

CanteraDouble PDSS_HKFT::bg(const CanteraDouble temp, const int ifunc) const
{
    static CanteraDouble bg_coeff[3] = { 6.107361, -1.074377E-2, 1.268348E-5};
    if (ifunc == 0) {
        return bg_coeff[0] + bg_coeff[1] * temp + bg_coeff[2] * temp * temp;
    }   else if (ifunc == 1) {
        return bg_coeff[1] + bg_coeff[2] * 2.0 * temp;
    }
    if (ifunc != 2) {
        return 0.0;
    }
    return bg_coeff[2] * 2.0;
}

CanteraDouble PDSS_HKFT::f(const CanteraDouble temp, const CanteraDouble pres, const int ifunc) const
{
    static CanteraDouble af_coeff[3] = { 3.666666E1, -0.1504956E-9, 0.5107997E-13};
    CanteraDouble TC = temp - 273.15;
    CanteraDouble presBar = pres / 1.0E5;

    if (TC < 155.0) {
        return 0.0;
    }
    TC = std::min(TC, 355.0);
    if (presBar > 1000.) {
        return 0.0;
    }

    CanteraDouble T1 = (TC-155.0)/300.;
    CanteraDouble p2 = (1000. - presBar) * (1000. - presBar);
    CanteraDouble p3 = (1000. - presBar) * p2;
    CanteraDouble p4 = p2 * p2;
    CanteraDouble fac2 = af_coeff[1] * p3 + af_coeff[2] * p4;
    if (ifunc == 0) {
        return pow(T1,4.8) + af_coeff[0] * pow(T1, 16.0) * fac2;
    } else if (ifunc == 1) {
        return (4.8 * pow(T1,3.8) + 16.0 * af_coeff[0] * pow(T1, 15.0)) / 300. * fac2;
    } else if (ifunc == 2) {
        return (4.8 * 3.8 * pow(T1,2.8) + 16.0 * 15.0 * af_coeff[0] * pow(T1, 14.0)) / (300. * 300.) * fac2;
    } else if (ifunc == 3) {
        CanteraDouble fac1 = pow(T1,4.8) + af_coeff[0] * pow(T1, 16.0);
        fac2 = - (3.0 * af_coeff[1] * p2 + 4.0 * af_coeff[2] * p3)/ 1.0E5;
        return fac1 * fac2;
    } else {
        throw NotImplementedError("PDSS_HKFT::f");
    }
}

CanteraDouble PDSS_HKFT::g(const CanteraDouble temp, const CanteraDouble pres, const int ifunc) const
{
    CanteraDouble afunc = ag(temp, 0);
    CanteraDouble bfunc = bg(temp, 0);
    m_waterSS->setState_TP(temp, pres);
    m_densWaterSS = m_waterSS->density();
    // density in gm cm-3
    CanteraDouble dens = m_densWaterSS * 1.0E-3;
    CanteraDouble gval = afunc * pow((1.0-dens), bfunc);
    if (dens >= 1.0) {
        return 0.0;
    }
    if (ifunc == 0) {
        return gval;
    } else if (ifunc == 1 || ifunc == 2) {
        CanteraDouble afuncdT = ag(temp, 1);
        CanteraDouble bfuncdT = bg(temp, 1);
        CanteraDouble alpha = m_waterSS->thermalExpansionCoeff();

        CanteraDouble fac1 = afuncdT * gval / afunc;
        CanteraDouble fac2 = bfuncdT * gval * log(1.0 - dens);
        CanteraDouble fac3 = gval * alpha * bfunc * dens / (1.0 - dens);

        CanteraDouble dgdt = fac1 + fac2 + fac3;
        if (ifunc == 1) {
            return dgdt;
        }

        CanteraDouble afuncdT2 = ag(temp, 2);
        CanteraDouble bfuncdT2 = bg(temp, 2);
        CanteraDouble dfac1dT = dgdt * afuncdT / afunc + afuncdT2 * gval / afunc
                             -  afuncdT * afuncdT * gval / (afunc * afunc);
        CanteraDouble ddensdT = - alpha * dens;
        CanteraDouble dfac2dT = bfuncdT2 * gval * log(1.0 - dens)
                              + bfuncdT * dgdt * log(1.0 - dens)
                              - bfuncdT * gval /(1.0 - dens) * ddensdT;
        CanteraDouble dalphadT = m_waterSS->dthermalExpansionCoeffdT();
        CanteraDouble dfac3dT = dgdt * alpha * bfunc * dens / (1.0 - dens)
                             + gval * dalphadT * bfunc * dens / (1.0 - dens)
                             + gval * alpha * bfuncdT * dens / (1.0 - dens)
                             + gval * alpha * bfunc * ddensdT / (1.0 - dens)
                             + gval * alpha * bfunc * dens / ((1.0 - dens) * (1.0 - dens)) * ddensdT;

        return dfac1dT + dfac2dT + dfac3dT;
    } else if (ifunc == 3) {
        CanteraDouble beta = m_waterSS->isothermalCompressibility();
        return - bfunc * gval * dens * beta / (1.0 - dens);
    } else {
        throw NotImplementedError("PDSS_HKFT::g");
    }
    return 0.0;
}

CanteraDouble PDSS_HKFT::gstar(const CanteraDouble temp, const CanteraDouble pres, const int ifunc) const
{
    CanteraDouble gval = g(temp, pres, ifunc);
    CanteraDouble fval = f(temp, pres, ifunc);
    CanteraDouble res = gval - fval;
    return res;
}

CanteraDouble PDSS_HKFT::LookupGe(const string& elemName)
{
    size_t iE = m_tp->elementIndex(elemName);
    if (iE == npos) {
        throw CanteraError("PDSS_HKFT::LookupGe", "element " + elemName + " not found");
    }
    CanteraDouble geValue = m_tp->entropyElement298(iE);
    if (geValue == ENTROPY298_UNKNOWN) {
        throw CanteraError("PDSS_HKFT::LookupGe",
                           "element " + elemName + " does not have a supplied entropy298");
    }
    return geValue * -298.15;
}

void PDSS_HKFT::convertDGFormation()
{
    // Ok let's get the element compositions and conversion factors.
    CanteraDouble totalSum = 0.0;
    for (size_t m = 0; m < m_tp->nElements(); m++) {
        CanteraDouble na = m_tp->nAtoms(m_spindex, m);
        if (na > 0.0) {
            totalSum += na * LookupGe(m_tp->elementName(m));
        }
    }
    // Add in the charge
    if (m_charge_j != 0.0) {
        totalSum -= m_charge_j * LookupGe("H");
    }
    CanteraDouble dg = m_units.convertTo(m_deltaG_formation_tr_pr, "J/kmol");
    //! Store the result into an internal variable.
    m_Mu0_tr_pr = dg + totalSum;
}

void PDSS_HKFT::reportParams(size_t& kindex, int& type, CanteraDouble* const c,
                             CanteraDouble& minTemp_, CanteraDouble& maxTemp_,
                             CanteraDouble& refPressure_) const
{
    // Fill in the first part
    PDSS::reportParams(kindex, type, c, minTemp_, maxTemp_,
                       refPressure_);

    c[0] = m_deltaG_formation_tr_pr;
    c[1] = m_deltaH_formation_tr_pr;
    c[2] = m_Mu0_tr_pr;
    c[3] = m_Entrop_tr_pr;
    c[4] = m_a1;
    c[5] = m_a2;
    c[6] = m_a3;
    c[7] = m_a4;
    c[8] = m_c1;
    c[9] = m_c2;
    c[10] = m_omega_pr_tr;
}

}
