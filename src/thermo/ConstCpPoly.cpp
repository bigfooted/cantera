/**
 *  @file ConstCpPoly.cpp
 * Declarations for the @link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType @endlink object that
 * employs a constant heat capacity assumption (see @ref spthermo and
 * @link Cantera::ConstCpPoly ConstCpPoly @endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/ConstCpPoly.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

ConstCpPoly::ConstCpPoly()
    : SpeciesThermoInterpType(0.0, std::numeric_limits<CanteraDouble>::infinity(), 0.0)
{
}

ConstCpPoly::ConstCpPoly(CanteraDouble tlow, CanteraDouble thigh, CanteraDouble pref,
                         const CanteraDouble* coeffs) :
    SpeciesThermoInterpType(tlow, thigh, pref)
{
    setParameters(coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
}

void ConstCpPoly::setParameters(CanteraDouble t0, CanteraDouble h0, CanteraDouble s0, CanteraDouble cp0)
{
    m_t0 = t0;
    m_logt0 = log(m_t0);
    m_cp0_R = cp0 / GasConstant;
    m_h0_R = h0 / GasConstant;
    m_s0_R = s0 / GasConstant;
}

void ConstCpPoly::updateProperties(const CanteraDouble* tt,
                                   CanteraDouble* cp_R,
                                   CanteraDouble* h_RT,
                                   CanteraDouble* s_R) const
{
    CanteraDouble t = *tt;
    CanteraDouble logt = log(t);
    CanteraDouble rt = 1.0/t;
    *cp_R = m_cp0_R;
    *h_RT = rt*(m_h0_R + (t - m_t0) * m_cp0_R);
    *s_R = m_s0_R + m_cp0_R * (logt - m_logt0);
}

void ConstCpPoly::updatePropertiesTemp(const CanteraDouble temp,
                                       CanteraDouble* cp_R,
                                       CanteraDouble* h_RT,
                                       CanteraDouble* s_R) const
{
    CanteraDouble logt = log(temp);
    CanteraDouble rt = 1.0/temp;
    *cp_R = m_cp0_R;
    *h_RT = rt*(m_h0_R + (temp - m_t0) * m_cp0_R);
    *s_R = m_s0_R + m_cp0_R * (logt - m_logt0);
}

void ConstCpPoly::reportParameters(size_t& n, int& type, CanteraDouble& tlow, CanteraDouble& thigh,
                                   CanteraDouble& pref, CanteraDouble* const coeffs) const
{
    n = 0;
    type = CONSTANT_CP;
    tlow = m_lowT;
    thigh = m_highT;
    pref = m_Pref;
    coeffs[0] = m_t0;
    coeffs[1] = m_h0_R * GasConstant;
    coeffs[2] = m_s0_R * GasConstant;
    coeffs[3] = m_cp0_R * GasConstant;
}

void ConstCpPoly::getParameters(AnyMap& thermo) const
{
    thermo["model"] = "constant-cp";
    SpeciesThermoInterpType::getParameters(thermo);
    thermo["T0"].setQuantity(m_t0, "K");
    thermo["h0"].setQuantity(m_h0_R * GasConstant, "J/kmol");
    thermo["s0"].setQuantity(m_s0_R * GasConstant, "J/kmol/K");
    thermo["cp0"].setQuantity(m_cp0_R * GasConstant, "J/kmol/K");
}

CanteraDouble ConstCpPoly::reportHf298(CanteraDouble* const h298) const
{
    CanteraDouble temp = 298.15;
    CanteraDouble h = GasConstant * (m_h0_R + (temp - m_t0) * m_cp0_R);
    if (h298) {
        *h298 = h;
    }
    return h;
}

void ConstCpPoly::modifyOneHf298(const size_t k, const CanteraDouble Hf298New)
{
    CanteraDouble hnow = reportHf298();
    CanteraDouble delH = Hf298New - hnow;
    m_h0_R += delH / GasConstant;
}

void ConstCpPoly::resetHf298() {
    m_h0_R = m_h0_R_orig;
}

}
