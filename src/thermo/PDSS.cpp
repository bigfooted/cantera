/**
 * @file PDSS.cpp
 * Implementation of a pressure dependent standard state
 * virtual function
 * (see class @link Cantera::PDSS PDSS@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/global.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/VPStandardStateTP.h"

namespace Cantera
{

CanteraDouble PDSS::enthalpy_mole() const
{
    throw NotImplementedError("PDSS::enthalpy_mole");
}

CanteraDouble PDSS::enthalpy_RT() const
{
    throw NotImplementedError("PDSS::enthalpy_RT");
}

CanteraDouble PDSS::intEnergy_mole() const
{
    throw NotImplementedError("PDSS::intEnergy_mole");
}

CanteraDouble PDSS::entropy_mole() const
{
    throw NotImplementedError("PDSS::entropy_mole");
}

CanteraDouble PDSS::entropy_R() const
{
    throw NotImplementedError("PDSS::entropy_R");
}

CanteraDouble PDSS::gibbs_mole() const
{
    throw NotImplementedError("PDSS::gibbs_mole");
}

CanteraDouble PDSS::gibbs_RT() const
{
    throw NotImplementedError("PDSS::gibbs_RT");
}

CanteraDouble PDSS::cp_mole() const
{
    throw NotImplementedError("PDSS::cp_mole");
}

CanteraDouble PDSS::cp_R() const
{
    throw NotImplementedError("PDSS::cp_R");
}

CanteraDouble PDSS::molarVolume() const
{
    throw NotImplementedError("PDSS::molarVolume");
}

CanteraDouble PDSS::density() const
{
    throw NotImplementedError("PDSS::density");
}

CanteraDouble PDSS::cv_mole() const
{
    throw NotImplementedError("PDSS::cv_mole");
}

CanteraDouble PDSS::gibbs_RT_ref() const
{
    throw NotImplementedError("PDSS::gibbs_RT_ref");
}

CanteraDouble PDSS::enthalpy_RT_ref() const
{
    throw NotImplementedError("PDSS::enthalpy_RT_ref");
}

CanteraDouble PDSS::entropy_R_ref() const
{
    throw NotImplementedError("PDSS::entropy_RT_ref");
}

CanteraDouble PDSS::cp_R_ref() const
{
    throw NotImplementedError("PDSS::entropy_RT_ref");
}

CanteraDouble PDSS::molarVolume_ref() const
{
    throw NotImplementedError("PDSS::molarVolume_ref");
}

CanteraDouble PDSS::pressure() const
{
    return m_pres;
}

CanteraDouble PDSS::thermalExpansionCoeff() const
{
    throw NotImplementedError("PDSS::thermalExpansionCoeff");
}

CanteraDouble PDSS::critTemperature() const
{
    throw NotImplementedError("PDSS::critTemperature");
}

CanteraDouble PDSS::critPressure() const
{
    throw NotImplementedError("PDSS::critPressure");
}

CanteraDouble PDSS::critDensity() const
{
    throw NotImplementedError("PDSS::critDensity");
}

void PDSS::setPressure(CanteraDouble pres)
{
    m_pres = pres;
}

CanteraDouble PDSS::temperature() const
{
    return m_temp;
}

void PDSS::setTemperature(CanteraDouble temp)
{
    m_temp = temp;
}

CanteraDouble PDSS::molecularWeight() const
{
    return m_mw;
}
void PDSS::setMolecularWeight(CanteraDouble mw)
{
    m_mw = mw;
}

void PDSS::setState_TP(CanteraDouble temp, CanteraDouble pres)
{
    throw NotImplementedError("PDSS::setState_TP");
}

CanteraDouble PDSS::satPressure(CanteraDouble t)
{
    throw NotImplementedError("PDSS::satPressure");
}

// PDSS_Molar methods

CanteraDouble PDSS_Molar::enthalpy_RT() const
{
    return enthalpy_mole() / (GasConstant * temperature());
}

CanteraDouble PDSS_Molar::entropy_R() const
{
    return entropy_mole() / GasConstant;
}

CanteraDouble PDSS_Molar::gibbs_RT() const
{
    return gibbs_mole() / (GasConstant * temperature());
}

CanteraDouble PDSS_Molar::cp_R() const
{
    return cp_mole() / GasConstant;
}

// PDSS_Nondimensional methods

PDSS_Nondimensional::PDSS_Nondimensional()
    : m_h0_RT(0.0)
    , m_cp0_R(0.0)
    , m_s0_R(0.0)
    , m_g0_RT(0.0)
    , m_V0(0.0)
    , m_hss_RT(0.0)
    , m_cpss_R(0.0)
    , m_sss_R(0.0)
    , m_gss_RT(0.0)
    , m_Vss(0.0)
{
}

CanteraDouble PDSS_Nondimensional::enthalpy_mole() const
{
    return enthalpy_RT() * GasConstant * temperature();
}

CanteraDouble PDSS_Nondimensional::entropy_mole() const
{
    return entropy_R() * GasConstant;
}

CanteraDouble PDSS_Nondimensional::gibbs_mole() const
{
    return gibbs_RT() * GasConstant * temperature();
}

CanteraDouble PDSS_Nondimensional::cp_mole() const
{
    return cp_R() * GasConstant;
}

CanteraDouble PDSS_Nondimensional::gibbs_RT_ref() const
{
    return m_g0_RT;
}

CanteraDouble PDSS_Nondimensional::enthalpy_RT_ref() const
{
    return m_h0_RT;
}

CanteraDouble PDSS_Nondimensional::entropy_R_ref() const
{
    return m_s0_R;
}

CanteraDouble PDSS_Nondimensional::cp_R_ref() const
{
    return m_cp0_R;
}

CanteraDouble PDSS_Nondimensional::molarVolume_ref() const
{
    return m_V0;
}

CanteraDouble PDSS_Nondimensional::enthalpy_RT() const
{
    return m_hss_RT;
}

CanteraDouble PDSS_Nondimensional::entropy_R() const
{
    return m_sss_R;
}

CanteraDouble PDSS_Nondimensional::gibbs_RT() const
{
    return m_gss_RT;
}

CanteraDouble PDSS_Nondimensional::cp_R() const
{
    return m_cpss_R;
}

CanteraDouble PDSS_Nondimensional::molarVolume() const
{
    return m_Vss;
}

CanteraDouble PDSS_Nondimensional::density() const
{
    return m_mw / m_Vss;
}

}
