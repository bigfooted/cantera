/**
 *  @file SpeciesThermoInterpType.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/thermo/PDSS.h"

namespace Cantera
{

SpeciesThermoInterpType::SpeciesThermoInterpType(CanteraDouble tlow,
                                                 CanteraDouble thigh,
                                                 CanteraDouble pref) :
    m_lowT(tlow),
    m_highT(thigh),
    m_Pref(pref)
{
}

void SpeciesThermoInterpType::updateProperties(const CanteraDouble* tempPoly,
        CanteraDouble* cp_R, CanteraDouble* h_RT, CanteraDouble* s_R) const
{
    CanteraDouble T = tempPoly[0];
    updatePropertiesTemp(T, cp_R, h_RT, s_R);
}

void SpeciesThermoInterpType::updatePropertiesTemp(const CanteraDouble temp,
        CanteraDouble* cp_R, CanteraDouble* h_RT, CanteraDouble* s_R) const
{
    throw NotImplementedError("SpeciesThermoInterpType::updatePropertiesTemp");
}

size_t SpeciesThermoInterpType::nCoeffs() const
{
    throw NotImplementedError("SpeciesThermoInterpType::nCoeffs");
}

void SpeciesThermoInterpType::reportParameters(size_t& index, int& type,
        CanteraDouble& minTemp, CanteraDouble& maxTemp, CanteraDouble& refPressure,
        CanteraDouble* const coeffs) const
{
    throw NotImplementedError("SpeciesThermoInterpType::reportParameters");
}

AnyMap SpeciesThermoInterpType::parameters(bool withInput) const
{
    AnyMap out;
    getParameters(out);
    if (withInput) {
        out.update(m_input);
    }
    return out;
}

void SpeciesThermoInterpType::getParameters(AnyMap& thermo) const
{
    if (m_Pref != OneAtm && reportType() != 0) {
        thermo["reference-pressure"].setQuantity(m_Pref, "Pa");
    }
}

CanteraDouble SpeciesThermoInterpType::reportHf298(CanteraDouble* const h298) const
{
    throw NotImplementedError("SpeciesThermoInterpType::reportHf298");
}

void SpeciesThermoInterpType::modifyOneHf298(const size_t k, const CanteraDouble Hf298New)
{
    throw NotImplementedError("SpeciesThermoInterpType::modifyOneHf298");
}

const AnyMap& SpeciesThermoInterpType::input() const
{
    return m_input;
}

AnyMap& SpeciesThermoInterpType::input()
{
    return m_input;
}

}
