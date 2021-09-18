//! @file ReactionData.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

void ArrheniusData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    update(bulk.temperature());
}

void BlowersMaselData::update(double T)
{
    temperature = T;
    logT = std::log(T);
    recipT = 1./T;
}

void BlowersMaselData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    update(bulk.temperature());
    bulk.getPartialMolarEnthalpies(m_grt.data());
    kin.getReactionDelta(m_grt.data(), dH.data());
}

void FalloffData::update(double T)
{
    temperature = T;
    logT = std::log(T);
    recipT = 1./T;
}

void FalloffData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    update(bulk.temperature());
}

void PlogData::update(double T)
{
    throw CanteraError("PlogData::update",
        "Missing state information: reaction type requires pressure.");
}

void PlogData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    update(bulk.temperature(), bulk.pressure());
}

void ChebyshevData::update(double T)
{
    throw CanteraError("ChebyshevData::update",
        "Missing state information: reaction type requires pressure.");
}

void ChebyshevData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    update(bulk.temperature(), bulk.pressure());
}

void CustomFunc1Data::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    temperature = bulk.temperature();
}

}
