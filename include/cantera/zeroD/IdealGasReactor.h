//! @file IdealGasReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASREACTOR_H
#define CT_IDEALGASREACTOR_H

#include "Reactor.h"

namespace Cantera
{

/**
 * Class IdealGasReactor is a class for stirred reactors that is specifically
 * optimized for ideal gases. In this formulation, temperature replaces the
 * total internal energy as a state variable.
 * @ingroup reactorGroup
 */
class IdealGasReactor : public Reactor
{
public:
    IdealGasReactor() {}

    string type() const override {
        return "IdealGasReactor";
    }

    void setThermoMgr(ThermoPhase& thermo) override;

    void getState(CanteraDouble* y) override;

    void initialize(CanteraDouble t0=0.0) override;

    void eval(CanteraDouble t, CanteraDouble* LHS, CanteraDouble* RHS) override;

    void updateState(CanteraDouble* y) override;

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass",
    //! "volume", "temperature", the name of a homogeneous phase species, or the
    //! name of a surface species.
    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;

protected:
    vector<CanteraDouble> m_uk; //!< Species molar internal energies
};

}

#endif
