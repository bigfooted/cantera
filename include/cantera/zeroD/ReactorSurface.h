//! @file ReactorSurface.h Header file for class ReactorSurface

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTOR_SURFACE_H
#define CT_REACTOR_SURFACE_H

#include "cantera/zeroD/ReactorBase.h"

namespace Cantera
{

class Kinetics;
class SurfPhase;

//! A surface where reactions can occur that is in contact with the bulk fluid of a
//! Reactor.
//! @ingroup wallGroup
class ReactorSurface
{
public:
    ReactorSurface() = default;
    virtual ~ReactorSurface() = default;
    ReactorSurface(const ReactorSurface&) = delete;
    ReactorSurface& operator=(const ReactorSurface&) = delete;

    //! Returns the surface area [m^2]
    CanteraDouble area() const;

    //! Set the surface area [m^2]
    void setArea(CanteraDouble a);

    //! Accessor for the SurfPhase object
    SurfPhase* thermo() {
        return m_thermo;
    }

    //! Accessor for the InterfaceKinetics object
    Kinetics* kinetics() {
        return m_kinetics;
    }

    //! Set the InterfaceKinetics object for this surface
    void setKinetics(Kinetics* kin);

    //! Set the reactor that this Surface interacts with
    void setReactor(ReactorBase* reactor);

    //! Number of sensitivity parameters associated with reactions on this
    //! surface
    size_t nSensParams() const {
        return m_params.size();
    }

    //! Set the surface coverages. Array `cov` has length equal to the number of
    //! surface species.
    void setCoverages(const CanteraDouble* cov);

    //! Set the surface coverages by name
    void setCoverages(const Composition& cov);

    //! Set the surface coverages by name
    void setCoverages(const string& cov);

    //! Get the surface coverages. Array `cov` should have length equal to the
    //! number of surface species.
    void getCoverages(CanteraDouble* cov) const;

    //! Set the coverages and temperature in the surface phase object to the
    //! values for this surface. The temperature is set to match the bulk phase
    //! of the attached Reactor.
    void syncState();

    //! Enable calculation of sensitivities with respect to the rate constant
    //! for reaction `i`.
    void addSensitivityReaction(size_t i);

    //! Set reaction rate multipliers. `params` is the global vector of
    //! sensitivity parameters. This function is called within
    //! ReactorNet::eval() before the reaction rates are evaluated.
    void setSensitivityParameters(const CanteraDouble* params);

    //! Set reaction rate multipliers back to their initial values. This
    //! function is called within ReactorNet::eval() after all rates have been
    //! evaluated.
    void resetSensitivityParameters();

protected:
    CanteraDouble m_area = 1.0;

    SurfPhase* m_thermo = nullptr;
    Kinetics* m_kinetics = nullptr;
    ReactorBase* m_reactor = nullptr;
    vector<CanteraDouble> m_cov;
    vector<SensitivityParameter> m_params;
};

}

#endif
