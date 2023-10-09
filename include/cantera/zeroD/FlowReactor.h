//! @file FlowReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FLOWREACTOR_H
#define CT_FLOWREACTOR_H

#include "IdealGasReactor.h"

namespace Cantera
{

//! Adiabatic flow in a constant-area duct with homogeneous and heterogeneous reactions
//! @ingroup reactorGroup
class FlowReactor : public IdealGasReactor
{
public:
    FlowReactor() = default;

    string type() const override {
        return "FlowReactor";
    }

    bool isOde() const override {
        return false;
    }

    bool timeIsIndependent() const override {
        return false;
    }

    //! Not implemented; FlowReactor implements getStateDAE() instead.
    void getState(CanteraDouble* y) override {
        throw NotImplementedError("FlowReactor::getState");
    }

    void getStateDae(CanteraDouble* y, CanteraDouble* ydot) override;
    void initialize(CanteraDouble t0=0.0) override;
    void syncState() override;
    void updateState(CanteraDouble* y) override;

    //! Not implemented; FlowReactor implements evalDae() instead.
    void eval(CanteraDouble t, CanteraDouble* LHS, CanteraDouble* RHS) override {
        throw NotImplementedError("FlowReactor::eval");
    }

    void evalDae(CanteraDouble t, CanteraDouble* y, CanteraDouble* ydot, CanteraDouble* residual) override;

    void getConstraints(CanteraDouble* constraints) override;

    //! Set the mass flow rate through the reactor [kg/s]
    void setMassFlowRate(CanteraDouble mdot);

    //! The current gas speed in the reactor [m/s]
    CanteraDouble speed() const {
        return m_u;
    }

    //! The cross-sectional area of the reactor [m^2]
    CanteraDouble area() const {
        return m_area;
    }

    //! @deprecated To be removed after %Cantera 3.0. Access distance through the
    //!     ReactorNet object
    CanteraDouble distance() const;

    //! Sets the area of the reactor [m^2]
    void setArea(CanteraDouble area);

    //! The ratio of the reactor's surface area to volume ratio [m^-1]
    //! @note If the surface area to volume ratio is unspecified by the user,
    //!       this will be calculated assuming the reactor is a cylinder.
    CanteraDouble surfaceAreaToVolumeRatio() const;

    //! Set the reactor's surface area to volume ratio [m^-1]
    void setSurfaceAreaToVolumeRatio(CanteraDouble sa_to_vol) {
        m_sa_to_vol = sa_to_vol;
    }

    //! Get the steady state tolerances used to determine the initial state for
    //! surface coverages
    CanteraDouble inletSurfaceAtol() const {
        return m_ss_atol;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setInletSurfaceAtol(CanteraDouble atol) {
        m_ss_atol = atol;
    }

    //! Get the steady state tolerances used to determine the initial state for
    //! surface coverages
    CanteraDouble inletSurfaceRtol() const {
        return m_ss_rtol;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setInletSurfaceRtol(CanteraDouble rtol) {
        m_ss_rtol = rtol;
    }

    //! Get the steady state tolerances used to determine the initial state for
    //! surface coverages
    CanteraDouble inletSurfaceMaxSteps() const {
        return m_max_ss_steps;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setInletSurfaceMaxSteps(int max_steps) {
        m_max_ss_steps = max_steps;
    }

    //! Get the steady state tolerances used to determine the initial state for
    //! surface coverages
    CanteraDouble inletSurfaceMaxErrorFailures() const {
        return m_max_ss_error_fails;
    }

    //! Set the steady state tolerances used to determine the initial state for
    //! surface coverages
    void setInletSurfaceMaxErrorFailures(int max_fails) {
        m_max_ss_error_fails = max_fails;
    }

    //! Return the index in the solution vector for this reactor of the component named
    //! *nm*. Possible values for *nm* are "density", "speed", "pressure",
    //! "temperature", the name of a homogeneous phase species, or the name of a surface
    //! species.
    size_t componentIndex(const string& nm) const override;

    string componentName(size_t k) override;

    void updateSurfaceState(CanteraDouble* y) override;

protected:
    //! Density [kg/m^3]. First component of the state vector.
    CanteraDouble m_rho = NAN;
    //! Axial velocity [m/s]. Second component of the state vector.
    CanteraDouble m_u = -1.0;
    //! Pressure [Pa]. Third component of the state vector.
    CanteraDouble m_P = NAN;
    //! Temperature [K]. Fourth component of the state vector.
    CanteraDouble m_T = NAN;
    //! offset to the species equations
    const size_t m_offset_Y = 4;
    //! reactor area [m^2]
    CanteraDouble m_area = 1.0;
    //! reactor surface area to volume ratio [m^-1]
    CanteraDouble m_sa_to_vol = -1.0;
    //! temporary storage for surface species production rates
    vector<CanteraDouble> m_sdot_temp;
    //! temporary storage for species partial molar enthalpies
    vector<CanteraDouble> m_hk;
    //! steady-state relative tolerance, used to determine initial surface coverages
    CanteraDouble m_ss_rtol = 1e-7;
    //! steady-state absolute tolerance, used to determine initial surface coverages
    CanteraDouble m_ss_atol = 1e-14;
    //! maximum number of steady-state coverage integrator-steps
    int m_max_ss_steps = 20000;
    //! maximum number of steady-state integrator error test failures
    int m_max_ss_error_fails = 10;
};
}

#endif
