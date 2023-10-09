/**
 * @file ctreactor.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CTC_REACTOR_H
#define CTC_REACTOR_H

#include "clib_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

    CANTERA_CAPI int reactor_new(const char* type);
    CANTERA_CAPI int reactor_del(int i);
    CANTERA_CAPI int reactor_setInitialVolume(int i, CanteraDouble v);
    CANTERA_CAPI int reactor_setChemistry(int i, int cflag);
    CANTERA_CAPI int reactor_setEnergy(int i, int eflag);
    CANTERA_CAPI int reactor_setThermoMgr(int i, int n);
    CANTERA_CAPI int reactor_setKineticsMgr(int i, int n);
    CANTERA_CAPI int reactor_insert(int i, int n);
    CANTERA_CAPI CanteraDouble reactor_mass(int i);
    CANTERA_CAPI CanteraDouble reactor_volume(int i);
    CANTERA_CAPI CanteraDouble reactor_density(int i);
    CANTERA_CAPI CanteraDouble reactor_temperature(int i);
    CANTERA_CAPI CanteraDouble reactor_enthalpy_mass(int i);
    CANTERA_CAPI CanteraDouble reactor_intEnergy_mass(int i);
    CANTERA_CAPI CanteraDouble reactor_pressure(int i);
    CANTERA_CAPI CanteraDouble reactor_massFraction(int i, int k);
    CANTERA_CAPI size_t reactor_nSensParams(int i);
    CANTERA_CAPI int reactor_addSensitivityReaction(int i, int rxn);
    CANTERA_CAPI int flowReactor_setMassFlowRate(int i, CanteraDouble mdot);

    CANTERA_CAPI int reactornet_new();
    CANTERA_CAPI int reactornet_del(int i);
    CANTERA_CAPI int reactornet_setInitialTime(int i, CanteraDouble t);
    CANTERA_CAPI int reactornet_setMaxTimeStep(int i, CanteraDouble maxstep);
    CANTERA_CAPI int reactornet_setTolerances(int i, CanteraDouble rtol, CanteraDouble atol);
    CANTERA_CAPI int reactornet_setSensitivityTolerances(int i, CanteraDouble rtol, CanteraDouble atol);
    CANTERA_CAPI int reactornet_addreactor(int i, int n);
    CANTERA_CAPI int reactornet_advance(int i, CanteraDouble t);
    CANTERA_CAPI CanteraDouble reactornet_step(int i);
    CANTERA_CAPI CanteraDouble reactornet_time(int i);
    CANTERA_CAPI CanteraDouble reactornet_rtol(int i);
    CANTERA_CAPI CanteraDouble reactornet_atol(int i);
    CANTERA_CAPI CanteraDouble reactornet_sensitivity(int i, const char* v, int p, int r);

    CANTERA_CAPI int flowdev_new(const char* type);
    CANTERA_CAPI int flowdev_del(int i);
    CANTERA_CAPI int flowdev_install(int i, int n, int m);
    //! @deprecated To be removed after %Cantera 3.0; replaced by flowdev_setPrimary
    CANTERA_CAPI int flowdev_setMaster(int i, int n);
    CANTERA_CAPI int flowdev_setPrimary(int i, int n);
    CANTERA_CAPI CanteraDouble flowdev_massFlowRate(int i);
    CANTERA_CAPI int flowdev_setMassFlowCoeff(int i, CanteraDouble v);
    CANTERA_CAPI int flowdev_setValveCoeff(int i, CanteraDouble v);
    CANTERA_CAPI int flowdev_setPressureCoeff(int i, CanteraDouble v);
    CANTERA_CAPI int flowdev_setPressureFunction(int i, int n);
    CANTERA_CAPI int flowdev_setTimeFunction(int i, int n);

    CANTERA_CAPI int wall_new(const char* type);
    CANTERA_CAPI int wall_del(int i);
    CANTERA_CAPI int wall_install(int i, int n, int m);
    //! @deprecated Only used by traditional MATLAB toolbox
    CANTERA_CAPI CanteraDouble wall_vdot(int i, CanteraDouble t);
    CANTERA_CAPI CanteraDouble wall_expansionRate(int i);
    //! @deprecated Only used by traditional MATLAB toolbox
    CANTERA_CAPI CanteraDouble wall_Q(int i, CanteraDouble t);
    CANTERA_CAPI CanteraDouble wall_heatRate(int i);
    CANTERA_CAPI CanteraDouble wall_area(int i);
    CANTERA_CAPI int wall_setArea(int i, CanteraDouble v);
    CANTERA_CAPI int wall_setThermalResistance(int i, CanteraDouble rth);
    CANTERA_CAPI int wall_setHeatTransferCoeff(int i, CanteraDouble u);
    CANTERA_CAPI int wall_setHeatFlux(int i, int n);
    CANTERA_CAPI int wall_setExpansionRateCoeff(int i, CanteraDouble k);
    CANTERA_CAPI int wall_setVelocity(int i, int n);
    CANTERA_CAPI int wall_setEmissivity(int i, CanteraDouble epsilon);
    CANTERA_CAPI int wall_ready(int i);

    CANTERA_CAPI int reactorsurface_new(int type);
    CANTERA_CAPI int reactorsurface_del(int i);
    CANTERA_CAPI int reactorsurface_install(int i, int n);
    CANTERA_CAPI int reactorsurface_setkinetics(int i, int n);
    CANTERA_CAPI CanteraDouble reactorsurface_area(int i);
    CANTERA_CAPI int reactorsurface_setArea(int i, CanteraDouble v);
    CANTERA_CAPI int reactorsurface_addSensitivityReaction(int i, int rxn);

    CANTERA_CAPI int ct_clearReactors();

#ifdef __cplusplus
}
#endif

#endif
