/**
 * @file ctmultiphase.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CTC_MULTIPHASE_H
#define CTC_MULTIPHASE_H

#include "clib_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

    CANTERA_CAPI int mix_new();
    CANTERA_CAPI int mix_del(int i);
    CANTERA_CAPI int ct_clearMix();
    CANTERA_CAPI int mix_addPhase(int i, int j, CanteraDouble moles);
    CANTERA_CAPI int mix_init(int i);
    CANTERA_CAPI int mix_updatePhases(int i);
    CANTERA_CAPI size_t mix_nElements(int i);
    CANTERA_CAPI size_t mix_elementIndex(int i, const char* name);
    CANTERA_CAPI size_t mix_speciesIndex(int i, int k, int p);
    CANTERA_CAPI size_t mix_nSpecies(int i);
    CANTERA_CAPI int mix_setTemperature(int i, CanteraDouble t);
    CANTERA_CAPI CanteraDouble mix_temperature(int i);
    CANTERA_CAPI CanteraDouble mix_minTemp(int i);
    CANTERA_CAPI CanteraDouble mix_maxTemp(int i);
    CANTERA_CAPI CanteraDouble mix_charge(int i);
    CANTERA_CAPI CanteraDouble mix_phaseCharge(int i, int p);
    CANTERA_CAPI int mix_setPressure(int i, CanteraDouble p);
    CANTERA_CAPI CanteraDouble mix_pressure(int i);
    CANTERA_CAPI CanteraDouble mix_nAtoms(int i, int k, int m);
    CANTERA_CAPI size_t mix_nPhases(int i);
    CANTERA_CAPI CanteraDouble mix_phaseMoles(int i, int n);
    CANTERA_CAPI int mix_setPhaseMoles(int i, int n, CanteraDouble v);
    CANTERA_CAPI int mix_setMoles(int i, size_t nlen, const CanteraDouble* n);
    CANTERA_CAPI int mix_setMolesByName(int i, const char* n);
    CANTERA_CAPI CanteraDouble mix_speciesMoles(int i, int k);
    CANTERA_CAPI CanteraDouble mix_elementMoles(int i, int m);
    CANTERA_CAPI CanteraDouble mix_equilibrate(int i, const char* XY, CanteraDouble err,
                                        int maxsteps, int maxiter, int loglevel);
    CANTERA_CAPI int mix_getChemPotentials(int i, size_t lenmu, CanteraDouble* mu);
    CANTERA_CAPI CanteraDouble mix_enthalpy(int i);
    CANTERA_CAPI CanteraDouble mix_entropy(int i);
    CANTERA_CAPI CanteraDouble mix_gibbs(int i);
    CANTERA_CAPI CanteraDouble mix_cp(int i);
    CANTERA_CAPI CanteraDouble mix_volume(int i);

    CANTERA_CAPI size_t mix_speciesPhaseIndex(int i, int k);
    CANTERA_CAPI CanteraDouble mix_moleFraction(int i, int k);

#ifdef __cplusplus
}
#endif

#endif
