/**
 * @file ct.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CTC_CT_H
#define CTC_CT_H

#include "clib_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

    CANTERA_CAPI int ct_appdelete();

    CANTERA_CAPI int soln_newSolution(const char* infile,
                                      const char* name,
                                      const char* transport);
    CANTERA_CAPI int soln_newInterface(const char* infile,
                                       const char* name,
                                       int na,
                                       const int* adjacent);
    CANTERA_CAPI int soln_del(int n); //!< note that linked objects are deleted as well
    CANTERA_CAPI int soln_name(int n, int buflen, char* buf);
    CANTERA_CAPI int soln_thermo(int n);
    CANTERA_CAPI int soln_kinetics(int n);
    CANTERA_CAPI int soln_transport(int n);
    //! note that soln_setTransportModel deletes the previous transport model
    CANTERA_CAPI int soln_setTransportModel(int n, const char* model);
    CANTERA_CAPI size_t soln_nAdjacent(int n);
    CANTERA_CAPI int soln_adjacent(int n, int a);

    CANTERA_CAPI int thermo_newFromFile(const char* filename, const char* phasename);
    CANTERA_CAPI int thermo_del(int n);
    CANTERA_CAPI size_t thermo_nElements(int n);
    CANTERA_CAPI size_t thermo_nSpecies(int n);
    CANTERA_CAPI CanteraDouble thermo_temperature(int n);
    CANTERA_CAPI int thermo_setTemperature(int n, CanteraDouble t);
    CANTERA_CAPI CanteraDouble thermo_density(int n);
    CANTERA_CAPI int thermo_setDensity(int n, CanteraDouble rho);
    CANTERA_CAPI CanteraDouble thermo_molarDensity(int n);
    CANTERA_CAPI int thermo_setMolarDensity(int n, CanteraDouble ndens);
    CANTERA_CAPI CanteraDouble thermo_meanMolecularWeight(int n);
    CANTERA_CAPI CanteraDouble thermo_moleFraction(int n, size_t k);
    CANTERA_CAPI CanteraDouble thermo_massFraction(int n, size_t k);
    CANTERA_CAPI int thermo_getMoleFractions(int n, size_t lenx, CanteraDouble* x);
    CANTERA_CAPI int thermo_getMassFractions(int n, size_t leny, CanteraDouble* y);
    CANTERA_CAPI int thermo_setMoleFractions(int n, size_t lenx,
                                            CanteraDouble* x, int norm);
    CANTERA_CAPI int thermo_setMassFractions(int n, size_t leny,
                                            CanteraDouble* y, int norm);
    CANTERA_CAPI int thermo_setMoleFractionsByName(int n, const char* x);
    CANTERA_CAPI int thermo_setMassFractionsByName(int n, const char* y);
    CANTERA_CAPI int thermo_getAtomicWeights(int n, size_t lenm, CanteraDouble* atw);
    CANTERA_CAPI int thermo_getMolecularWeights(int n, size_t lenm, CanteraDouble* mw);
    CANTERA_CAPI int thermo_getCharges(int n, size_t lenm, CanteraDouble* sc);
    CANTERA_CAPI int thermo_getElementName(int n, size_t k, size_t lennm, char* nm);
    CANTERA_CAPI int thermo_getSpeciesName(int n, size_t m, size_t lennm, char* nm);
    CANTERA_CAPI int thermo_getName(int n, size_t lennm, char* nm);
    CANTERA_CAPI int thermo_setName(int n, const char* nm);
    CANTERA_CAPI size_t thermo_elementIndex(int n, const char* nm);
    CANTERA_CAPI size_t thermo_speciesIndex(int n, const char* nm);
    CANTERA_CAPI int thermo_report(int nth,
                                  int ibuf, char* buf, int show_thermo);
    CANTERA_CAPI int thermo_print(int nth, int show_thermo, CanteraDouble threshold);
    CANTERA_CAPI CanteraDouble thermo_nAtoms(int n, size_t k, size_t m);
    CANTERA_CAPI int thermo_addElement(int n, const char* name, CanteraDouble weight);
    CANTERA_CAPI int thermo_getEosType(int n, size_t leneos, char* eos);
    CANTERA_CAPI CanteraDouble thermo_refPressure(int n);
    CANTERA_CAPI CanteraDouble thermo_minTemp(int n, int k);
    CANTERA_CAPI CanteraDouble thermo_maxTemp(int n, int k);
    CANTERA_CAPI CanteraDouble thermo_enthalpy_mole(int n);
    CANTERA_CAPI CanteraDouble thermo_intEnergy_mole(int n);
    CANTERA_CAPI CanteraDouble thermo_entropy_mole(int n);
    CANTERA_CAPI CanteraDouble thermo_gibbs_mole(int n);
    CANTERA_CAPI CanteraDouble thermo_cp_mole(int n);
    CANTERA_CAPI CanteraDouble thermo_cv_mole(int n);
    CANTERA_CAPI CanteraDouble thermo_pressure(int n);
    CANTERA_CAPI int thermo_setPressure(int n, CanteraDouble p);
    CANTERA_CAPI CanteraDouble thermo_enthalpy_mass(int n);
    CANTERA_CAPI CanteraDouble thermo_intEnergy_mass(int n);
    CANTERA_CAPI CanteraDouble thermo_entropy_mass(int n);
    CANTERA_CAPI CanteraDouble thermo_gibbs_mass(int n);
    CANTERA_CAPI CanteraDouble thermo_cp_mass(int n);
    CANTERA_CAPI CanteraDouble thermo_cv_mass(int n);
    CANTERA_CAPI CanteraDouble thermo_electricPotential(int n);
    CANTERA_CAPI CanteraDouble thermo_thermalExpansionCoeff(int n);
    CANTERA_CAPI CanteraDouble thermo_isothermalCompressibility(int n);
    CANTERA_CAPI int thermo_chemPotentials(int n, size_t lenm, CanteraDouble* murt);
    CANTERA_CAPI int thermo_getEnthalpies_RT(int n, size_t lenm, CanteraDouble* h_rt);
    CANTERA_CAPI int thermo_getEntropies_R(int n, size_t lenm, CanteraDouble* s_r);
    CANTERA_CAPI int thermo_getCp_R(int n, size_t lenm, CanteraDouble* cp_r);
    CANTERA_CAPI int thermo_setElectricPotential(int n, CanteraDouble v);
    CANTERA_CAPI int thermo_set_TP(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_set_TD(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_set_RP(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_set_DP(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_set_HP(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_set_UV(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_set_SV(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_set_SP(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_set_ST(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_set_TV(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_set_PV(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_set_UP(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_set_VH(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_set_TH(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_set_SH(int n, CanteraDouble* vals);
    CANTERA_CAPI int thermo_equilibrate(int n, const char* XY, int solver,
                                        CanteraDouble rtol, int maxsteps, int maxiter,
                                        int loglevel);

    CANTERA_CAPI CanteraDouble thermo_critTemperature(int n);
    CANTERA_CAPI CanteraDouble thermo_critPressure(int n);
    CANTERA_CAPI CanteraDouble thermo_critDensity(int n);
    CANTERA_CAPI CanteraDouble thermo_vaporFraction(int n);
    CANTERA_CAPI CanteraDouble thermo_satTemperature(int n, CanteraDouble p);
    CANTERA_CAPI CanteraDouble thermo_satPressure(int n, CanteraDouble t);
    CANTERA_CAPI int thermo_setState_Psat(int n, CanteraDouble p, CanteraDouble x);
    CANTERA_CAPI int thermo_setState_Tsat(int n, CanteraDouble t, CanteraDouble x);

    //! @since Starting in %Cantera 3.0, the "phasename" argument should be blank
    CANTERA_CAPI int kin_newFromFile(const char* filename, const char* phasename,
                                     int reactingPhase, int neighbor1, int neighbor2,
                                     int neighbor3, int neighbor4);
    CANTERA_CAPI int kin_del(int n);
    CANTERA_CAPI size_t kin_nSpecies(int n);
    CANTERA_CAPI size_t kin_nReactions(int n);
    CANTERA_CAPI size_t kin_nPhases(int n);
    CANTERA_CAPI size_t kin_phaseIndex(int n, const char* ph);
    CANTERA_CAPI size_t kin_reactionPhaseIndex(int n);
    CANTERA_CAPI CanteraDouble kin_reactantStoichCoeff(int n, int i, int k);
    CANTERA_CAPI CanteraDouble kin_productStoichCoeff(int n, int i, int k);
    CANTERA_CAPI int kin_getReactionType(int n, int i, size_t len, char* name);
    CANTERA_CAPI int kin_getFwdRatesOfProgress(int n, size_t len, CanteraDouble* fwdROP);
    CANTERA_CAPI int kin_getRevRatesOfProgress(int n, size_t len, CanteraDouble* revROP);
    CANTERA_CAPI int kin_getNetRatesOfProgress(int n, size_t len, CanteraDouble* netROP);
    CANTERA_CAPI int kin_getEquilibriumConstants(int n, size_t len, CanteraDouble* kc);

    CANTERA_CAPI int kin_getFwdRateConstants(int n, size_t len, CanteraDouble* kfwd);
    CANTERA_CAPI int kin_getRevRateConstants(int n, int doIrreversible, size_t len, CanteraDouble* krev);
    CANTERA_CAPI int kin_getDelta(int n, int job, size_t len, CanteraDouble* delta);
    CANTERA_CAPI int kin_getCreationRates(int n, size_t len, CanteraDouble* cdot);
    CANTERA_CAPI int kin_getDestructionRates(int n, size_t len, CanteraDouble* ddot);
    CANTERA_CAPI int kin_getNetProductionRates(int n, size_t len, CanteraDouble* wdot);
    CANTERA_CAPI int kin_getSourceTerms(int n, size_t len, CanteraDouble* ydot);
    CANTERA_CAPI CanteraDouble kin_multiplier(int n, int i);
    CANTERA_CAPI int kin_getReactionString(int n, int i, int len, char* buf);
    CANTERA_CAPI int kin_setMultiplier(int n, int i, CanteraDouble v);

    CANTERA_CAPI int kin_isReversible(int n, int i);
    CANTERA_CAPI int kin_getType(int n, size_t len, char* name);
    CANTERA_CAPI size_t kin_start(int n, int p);
    CANTERA_CAPI size_t kin_speciesIndex(int n, const char* nm, const char* ph);
    CANTERA_CAPI int kin_advanceCoverages(int n, CanteraDouble tstep);
    CANTERA_CAPI size_t kin_phase(int n, size_t i);

    CANTERA_CAPI int trans_newDefault(int th, int loglevel);
    CANTERA_CAPI int trans_new(const char* model, int th, int loglevel);
    CANTERA_CAPI int trans_del(int n);
    CANTERA_CAPI int trans_transportModel(int n, int lennm, char* nm);
    CANTERA_CAPI CanteraDouble trans_viscosity(int n);
    CANTERA_CAPI CanteraDouble trans_electricalConductivity(int n);
    CANTERA_CAPI CanteraDouble trans_thermalConductivity(int n);
    CANTERA_CAPI int trans_getThermalDiffCoeffs(int n, int ldt, CanteraDouble* dt);
    CANTERA_CAPI int trans_getMixDiffCoeffs(int n, int ld, CanteraDouble* d);
    CANTERA_CAPI int trans_getBinDiffCoeffs(int n, int ld, CanteraDouble* d);
    CANTERA_CAPI int trans_getMultiDiffCoeffs(int n, int ld, CanteraDouble* d);
    CANTERA_CAPI int trans_setParameters(int n, int type, int k, CanteraDouble* d);
    CANTERA_CAPI int trans_getMolarFluxes(int n, const CanteraDouble* state1,
                                          const CanteraDouble* state2, CanteraDouble delta, CanteraDouble* fluxes);
    CANTERA_CAPI int trans_getMassFluxes(int n, const CanteraDouble* state1,
                                         const CanteraDouble* state2, CanteraDouble delta, CanteraDouble* fluxes);

    CANTERA_CAPI int ct_getCanteraError(int buflen, char* buf);
    CANTERA_CAPI int ct_setLogWriter(void* logger);
    CANTERA_CAPI int ct_setLogCallback(LogCallback writer);
    CANTERA_CAPI int ct_addCanteraDirectory(size_t buflen, const char* buf);
    CANTERA_CAPI int ct_getDataDirectories(int buflen, char* buf, const char* sep);
    CANTERA_CAPI int ct_getCanteraVersion(int buflen, char* buf);
    CANTERA_CAPI int ct_getGitCommit(int buflen, char* buf);
    CANTERA_CAPI int ct_suppress_thermo_warnings(int suppress);
    CANTERA_CAPI int ct_use_legacy_rate_constants(int legacy);
    CANTERA_CAPI int ct_clearStorage();
    CANTERA_CAPI int ct_resetStorage();

#ifdef __cplusplus
}
#endif

#endif
