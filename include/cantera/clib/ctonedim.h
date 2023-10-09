/**
 * @file ctonedim.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CTC_ONEDIM_H
#define CTC_ONEDIM_H

#include "clib_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

    CANTERA_CAPI int ct_clearOneDim();
    CANTERA_CAPI int domain_new(const char* type, int i, const char* id);
    CANTERA_CAPI int domain_del(int i);
    CANTERA_CAPI int domain_type(int i);
    CANTERA_CAPI int domain_type3(int i, size_t lennm, char* nm);
    CANTERA_CAPI size_t domain_index(int i);
    CANTERA_CAPI size_t domain_nComponents(int i);
    CANTERA_CAPI size_t domain_nPoints(int i);
    CANTERA_CAPI int domain_componentName(int i, int n, int sz, char* nameout);
    CANTERA_CAPI size_t domain_componentIndex(int i, const char* name);
    CANTERA_CAPI int domain_setBounds(int i, int n, CanteraDouble lower,
                                      CanteraDouble upper);
    CANTERA_CAPI CanteraDouble domain_lowerBound(int i, int n);
    CANTERA_CAPI CanteraDouble domain_upperBound(int i, int n);
    CANTERA_CAPI int domain_setSteadyTolerances(int i, int n, CanteraDouble rtol,
                                                CanteraDouble atol);
    CANTERA_CAPI int domain_setTransientTolerances(int i, int n, CanteraDouble rtol,
                                                   CanteraDouble atol);
    CANTERA_CAPI CanteraDouble domain_rtol(int i, int n);
    CANTERA_CAPI CanteraDouble domain_atol(int i, int n);
    CANTERA_CAPI int domain_setupGrid(int i, size_t npts, const CanteraDouble* grid);
    CANTERA_CAPI int domain_setID(int i, const char* id);
    CANTERA_CAPI CanteraDouble domain_grid(int i, int n);

    CANTERA_CAPI int bdry_setMdot(int i, CanteraDouble mdot);
    CANTERA_CAPI int bdry_setTemperature(int i, CanteraDouble t);
    CANTERA_CAPI int bdry_setSpreadRate(int i, CanteraDouble v);
    CANTERA_CAPI int bdry_setMoleFractions(int i, const char* x);
    CANTERA_CAPI CanteraDouble bdry_temperature(int i);
    CANTERA_CAPI CanteraDouble bdry_spreadRate(int i);
    CANTERA_CAPI CanteraDouble bdry_massFraction(int i, int k);
    CANTERA_CAPI CanteraDouble bdry_mdot(int i);

    CANTERA_CAPI int reactingsurf_setkineticsmgr(int i, int j);
    CANTERA_CAPI int reactingsurf_enableCoverageEqs(int i, int onoff);

    CANTERA_CAPI int inlet_new();
    CANTERA_CAPI int outlet_new();
    CANTERA_CAPI int outletres_new();
    CANTERA_CAPI int symm_new();
    CANTERA_CAPI int surf_new();
    CANTERA_CAPI int reactingsurf_new();

    CANTERA_CAPI int inlet_setSpreadRate(int i, CanteraDouble v);

    CANTERA_CAPI int stflow_new(int iph, int ikin, int itr, int itype);
    CANTERA_CAPI int stflow_setTransport(int i, int itr);
    CANTERA_CAPI int stflow_enableSoret(int i, int iSoret);
    CANTERA_CAPI int stflow_setPressure(int i, CanteraDouble p);
    CANTERA_CAPI CanteraDouble stflow_pressure(int i);
    CANTERA_CAPI int stflow_setFixedTempProfile(int i, size_t n, const CanteraDouble* pos,
            size_t m, const CanteraDouble* temp);
    CANTERA_CAPI int stflow_solveEnergyEqn(int i, int flag);

    CANTERA_CAPI int sim1D_new(size_t nd, const int* domains);
    CANTERA_CAPI int sim1D_del(int i);
    CANTERA_CAPI int sim1D_setValue(int i, int dom, int comp, int localPoint, CanteraDouble value);
    CANTERA_CAPI int sim1D_setProfile(int i, int dom, int comp,
                                      size_t np, const CanteraDouble* pos, size_t nv, const CanteraDouble* v);
    CANTERA_CAPI int sim1D_setFlatProfile(int i, int dom, int comp, CanteraDouble v);
    CANTERA_CAPI int sim1D_show(int i, const char* fname);
    CANTERA_CAPI int sim1D_showSolution(int i, const char* fname);
    CANTERA_CAPI int sim1D_setTimeStep(int i, CanteraDouble stepsize, size_t ns, const int* nsteps);
    CANTERA_CAPI int sim1D_getInitialSoln(int i);
    CANTERA_CAPI int sim1D_solve(int i, int loglevel, int refine_grid);
    CANTERA_CAPI int sim1D_refine(int i, int loglevel);
    CANTERA_CAPI int sim1D_setRefineCriteria(int i, int dom, CanteraDouble ratio,
            CanteraDouble slope, CanteraDouble curve, CanteraDouble prune);
    CANTERA_CAPI int sim1D_setGridMin(int i, int dom, CanteraDouble gridmin);
    CANTERA_CAPI int sim1D_save(int i, const char* fname, const char* id,
                                const char* desc);
    CANTERA_CAPI int sim1D_restore(int i, const char* fname, const char* id);
    CANTERA_CAPI int sim1D_writeStats(int i, int printTime);
    CANTERA_CAPI int sim1D_domainIndex(int i, const char* name);
    CANTERA_CAPI CanteraDouble sim1D_value(int i, int idom, int icomp, int localPoint);
    CANTERA_CAPI CanteraDouble sim1D_workValue(int i, int idom,
                                        int icomp, int localPoint);
    CANTERA_CAPI int sim1D_eval(int i, CanteraDouble rdt, int count);
    CANTERA_CAPI int sim1D_setMaxJacAge(int i, int ss_age, int ts_age);
    CANTERA_CAPI int sim1D_setFixedTemperature(int i, CanteraDouble temp);

#ifdef __cplusplus
}
#endif

#endif
