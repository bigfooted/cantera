/**
 * @file ctsurf.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CTC_SURF_H
#define CTC_SURF_H

#include "clib_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

    CANTERA_CAPI int surf_setCoverages(int i, const CanteraDouble* c, int norm);
    CANTERA_CAPI int surf_getCoverages(int i, CanteraDouble* c);
    CANTERA_CAPI int surf_setConcentrations(int i, const CanteraDouble* c);
    CANTERA_CAPI int surf_getConcentrations(int i, CanteraDouble* c);
    CANTERA_CAPI int surf_setSiteDensity(int i, CanteraDouble s0);
    CANTERA_CAPI CanteraDouble surf_siteDensity(int i);
    CANTERA_CAPI int surf_setCoveragesByName(int i, const char* c);

#ifdef __cplusplus
}
#endif

#endif
