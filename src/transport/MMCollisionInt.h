/**
 * @file MMCollisionInt.h
 *  Monchick and Mason collision integrals
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MMCOLLISIONINT_H
#define CT_MMCOLLISIONINT_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

//! Calculation of Collision integrals
/*!
 * This class provides functions that interpolate the tabulated collision
 * integrals in Monchick and Mason @cite monchick1961.
 *
 * @ingroup tranprops
 */
class MMCollisionInt
{
public:
    MMCollisionInt() {}
    virtual ~MMCollisionInt() {}

    //! Initialize the object for calculation
    /*!
     *  @param tsmin       Minimum value of Tstar to carry out the fitting
     *  @param tsmax       Maximum value of Tstar to carry out the fitting
     *  @param loglevel    Set the loglevel for the object. The default
     *                     loglevel is zero, indicating no output.
     */
    void init(CanteraDouble tsmin, CanteraDouble tsmax, int loglevel = 0);

    CanteraDouble omega22(CanteraDouble ts, CanteraDouble deltastar);
    CanteraDouble astar(CanteraDouble ts, CanteraDouble deltastar);
    CanteraDouble bstar(CanteraDouble ts, CanteraDouble deltastar);
    CanteraDouble cstar(CanteraDouble ts, CanteraDouble deltastar);
    void fit(int degree, CanteraDouble deltastar, CanteraDouble* astar, CanteraDouble* bstar, CanteraDouble* cstar);
    void fit_omega22(int degree, CanteraDouble deltastar, CanteraDouble* om22);
    CanteraDouble omega11(CanteraDouble ts, CanteraDouble deltastar) {
        return omega22(ts, deltastar)/astar(ts, deltastar);
    }

private:
    CanteraDouble fitDelta(int table, int ntstar, int degree, CanteraDouble* c);

    vector<vector<CanteraDouble>> m_o22poly;
    vector<vector<CanteraDouble>> m_apoly;
    vector<vector<CanteraDouble>> m_bpoly;
    vector<vector<CanteraDouble>> m_cpoly;

    static CanteraDouble delta[8];
    static CanteraDouble tstar22[37];

    //! Table of omega22 values from MM
    static CanteraDouble omega22_table[37*8];

    //! table of tstar values
    static CanteraDouble tstar[39];

    //! astar table from MM
    static CanteraDouble astar_table[39*8];

    //! bstar table from MM
    static CanteraDouble bstar_table[39*8];

    //! cstar table from MM
    static CanteraDouble cstar_table[39*8];

    //! Log temp
    vector<CanteraDouble> m_logTemp;

    int m_nmin;
    int m_nmax;

    //! loglevel
    int m_loglevel;
};
}
#endif
