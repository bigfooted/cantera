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
/**
 * This class provides functions that interpolate the tabulated collision integrals in
 * Monchick and Mason @cite monchick1961.
 *
 * The collision integrals computed by Monchick and Mason use the Stockmayer potential,
 * which models a polar molecule as a spherical potential with a point dipole at the
 * center). Equation 16 of Monchick and Mason @cite monchick1961 gives the potential
 * as:
 *
 * @f[
 *  \phi(r) = 4 \epsilon_0 \left[ \left(\frac{\sigma_0}{r}\right)^{12} - \left(\frac{\sigma_0}{r}\right)^6 + \delta \left(\frac{\sigma_0}{r}\right)^3 \right]
 * @f]
 *
 * Where @f$ \epsilon_0 @f$ is the depth of the potential well, @f$ \sigma_0 @f$ is the
 * distance at which the potential between the two molecules is zero, @f$ \delta @f$ is
 * defined as:
 *
 * @f[
 *  \delta = \frac{1}{4} (\mu^*)^2 \zeta \left( \theta_1, \theta_2, \phi \right)
 * @f]
 *
 * @f$ \mu^* @f$ is the reduced dipole moment. @f$ \theta_1 @f$ , @f$ \theta_2 @f$ ,
 * and @f$ \phi @f$ are angles related to trajectories of colliding molecules. In the
 * work of Monchick and Mason, these details are not what is presented in the tables.
 * Instead, the tables are presented as being functions of the reduced temperature,
 * @f$ T^* @f$, and the @f$ \delta @f$ parameter. The reduced dipole moment,
 * @f$ \mu^* @f$ is defined as:
 *
 * @f[
 *  \mu^* = \frac{\mu}{\sqrt{\epsilon_0 \sigma_0^3}}
 * @f]
 *
 * Where @f$ \mu @f$ is the dipole moment of the molecule and the other parameters
 * have been defined earlier. This work considers only the collisions of like
 * molecules, so only a single value is needed.
 *
 * The tabulated data comes from the averaged collision integrals in tables
 * V through VIII of Monchick and Mason @cite monchick1961.
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
     */
    void init(CanteraDouble tsmin, CanteraDouble tsmax);

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
    CanteraDouble quadInterp(CanteraDouble x0, CanteraDouble* x, CanteraDouble* y);

    vector<vector<CanteraDouble>> m_o22poly;
    vector<vector<CanteraDouble>> m_apoly;
    vector<vector<CanteraDouble>> m_bpoly;
    vector<vector<CanteraDouble>> m_cpoly;

    static CanteraDouble delta[8];
    static CanteraDouble tstar22[37];

    //! Table of omega22 values
    static CanteraDouble omega22_table[37*8];

    //! T* values (reduced temperature)
    static CanteraDouble tstar[39];

    //! astar table
    static CanteraDouble astar_table[39*8];

    //! bstar table
    static CanteraDouble bstar_table[39*8];

    //! cstar table
    static CanteraDouble cstar_table[39*8];

    //! Log temp
    vector<CanteraDouble> m_logTemp;

    //! Index of the tstar array that encompasses the minimum temperature
    //! fitting range value of tsmin.
    int m_nmin;

    //! Index of the tstar array that encompasses the maximum temperature
    //! fitting range value of tsmax.
    int m_nmax;
};

}
#endif
