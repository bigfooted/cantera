//! @file MultiNewton.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTINEWTON_H
#define CT_MULTINEWTON_H

#include "MultiJac.h"

namespace Cantera
{

//! @defgroup onedUtilsGroup Utilities
//! Utility classes and functions for one-dimensional problems.
//! @ingroup onedGroup

/**
 * Newton iterator for multi-domain, one-dimensional problems.
 * Used by class OneDim.
 * @ingroup onedUtilsGroup
 */
class MultiNewton
{
public:
    MultiNewton(int sz);
    virtual ~MultiNewton() {};
    MultiNewton(const MultiNewton&) = delete;
    MultiNewton& operator=(const MultiNewton&) = delete;

    size_t size() {
        return m_n;
    }

    //! Compute the undamped Newton step.  The residual function is evaluated
    //! at `x`, but the Jacobian is not recomputed.
    void step(CanteraDouble* x, CanteraDouble* step, OneDim& r, MultiJac& jac, int loglevel);

    /**
     * Return the factor by which the undamped Newton step 'step0'
     * must be multiplied in order to keep all solution components in
     * all domains between their specified lower and upper bounds.
     */
    CanteraDouble boundStep(const CanteraDouble* x0, const CanteraDouble* step0,
                     const OneDim& r, int loglevel);

    /**
     * On entry, step0 must contain an undamped Newton step for the solution x0.
     * This method attempts to find a damping coefficient such that the next
     * undamped step would have a norm smaller than that of step0. If
     * successful, the new solution after taking the damped step is returned in
     * x1, and the undamped step at x1 is returned in step1.
     */
    int dampStep(const CanteraDouble* x0, const CanteraDouble* step0, CanteraDouble* x1, CanteraDouble* step1,
                 CanteraDouble& s1, OneDim& r, MultiJac& jac, int loglevel, bool writetitle);

    //! Compute the weighted 2-norm of `step`.
    CanteraDouble norm2(const CanteraDouble* x, const CanteraDouble* step, OneDim& r) const;

    /**
     * Find the solution to F(X) = 0 by damped Newton iteration. On entry, x0
     * contains an initial estimate of the solution. On successful return, x1
     * contains the converged solution.
     */
    int solve(CanteraDouble* x0, CanteraDouble* x1, OneDim& r, MultiJac& jac, int loglevel);

    //! Set options.
    void setOptions(int maxJacAge = 5) {
        m_maxAge = maxJacAge;
    }

    //! Change the problem size.
    void resize(size_t points);

protected:
    //! Work arrays of size #m_n used in solve().
    vector<CanteraDouble> m_x, m_stp, m_stp1;

    int m_maxAge = 5;

    //! number of variables
    size_t m_n;

    CanteraDouble m_elapsed = 0.0;
};
}

#endif
