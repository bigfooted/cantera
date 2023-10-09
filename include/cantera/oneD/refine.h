// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REFINE_H
#define CT_REFINE_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

class Domain1D;

//! Refine Domain1D grids so that profiles satisfy adaptation tolerances
//! @ingroup onedUtilsGroup
 class Refiner
{
public:
    Refiner(Domain1D& domain);
    virtual ~Refiner() {}
    Refiner(const Refiner&) = delete;
    Refiner& operator=(const Refiner&) = delete;

    //! Set grid refinement criteria
    /*!
     *  @param ratio Maximum ratio between grid spacing at adjacent intervals.
     *      That is, `(x[j+1] - x[j]) / (x[j] - x[j-1]) < ratio`
     *  @param slope Maximum fractional change in the value of each solution
     *      component between adjacent grid points
     *  @param curve Maximum fractional change in the derivative of each
     *      solution component between adjacent grid points.
     *  @param prune Threshold for removing unnecessary grid points. `prune`
     *      should be smaller than both `slope` and `curve`. Set `prune <= 0`
     *      to disable pruning.
     */
    void setCriteria(CanteraDouble ratio = 10.0,
                     CanteraDouble slope = 0.8,
                     CanteraDouble curve = 0.8,
                     CanteraDouble prune = -0.1);

    //! Get the grid refinement criteria. @see Refiner::setCriteria
    vector<CanteraDouble> getCriteria()
    {
        return {m_ratio, m_slope, m_curve, m_prune};
    }

    void setActive(int comp, bool state = true) {
        m_active[comp] = state;
    }

    //! Set the maximum number of points allowed in the domain
    void setMaxPoints(int npmax) {
        m_npmax = npmax;
    }

    //! Returns the maximum number of points allowed in the domain
    size_t maxPoints() const {
        return m_npmax;
    }

    //! Set the minimum allowable spacing between adjacent grid points [m].
    void setGridMin(CanteraDouble gridmin) {
        m_gridmin = gridmin;
    }

    //! Returns the the minimum allowable spacing between adjacent
    //! grid points [m].
    CanteraDouble gridMin() const {
        return m_gridmin;
    }

    int analyze(size_t n, const CanteraDouble* z, const CanteraDouble* x);
    int getNewGrid(int n, const CanteraDouble* z, int nn, CanteraDouble* znew);
    int nNewPoints() {
        return static_cast<int>(m_loc.size());
    }
    void show();
    bool newPointNeeded(size_t j) {
        return m_loc.find(j) != m_loc.end();
    }
    bool keepPoint(size_t j) {
        return (m_keep[j] != -1);
    }
    CanteraDouble value(const CanteraDouble* x, size_t i, size_t j);

    CanteraDouble maxRatio() {
        return m_ratio;
    }
    CanteraDouble maxDelta() {
        return m_slope;
    }
    CanteraDouble maxSlope() {
        return m_curve;
    }
    CanteraDouble prune() {
        return m_prune;
    }

protected:
    //! Indices of grid points that need new grid points added after them
    set<size_t> m_loc;
    map<size_t, int> m_keep;
    //! Names of components that require the addition of new grid points
    set<string> m_c;
    vector<bool> m_active;
    CanteraDouble m_ratio = 10.0;
    CanteraDouble m_slope = 0.8;
    CanteraDouble m_curve = 0.8;
    CanteraDouble m_prune = -0.001;
    CanteraDouble m_min_range = 0.01;
    Domain1D* m_domain;
    size_t m_nv;
    size_t m_npmax = 1000;
    CanteraDouble m_thresh = std::sqrt(std::numeric_limits<CanteraDouble>::epsilon());
    CanteraDouble m_gridmin = 1e-10; //!< minimum grid spacing [m]
};

}

#endif
