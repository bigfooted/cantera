//! @file ResidJacEval.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#define CT_SKIP_DEPRECATION_WARNINGS
#include "cantera/numerics/ResidJacEval.h"
#include "cantera/base/global.h"

namespace Cantera
{
ResidJacEval::ResidJacEval(CanteraDouble atol) :
    m_atol(atol)
{
    warn_deprecated("class ResidJacEval", "To be removed after Cantera 3.0");
}

int ResidJacEval::nEquations() const
{
    return neq_;
}

void ResidJacEval::setAtol(CanteraDouble atol)
{
    m_atol = atol;
    if (m_atol <= 0.0) {
        throw CanteraError("ResidJacEval::setAtol",
                           "atol must be greater than zero");
    }
}

int ResidJacEval::getInitialConditions(CanteraDouble t0, CanteraDouble* const y,
                                       CanteraDouble* const ydot)
{
    for (int i = 0; i < neq_; i++) {
        y[i] = 0.0;
    }
    if (ydot) {
        for (int i = 0; i < neq_; i++) {
            ydot[i] = 0.0;
        }
    }
    return 1;
}

void ResidJacEval::user_out2(const int ifunc, const CanteraDouble t,
                             const CanteraDouble deltaT, const CanteraDouble* y,
                             const CanteraDouble* ydot)
{
}

void ResidJacEval::user_out(const int ifunc, const CanteraDouble t,
                            const CanteraDouble* y, const CanteraDouble* ydot)
{
    user_out2(ifunc, t, 0.0, y, ydot);
}

int ResidJacEval::evalTimeTrackingEqns(const CanteraDouble t,
                                       const CanteraDouble delta_t,
                                       const CanteraDouble* y,
                                       const CanteraDouble* ydot)
{
    return 1;
}

int ResidJacEval::calcDeltaSolnVariables(const CanteraDouble t,
                                         const CanteraDouble* const ySoln,
                                         const CanteraDouble* const ySolnDot,
                                         CanteraDouble* const deltaYSoln,
                                         const CanteraDouble* const solnWeights)
{
    if (!solnWeights) {
        for (int i = 0; i < neq_; i++) {
            deltaYSoln[i] = m_atol + fabs(1.0E-6 * ySoln[i]);
        }
    } else {
        for (int i = 0; i < neq_; i++) {
            deltaYSoln[i] = std::max(1.0E-2 * solnWeights[i], 1.0E-6 * fabs(ySoln[i]));
        }
    }
    return 1;
}

void ResidJacEval::calcSolnScales(const CanteraDouble t,
                                  const CanteraDouble* const ysoln,
                                  const CanteraDouble* const ysolnOld,
                                  CanteraDouble* const ysolnScales)
{
    if (ysolnScales && ysolnScales[0] == 0.0) {
        for (int i = 0; i < neq_; i++) {
            ysolnScales[i] = 1.0;
        }
    }
}

CanteraDouble ResidJacEval::filterNewStep(CanteraDouble t, const CanteraDouble* const ybase, CanteraDouble* const step)
{
    return 0.0;
}

CanteraDouble ResidJacEval::filterSolnPrediction(CanteraDouble t, CanteraDouble* const y)
{
    return 0.0;
}

bool ResidJacEval::evalStoppingCritera(const CanteraDouble t,
                                       const CanteraDouble delta_t,
                                       const CanteraDouble* const y,
                                       const CanteraDouble* const ydot)
{
    return false;
}

int ResidJacEval::matrixConditioning(CanteraDouble* const matrix, const int nrows,
                                     CanteraDouble* const rhs)
{
    return 1;
}

int ResidJacEval::evalResidNJ(const CanteraDouble t, const CanteraDouble deltaT,
                              const CanteraDouble* y, const CanteraDouble* ydot,
                              CanteraDouble* const resid,
                              const ResidEval_Type_Enum evalType,
                              const int id_x, const CanteraDouble delta_x)
{
    throw NotImplementedError("ResidJacEval::evalResidNJ");
}

int ResidJacEval::eval(const CanteraDouble t, const CanteraDouble* const y, const CanteraDouble* const ydot,
                       CanteraDouble* const r)
{
    CanteraDouble deltaT = -1.0;
    return evalResidNJ(t, deltaT, y, ydot, r);
}

int ResidJacEval::evalJacobian(const CanteraDouble t, const CanteraDouble delta_t,
                               CanteraDouble cj, const CanteraDouble* const y,
                               const CanteraDouble* const ydot, DenseMatrix& J,
                               CanteraDouble* const resid)
{
    CanteraDouble* const* jac_colPts = J.colPts();
    return evalJacobianDP(t, delta_t, cj, y, ydot, jac_colPts, resid);
}

int ResidJacEval::evalJacobianDP(const CanteraDouble t, const CanteraDouble delta_t,
                                 const CanteraDouble c_j,
                                 const CanteraDouble* const y,
                                 const CanteraDouble* const ydot,
                                 CanteraDouble* const* jac_colPts,
                                 CanteraDouble* const resid)
{
    throw NotImplementedError("ResidJacEval::evalJacobianDP");
}

}
