//! @file polyfit.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/polyfit.h"
#include "cantera/numerics/eigen_dense.h"
#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

CanteraDouble polyfit(size_t n, size_t deg, const CanteraDouble* xp, const CanteraDouble* yp,
               const CanteraDouble* wp, CanteraDouble* pp)
{
    ConstMappedVector x(xp, n);
    Cantera::VectorXd y = ConstMappedVector(yp, n);
    MappedVector p(pp, deg+1);

    if (deg >= n) {
        throw CanteraError("polyfit", "Polynomial degree ({}) must be less "
            "than number of input data points ({})", deg, n);
    }

    // Construct A such that each row i of A has the elements
    // 1, x[i], x[i]^2, x[i]^3 ... + x[i]^deg
    Cantera::MatrixXd A(n, deg+1);
    A.col(0).setConstant(1.0);

    if (deg > 0) {
        A.col(1) = x;
    }
    for (size_t i = 1; i < deg; i++) {
        A.col(i+1) = A.col(i).array() * x.array();
    }

    if (wp != nullptr && wp[0] > 0) {
        // For compatibility with old Fortran dpolft, input weights are the
        // squares of the weight vector used in this algorithm
        Cantera::VectorXd w = ConstMappedVector(wp, n).cwiseSqrt().eval();

        // Multiply by the weights on both sides
        A = w.asDiagonal() * A;
        y.array() *= w.array();
    }

    // Solve W*A*p = W*y to find the polynomial coefficients
    p = A.colPivHouseholderQr().solve(y);

    // Evaluate the computed polynomial at the input x coordinates to compute
    // the RMS error as the return value
    return (A*p - y).eval().norm() / sqrt(n);
}

}
