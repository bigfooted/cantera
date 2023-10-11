//! @file funcs.cpp file containing miscellaneous numerical functions.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/funcs.h"
#include "cantera/numerics/polyfit.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

CanteraDouble linearInterp(CanteraDouble x, const vector<CanteraDouble>& xpts, const vector<CanteraDouble>& fpts)
{
    if (x <= xpts[0]) {
        return fpts[0];
    }
    if (x >= xpts.back()) {
        return fpts.back();
    }
    auto loc = lower_bound(xpts.begin(), xpts.end(), x);
    int iloc = int(loc - xpts.begin()) - 1;
    CanteraDouble ff = fpts[iloc] +
                (x - xpts[iloc])*(fpts[iloc + 1]
                                  - fpts[iloc])/(xpts[iloc + 1] - xpts[iloc]);
    return ff;
}

CanteraDouble trapezoidal(const Cantera::ArrayXd& f, const Cantera::ArrayXd& x)
{
    // check length
    if (f.size() != x.size()) {
        throw CanteraError("trapezoidal",
                           "Vector lengths need to be the same.");
    }
    // Vector of f(i+1) + f(i)
    Cantera::VectorXd f_av = f.tail(f.size() - 1) + f.head(f.size() - 1);
    // Vector of x(i+1) - x(i)
    Cantera::VectorXd x_diff = x.tail(x.size() - 1) - x.head(x.size() - 1);
    // check if the coordinate is a monotonically increase vector.
    if ((x_diff.array() <= 0.0).any()) {
        throw CanteraError("trapezoidal",
            "x (coordinate) needs to be the monotonically increasing.");
    }
    return f_av.dot(x_diff) / 2.0;
}

//! Numerical integration of a function using Simpson's rule.
//! Only for odd number of points. This function is used only
//! by calling simpson.
/*!
 * Vector x contains a monotonic sequence of grid points, and
 * Vector f contains function values defined at these points.
 * The size of x and f must be the same.
 *
 * @param  f vector of function value
 * @param  x vector of function coordinate
 */
CanteraDouble basicSimpson(const Cantera::ArrayXd& f, const Cantera::ArrayXd& x)
{
    if (f.size() < 2) {
        throw CanteraError("basicSimpson",
                           "Vector lengths need to be larger than two.");
    }
    if (f.size()%2 == 0) {
        throw CanteraError("basicSimpson",
                           "Vector lengths need to be an odd number.");
    }

    size_t N = f.size() - 1;
    Cantera::VectorXd h = x.tail(N) - x.head(N);

    CanteraDouble sum = 0.0;
    for (size_t i = 1; i < N; i+=2) {
        CanteraDouble h0 = h[i-1];
        CanteraDouble h1 = h[i];
        CanteraDouble hph = h1 + h0;
        CanteraDouble hdh = h1 / h0;
        CanteraDouble hmh = h1 * h0;
        sum += (hph / 6.0) * (
                    (2.0 - hdh) * f[i - 1] + (pow(hph, 2) / hmh) * f[i] +
                    (2.0 - 1.0 / hdh) * f[i + 1]);
    }
    return sum;
}

CanteraDouble simpson(const Cantera::ArrayXd& f, const Cantera::ArrayXd& x)
{
    Cantera::ArrayXd h = x.tail(x.size() - 1) - x.head(x.size() - 1);
    if ((h <= 0.0).any()) {
        throw CanteraError("simpson",
            "Values of x need to be positive and monotonically increasing.");
    }
    if (f.size() != x.size()) {
        throw CanteraError("simpson", "Vector lengths need to be the same.");
    }

    if (f.size()%2 == 1) {
        return basicSimpson(f, x);
    } else if (f.size() == 2) {
        return 0.5 * h[0] * (f[1] + f[0]);
    } else {
        size_t N = f.size() - 1;
        // pick first N-1 point for simpson
        CanteraDouble headSimps = basicSimpson(f.head(N), x.head(N));
        // Use trapezoidal rules for the last interval
        CanteraDouble tailTrap = 0.5 * h[N-1] * (f[N] + f[N-1]);
        return headSimps + tailTrap;
    }
}

CanteraDouble numericalQuadrature(const string& method,
                           const Cantera::ArrayXd& f,
                           const Cantera::ArrayXd& x)
{
    if (method == "simpson") {
        return simpson(f, x);
    } else if (method == "trapezoidal") {
        return trapezoidal(f, x);
    } else {
        throw CanteraError("numericalQuadrature",
                           "Unknown method of numerical quadrature. "
                           "Please use 'simpson' or 'trapezoidal'");
    }
}

}
