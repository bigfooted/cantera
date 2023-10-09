//! @file ChebyshevRate.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CHEBYSHEV_H
#define CT_CHEBYSHEV_H

#include "cantera/kinetics/ReactionRate.h"
#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/MultiRate.h"
#include "cantera/base/Array.h"
#include "cantera/base/global.h"

namespace Cantera
{

//! Data container holding shared data specific to ChebyshevRate
/**
 * The data container `ChebyshevData` holds precalculated data common to
 * all `ChebyshevRate` objects.
 */
struct ChebyshevData : public ReactionData
{
    ChebyshevData() = default;

    void update(CanteraDouble T) override;

    void update(CanteraDouble T, CanteraDouble P) override {
        ReactionData::update(T);
        pressure = P;
        log10P = std::log10(P);
    }

    bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    using ReactionData::update;

    //! Perturb pressure of data container
    /**
     * The method is used for the evaluation of numerical derivatives.
     * @param  deltaP  relative pressure perturbation
     */
    void perturbPressure(CanteraDouble deltaP);

    void restore() override;

    void invalidateCache() override {
        ReactionData::invalidateCache();
        pressure = NAN;
    }

    CanteraDouble pressure = NAN; //!< pressure
    CanteraDouble log10P = 0.0; //!< base 10 logarithm of pressure

protected:
    CanteraDouble m_pressure_buf = -1.0; //!< buffered pressure
};

//! Pressure-dependent rate expression where the rate coefficient is expressed
//! as a bivariate Chebyshev polynomial in temperature and pressure.
/*!
 * The rate constant can be written as:
 * @f[
 *     \log_{10} k(T,P) = \sum_{t=1}^{N_T} \sum_{p=1}^{N_P} \alpha_{tp}
 *                        \phi_t(\tilde{T}) \phi_p(\tilde{P})
 * @f]
 * where @f$ \alpha_{tp} @f$ are the constants defining the rate, @f$ \phi_n(x) @f$
 * is the Chebyshev polynomial of the first kind of degree *n* evaluated at
 * *x*, and
 * @f[
 *  \tilde{T} \equiv \frac{2T^{-1} - T_\mathrm{min}^{-1} - T_\mathrm{max}^{-1}}
 *                        {T_\mathrm{max}^{-1} - T_\mathrm{min}^{-1}}
 * @f]
 * @f[
 *  \tilde{P} \equiv \frac{2 \log_{10} P - \log_{10} P_\mathrm{min} - \log_{10} P_\mathrm{max}}
 *                        {\log_{10} P_\mathrm{max} - \log_{10} P_\mathrm{min}}
 * @f]
 * are reduced temperature and reduced pressures which map the ranges
 * @f$ (T_\mathrm{min}, T_\mathrm{max}) @f$ and
 * @f$ (P_\mathrm{min}, P_\mathrm{max}) @f$ to (-1, 1).
 *
 * A ChebyshevRate rate expression is specified in terms of the coefficient matrix
 * @f$ \alpha @f$ and the temperature and pressure ranges. Note that the
 * Chebyshev polynomials are not defined outside the interval (-1,1), and
 * therefore extrapolation of rates outside the range of temperatures and
 * pressures for which they are defined is strongly discouraged.
 *
 * @ingroup otherRateGroup
 */
class ChebyshevRate final : public ReactionRate
{
public:
    //! Default constructor.
    ChebyshevRate() = default;

    //! Constructor directly from coefficient array
    /*!
     *  @param Tmin    Minimum temperature [K]
     *  @param Tmax    Maximum temperature [K]
     *  @param Pmin    Minimum pressure [Pa]
     *  @param Pmax    Maximum pressure [Pa]
     *  @param coeffs  Coefficient array dimensioned `nT` by `nP` where `nT` and
     *      `nP` are the number of temperatures and pressures used in the fit,
     *      respectively.
     */
    ChebyshevRate(CanteraDouble Tmin, CanteraDouble Tmax, CanteraDouble Pmin, CanteraDouble Pmax,
                  const Array2D& coeffs);

    ChebyshevRate(const AnyMap& node, const UnitStack& rate_units={});

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<ChebyshevRate, ChebyshevData>>();
    }

    const string type() const override { return "Chebyshev"; }

    //! Perform object setup based on AnyMap node information
    /*!
     *  @param node  AnyMap containing rate information
     *  @param rate_units  Unit definitions specific to rate information
     */
    void setParameters(const AnyMap& node, const UnitStack& rate_units) override;

    void getParameters(AnyMap& rateNode) const override;

    //! @deprecated To be removed after %Cantera 3.0.
    void getParameters(AnyMap& rateNode, const Units& rate_units) const {
        warn_deprecated("ChebyshevRate:getParameters",
            "To be removed after Cantera 3.0. Second argument is no longer needed.");
        return getParameters(rateNode);
    }

    void validate(const string& equation, const Kinetics& kin) override;

    //! Update information specific to reaction
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    void updateFromStruct(const ChebyshevData& shared_data) {
        if (shared_data.log10P != m_log10P) {
            m_log10P = shared_data.log10P;
            CanteraDouble Pr = (2 * shared_data.log10P + PrNum_) * PrDen_;
            CanteraDouble Cnm1 = Pr;
            CanteraDouble Cn = 1;
            CanteraDouble Cnp1;
            for (size_t i = 0; i < m_coeffs.nRows(); i++) {
                dotProd_[i] = m_coeffs(i, 0);
            }
            for (size_t j = 1; j < m_coeffs.nColumns(); j++) {
                Cnp1 = 2 * Pr * Cn - Cnm1;
                for (size_t i = 0; i < m_coeffs.nRows(); i++) {
                    dotProd_[i] += Cnp1 * m_coeffs(i, j);
                }
                Cnm1 = Cn;
                Cn = Cnp1;
            }
        }
    }

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    CanteraDouble evalFromStruct(const ChebyshevData& shared_data) {
        CanteraDouble Tr = (2 * shared_data.recipT + TrNum_) * TrDen_;
        CanteraDouble Cnm1 = Tr;
        CanteraDouble Cn = 1;
        CanteraDouble Cnp1;
        CanteraDouble logk = dotProd_[0];
        for (size_t i = 1; i < m_coeffs.nRows(); i++) {
            Cnp1 = 2 * Tr * Cn - Cnm1;
            logk += Cnp1 * dotProd_[i];
            Cnm1 = Cn;
            Cn = Cnp1;
        }
        return std::pow(10, logk);
    }

    //! Set limits for ChebyshevRate object
    /*!
     *  @param Tmin    Minimum temperature [K]
     *  @param Tmax    Maximum temperature [K]
     *  @param Pmin    Minimum pressure [Pa]
     *  @param Pmax    Maximum pressure [Pa]
     */
    void setLimits(CanteraDouble Tmin, CanteraDouble Tmax, CanteraDouble Pmin, CanteraDouble Pmax);

    //! Minimum valid temperature [K]
    CanteraDouble Tmin() const {
        return Tmin_;
    }

    //! Maximum valid temperature [K]
    CanteraDouble Tmax() const {
        return Tmax_;
    }

    //! Minimum valid pressure [Pa]
    CanteraDouble Pmin() const {
        return Pmin_;
    }

    //! Maximum valid pressure [Pa]
    CanteraDouble Pmax() const {
        return Pmax_;
    }

    //! Number of points in the pressure direction
    size_t nPressure() const {
        return m_coeffs.nColumns();
    }

    //! Number of points in the temperature direction
    size_t nTemperature() const {
        return m_coeffs.nRows();
    }

    //! Access Chebyshev coefficients as 2-dimensional array with temperature and
    //! pressure dimensions corresponding to rows and columns, respectively.
    const Array2D& data() const {
        return m_coeffs;
    }

    //! Set the Chebyshev coefficients as 2-dimensional array.
    void setData(const Array2D& coeffs);

protected:
    CanteraDouble m_log10P = NAN; //!< value detecting updates
    CanteraDouble Tmin_, Tmax_; //!< valid temperature range
    CanteraDouble Pmin_, Pmax_; //!< valid pressure range
    CanteraDouble TrNum_, TrDen_; //!< terms appearing in the reduced temperature
    CanteraDouble PrNum_, PrDen_; //!< terms appearing in the reduced pressure

    Array2D m_coeffs; //!<< coefficient array
    vector<CanteraDouble> dotProd_; //!< dot product of coeffs with the reduced pressure polynomial
};

}

#endif
