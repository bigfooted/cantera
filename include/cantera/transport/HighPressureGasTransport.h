/**
 *  @file HighPressureGasTransport.h
 *  Interface for class HighPressureGasTransport
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_HIGHPRESSUREGASTRAN_H
#define CT_HIGHPRESSUREGASTRAN_H

// Cantera includes
#include "GasTransport.h"
#include "cantera/numerics/DenseMatrix.h"
#include "cantera/transport/MultiTransport.h"

namespace Cantera
{

//! Class MultiTransport implements transport properties for
//! high pressure gas mixtures.
/*!
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to %Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of  %Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 *
 * The implementation employs a method of corresponding states, using the Takahashi
 * @cite takahashi1975 approach for binary diffusion coefficients (using multicomponent
 * averaging rules for the mixture properties), and the Lucas method for the viscosity
 * of a high-pressure gas mixture. All methods are described in Poling et al.
 * @cite poling2001 (viscosity in Ch. 9, thermal conductivity in Ch. 10, and diffusion
 * coefficients in Ch. 11).
 *
 * @ingroup tranprops
 */
class HighPressureGasTransport : public MultiTransport
{
protected:
    //! default constructor
    /*!
     *   @param thermo  Optional parameter for the pointer to the ThermoPhase object
     */
    HighPressureGasTransport(ThermoPhase* thermo=0);

public:
    string transportModel() const override {
        return "HighPressureGas";
    }

    //! Return the thermal diffusion coefficients (kg/m/s)
    /*!
     *  Currently not implemented for this model
     */
    void getThermalDiffCoeffs(CanteraDouble* const dt) override;

    CanteraDouble thermalConductivity() override;

    /**
     * Returns the matrix of binary diffusion coefficients
     *
     *      d[ld*j +  i] = rp*m_bdiff(i,j)*(DP)_R;
     *
     * @param ld    offset of rows in the storage
     * @param d     output vector of diffusion coefficients.  Units of m**2 / s
     */
    void getBinaryDiffCoeffs(const size_t ld, CanteraDouble* const d) override;

    void getMultiDiffCoeffs(const size_t ld, CanteraDouble* const d) override;

    CanteraDouble viscosity() override;

    friend class TransportFactory;

protected:
    CanteraDouble Tcrit_i(size_t i);

    CanteraDouble Pcrit_i(size_t i);

    CanteraDouble Vcrit_i(size_t i);

    CanteraDouble Zcrit_i(size_t i);

    vector<CanteraDouble> store(size_t i, size_t nsp);

    CanteraDouble FQ_i(CanteraDouble Q, CanteraDouble Tr, CanteraDouble MW);

    CanteraDouble setPcorr(CanteraDouble Pr, CanteraDouble Tr);
};
}
#endif
