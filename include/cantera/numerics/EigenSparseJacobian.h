//! @file EigenSparseJacobian.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef EIGENSPARSEJACOBIAN_H
#define EIGENSPARSEJACOBIAN_H

#include "cantera/numerics/SystemJacobian.h"
#include "cantera/numerics/eigen_sparse.h"

namespace Cantera
{

//! System Jacobians that use Eigen sparse matrices for storage
class EigenSparseJacobian : public SystemJacobian
{
public:
    EigenSparseJacobian() = default;
    void initialize(size_t networkSize) override;
    void reset() override;
    void setValue(size_t row, size_t col, CanteraDouble value) override;
    void updatePreconditioner() override;
    void updateTransient(CanteraDouble rdt, int* mask) override;

    //! Return underlying Jacobian matrix
    //! @ingroup derivGroup
    Eigen::SparseMatrix<CanteraDouble> jacobian();

    //! Return the internal preconditioner matrix
    Eigen::SparseMatrix<CanteraDouble> matrix() {
        updatePreconditioner();
        return m_matrix;
    }

    void printPreconditioner() override;

    //! Print jacobian contents
    void printJacobian();

protected:
    //! Vector of triples representing the jacobian used in preconditioning
    vector<Eigen::Triplet<CanteraDouble>> m_jac_trips;

    //! Storage of appropriately sized identity matrix for making the preconditioner
    Eigen::SparseMatrix<CanteraDouble> m_identity;

    //! Container that is the sparse preconditioner
    Eigen::SparseMatrix<CanteraDouble> m_matrix;
};

}

#endif
