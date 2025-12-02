// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_EIGEN_SPARSE_H
#define CT_EIGEN_SPARSE_H

#include "cantera/base/config.h"
#include "cantera/base/ct_typedefs.h"
#if CT_USE_SYSTEM_EIGEN
    #if CT_USE_SYSTEM_EIGEN_PREFIXED
    #include <eigen3/Eigen/Sparse>
    #else
    #include <Eigen/Sparse>
    #endif
#else
#include "cantera/ext/Eigen/Sparse"
#endif

namespace Cantera
{
//! @ingroup matrices
typedef std::vector<Eigen::Triplet<CanteraDouble>> SparseTriplets;

typedef Eigen::Matrix<CanteraDouble,Eigen::Dynamic,1> VectorXd;
typedef Eigen::Matrix<CanteraDoublePassive,Eigen::Dynamic,1> PassiveVectorXd;

typedef Eigen::Array<CanteraDouble,Eigen::Dynamic,1> ArrayXd;
typedef Eigen::Array<CanteraDoublePassive,Eigen::Dynamic,1> PassiveArrayXd;

typedef Eigen::Matrix<CanteraDouble,Eigen::Dynamic,Eigen::Dynamic> MatrixXd;

}

#endif
