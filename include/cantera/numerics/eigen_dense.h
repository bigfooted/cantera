// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_EIGEN_DENSE_H
#define CT_EIGEN_DENSE_H

#include "cantera/base/config.h"
#if CT_USE_SYSTEM_EIGEN
    #if CT_USE_SYSTEM_EIGEN_PREFIXED
    #include <eigen3/Eigen/Dense>
    #else
    #include <Eigen/Dense>
    #endif
#else
#include "cantera/ext/Eigen/Dense"
#endif

namespace Cantera
{

//! @addtogroup matrices
//! @{

//typedef Eigen::Map<Eigen::MatrixXd> MappedMatrix;
//typedef Eigen::Map<const Eigen::MatrixXd> ConstMappedMatrix;
//typedef Eigen::Map<Eigen::VectorXd> MappedVector;
//typedef Eigen::Map<const Eigen::VectorXd> ConstMappedVector;
//typedef Eigen::Map<Eigen::RowVectorXd> MappedRowVector;
//typedef Eigen::Map<const Eigen::RowVectorXd> ConstMappedRowVector;

typedef Eigen::Map<Eigen::Matrix<CanteraDouble,Eigen::Dynamic,Eigen::Dynamic>> MappedMatrix;
typedef Eigen::Map<const Eigen::Matrix<CanteraDouble,Eigen::Dynamic,Eigen::Dynamic>> ConstMappedMatrix;

typedef Eigen::Map<Eigen::Matrix<CanteraDouble,Eigen::Dynamic,1>> MappedVector;
typedef Eigen::Map<const Eigen::Matrix<CanteraDouble,Eigen::Dynamic,1>> ConstMappedVector;

typedef Eigen::Map<Eigen::Matrix<CanteraDouble,1, Eigen::Dynamic>> MappedRowVector;
typedef Eigen::Map<const Eigen::Matrix<CanteraDouble,1, Eigen::Dynamic>> ConstMappedRowVector;

typedef Eigen::Matrix<CanteraDouble,Eigen::Dynamic,1> VectorXd;
typedef Eigen::Matrix<CanteraDoublePassive,Eigen::Dynamic,1> PassiveVectorXd;

typedef Eigen::Array<CanteraDouble,Eigen::Dynamic,1> ArrayXd;
typedef Eigen::Array<CanteraDoublePassive,Eigen::Dynamic,1> PassiveArrayXd;

typedef Eigen::Matrix<CanteraDouble,Eigen::Dynamic,Eigen::Dynamic> MatrixXd;

//! @}

}

#endif
