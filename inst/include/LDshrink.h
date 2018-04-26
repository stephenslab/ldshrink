#ifndef RSSR_EQTL_H
#define RSSR_EQTL_H
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>
//#include <unsupported/Eigen/CXX11/Tensor>
#include "EigenH5.h"
#include "LDshrink/LDshrink_types.h"
#include "LDshrink/LD.h"
#include "LDshrink/LD_p.h"

Eigen::MatrixXd calcLD(const c_Matrix_internal hmata,const c_arrayxd_internal mapa,const double m=85, const double Ne=11490.672741, const double cutoff=1e-3);


#endif
