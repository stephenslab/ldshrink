#ifndef RSSR_EQTL_H
#define RSSR_EQTL_H
#include <RcppEigen.h>

#include "LDshrink/LDshrink_types.h"
#include "LDshrink/LD.h"
#include "LDshrink/LD_p.h"

Eigen::MatrixXd calcLD(const c_Matrix_internal hmata,const c_arrayxd_internal mapa,const double m=85, const double Ne=11490.672741, const double cutoff=1e-3);
  

#endif
