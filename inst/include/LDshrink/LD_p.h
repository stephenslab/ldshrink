#ifndef LD_P_H
#define LD_P_H
#include <RcppEigen.h>


#include "LDshrink_types.h"

Eigen::MatrixXd calc_cov_p(c_Matrix_external mata, c_Matrix_external matb);
Eigen::ArrayXd calc_variance(c_Matrix_external mat);
std::pair<int,int> chunk_block(const int chunk, const int chunksize, const int p);
void cov_2_cor_p(Matrix_internal covmat,arrayxd_internal rowvar,arrayxd_internal colvar);
Eigen::MatrixXd calcLD_par(const Matrix_external hmat,const arrayxd_external map,const arrayxd_external ldparams,const arrayxi_external id);
Eigen::MatrixXd calcLD_d(const c_Matrix_external hmata,const c_arrayxd_external mapa,const c_Matrix_external hmatb,const c_arrayxd_external mapb,const arrayxd_external ldparams,const bool isDiag);
#endif