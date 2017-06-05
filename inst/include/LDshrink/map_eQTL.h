#ifndef MAP_EQTL_H
#define MAP_EQTL_H
#include <RcppEigen.h>
#include "RSSReQTL_types.h"

Eigen::MatrixXd eqtl_orthogonalize_covar( c_Matrix_internal covariates);

void eqtl_orthogonalize_data(Matrix_internal data, c_Matrix_internal ortho_covar);

Eigen::MatrixXd map_eqtl_lm(Matrix_internal genotype, arrayxd_internal expression);

void fast_eqtl_lm(c_arrayxd_internal genotype, c_arrayxd_internal expression, const int n,arrayxd_internal tresid, arrayxd_internal retvec );


#endif