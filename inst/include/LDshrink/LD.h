#ifndef LD_H
#define LD_H
#include <RcppEigen.h>


#include "LDshrink_types.h"


// double calc_nmsum(const double m);
// 
// Eigen::MatrixXd calc_cov( c_Matrix_internal mat);
// 
// Eigen::ArrayXd calc_variance(c_Matrix_internal mat);
// 
// void cov_2_cor(Matrix_internal covmat);
// 
// Eigen::MatrixXd cov_2_cor_exp(Matrix_external covmat);
// 
// Eigen::MatrixXd calcLD(c_Matrix_internal hmata,c_arrayxd_internal mapa,const double m, const double Ne, const double cutoff);
// 
// Eigen::SparseMatrix<double> sp_calcLD(c_Matrix_internal hmata,c_arrayxd_internal mapa,const double m, const double Ne, const double cutoff);
//   
void calcLD_pa(Eigen::MatrixXd &hmata,const std::vector<double> &mapa,Eigen::MatrixXd & S,const double m, const double Ne, const double cutoff);
  
#endif