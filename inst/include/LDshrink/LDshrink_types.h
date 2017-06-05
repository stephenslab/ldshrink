#ifndef RSSR_EQTL_TYPES_H
#define RSSR_EQTL_TYPES_H
#include <RcppEigen.h>


#include <RcppEigen.h>

typedef Eigen::Array<double, Eigen::Dynamic, 1> ColumnArray;
typedef Eigen::Array<double, 1,Eigen::Dynamic > RowArray;

typedef Eigen::Ref<const ColumnArray> c_column_internal;
typedef Eigen::Ref<const RowArray> c_row_internal;
typedef Eigen::Ref<const Eigen::ArrayXd> c_arrayxd_internal;
typedef Eigen::Ref<const Eigen::ArrayXXd> c_arrayxxd_internal;
typedef Eigen::Ref<const Eigen::VectorXd> c_vectorxd_internal;
typedef Eigen::MappedSparseMatrix<double> c_sparseMatrix_internal;

typedef Eigen::Ref<const Eigen::MatrixXd > c_Matrix_internal;
typedef Eigen::Ref<Eigen::MatrixXd > Matrix_internal;
typedef Eigen::Ref<Eigen::ArrayXd> arrayxd_internal;
typedef Eigen::Ref<Eigen::ArrayXXd> arrayxxd_internal;
typedef Eigen::Ref<ColumnArray> column_internal;
typedef Eigen::Ref<RowArray> row_internal;

typedef Eigen::Map<Eigen::ArrayXd> arrayxd_external;
typedef Eigen::MappedSparseMatrix<double> sparseMatrix_external;
typedef Eigen::Map<Eigen::MatrixXd > Matrix_external;
typedef Eigen::Map<Eigen::VectorXd> vectorxd_external;




#endif
