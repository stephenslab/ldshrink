#pragma once
#include <RcppEigen.h>

template <class ArgType> struct mapdiff_helper {

  typedef Eigen::Matrix<typename ArgType::Scalar, ArgType::SizeAtCompileTime,
                        ArgType::SizeAtCompileTime, Eigen::ColMajor,
                        ArgType::MaxSizeAtCompileTime,
                        ArgType::MaxSizeAtCompileTime>
      MatrixType;
};

template <class ArgType, typename Scalar = typename ArgType::Scalar>
class mapdiff_functor {
  const Eigen::Map<const ArgType> m_vec;
  const Scalar m;
  const Scalar Ne;
  const Scalar cutoff;

public:
  mapdiff_functor(const Eigen::Map<const ArgType> arg, const Scalar m_,
                  const Scalar Ne_, const Scalar cutoff_)
      : m_vec(arg), m(m_), Ne(Ne_), cutoff(cutoff_) {
    //    Rcpp::Rcerr<<"m:"<<m<<",Ne:"<<Ne<<"cutoff:"<<cutoff<<std::endl;
  }
  const Scalar operator()(Eigen::Index row, Eigen::Index col) const {

    using namespace Eigen;
    const Scalar map_dist = std::fabs(m_vec(row) - m_vec(col));
    Scalar rho = 4 * Ne * map_dist / 100;
    rho = -rho / (2 * m);
    rho = std::exp(rho);

    //   true_map_matrix(i,j)=
    // auto rho = std::exp(-(4*Ne*std::fabs(m_vec(row)-m_vec(col))/100)/(2*m));
    return (rho < cutoff ? 0 : rho);
  }
};

template <class ArgType, typename Scalar = typename ArgType::Scalar>
Eigen::CwiseNullaryOp<mapdiff_functor<ArgType>,
                      typename mapdiff_helper<ArgType>::MatrixType>
makeMapDiff(const Eigen::Map<const ArgType> &arg, const Scalar m,
            const Scalar Ne, const Scalar cutoff) {
  typedef typename mapdiff_helper<ArgType>::MatrixType MatrixType;
  return MatrixType::NullaryExpr(
      arg.size(), arg.size(),
      mapdiff_functor<ArgType>(arg.derived(), m, Ne, cutoff));
}
