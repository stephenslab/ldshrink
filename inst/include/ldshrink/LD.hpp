#pragma once
#include <RcppEigen.h>
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
//#include <memory>
#include "mapdiff.hpp"
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppParallel)]]

#include <RcppParallel.h>


inline double calc_nmsum(const double m) {
  int msize = (2 * (int)m - 1);
  Eigen::ArrayXd tx(msize);
  tx.setLinSpaced(msize, 1, (int)(2 * m - 1));
  return (1 / tx).sum();
}

inline double calc_theta(const double m){
  double nmsum=calc_nmsum(m);
  return((1/nmsum)/(2*m+1/nmsum));
}

template <typename T, int RN> class ldshrink_data {
public:
  std::vector<T> map;
  size_t p;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, RN> data;
  std::vector<T> data_vars;
  const bool SNPfirst;
  const size_t N;
  const bool isGeno = true;
  std::vector<size_t> indices;

  ldshrink_data(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, RN> data_,
                std::vector<T> map_, std::vector<size_t> indices_)
      : map(std::move(map_)), p(map.size()), data(std::move(data_)),
        SNPfirst(data.rows() == p), is_scaled(false), is_vars(false),
        N(SNPfirst ? data.cols() : data.rows()), indices(std::move(indices_)) {
    if (SNPfirst) {
      data.transposeInPlace();
    }
    if (indices.empty()) {
      indices.resize(p);
      std::iota(indices.begin(), indices.end(), 0);
    }
    if (indices.size() != p) {
      Rcpp::stop(
          "indices must be equal in size to the length of the genetic map");
    }

    if (p == N) {
      Rcpp::Rcerr << "Warning: genoype_data matrix is of dimension " << p << "x"
                  << p << ", assuming data is stored SNPxSample!";
    }
    if (!std::is_sorted(map.begin(), map.end(), std::less<T>())) {
      Rcpp::stop("Recombination map must be non-decreasing\n");
    }
  }
  bool is_scaled;
  bool is_vars;
  void scale_data() {
    if (!is_scaled) {
      data = data.rowwise() - data.colwise().mean();
      is_scaled = true;
    }
  }
  void calc_vars() {
    if (!is_vars) {
      if (!is_scaled) {
        scale_data();
      }
      data_vars.resize(p);
      tbb::parallel_for(tbb::blocked_range<size_t>(0, p),
                        [&](const tbb::blocked_range<size_t> &r) {
                          for (size_t i = r.begin(); i != r.end(); i++) {
                            data_vars[i] =
                                (data.col(i).array().square().sum() / (N - 1));
                          }
                        });
      is_vars = true;
    }
  }
};

template <typename T, int RN> class Sparse_cov {
  tbb::concurrent_vector<
      typename Eigen::Triplet<T, typename Eigen::SparseMatrix<T>::StorageIndex>>
      tripletList;
  T m = 0;
  T Ne = 0;
  T cutoff = 0;
  T r2cutoff = 0;
  bool progress = true;
  typedef Eigen::Triplet<T, typename Eigen::SparseMatrix<T>::StorageIndex> Trip;

public:
  Sparse_cov(bool useldshrink_ = true) { tripletList.reserve(10000); }
  Sparse_cov(const T m_, const T Ne_, const T cutoff_, const T r2cutoff_ = 0,
             const bool progress_ = true)
      : m(m_), Ne(Ne_), cutoff(cutoff_), r2cutoff(r2cutoff_),
        progress(progress_) {
    tripletList.reserve(10000);
  }
  std::pair<size_t, size_t> row_col_t(const size_t k, const size_t p) {
    size_t i = floor((2 * p + 1 - sqrt((2 * p + 1) * (2 * p + 1) - 8 * k)) / 2);
    size_t j = k - (2 * p - 1 - i) * i / 2;
    return (std::pair<size_t, size_t>(i, j));
  }
  void ldshrink(const ldshrink_data<T, RN> *a) {
    const size_t p = a->p;
    const size_t p_od = (((p * p) - p) / 2) + p;
    Progress pr(p_od, progress);
    if (!a->is_vars) {
      Rcpp::stop("data must be scaled before calling ldshrink");
    }
    const auto &map = a->map;
    const auto &scaled_data = a->data;
    const auto &data_vars = a->data_vars;
    const auto &indices = a->indices;
    const int num_ind = a->N;
    const bool isGeno = a->isGeno;
    const T GenoMult = isGeno ? 0.5 : 1;
    const T theta = calc_theta(m);
    const T pre_mult =
        GenoMult * (1 - theta) * (1 - theta) + 0.5 * theta * (1 - 0.5 * theta);
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, p_od),
        [&](const tbb::blocked_range<size_t> &r) {
          for (size_t ir = r.begin(); ir != r.end(); ir++) {
            auto idx = row_col_t(ir, p);
            const size_t i = idx.first;
            const size_t j = idx.second;
            const size_t index_i = indices[i];
            const size_t index_j = indices[j];
            if (i == j) {
              tripletList.push_back(Trip(index_i, index_j, 1.0));
            } else {
              T map_dist = map[j] - map[i];
              T rho = 4 * Ne * (map_dist) / 100;
              rho = -rho / (2 * m);
              rho = std::exp(rho);
              if (rho >= cutoff) {
                T cov = (scaled_data.col(j).dot(scaled_data.col(i)) /
                         (num_ind - 1)) *
                        rho * GenoMult * (1 - theta) * (1 - theta);
                cov = cov / (pre_mult * std::sqrt(data_vars[i] * data_vars[j]));
                if (((cov * cov) >= r2cutoff)) {
                  tripletList.push_back(Trip(index_i, index_j, cov));
                }
              }
            }
          }
          pr.increment(r.end() - r.begin());
        });
  }
  void cor(const ldshrink_data<T, RN> *a) {
    const size_t p = a->p;
    const size_t p_od = (((p * p) - p) / 2) + p;
    Progress pr(p_od, progress);
    if (!a->is_vars) {
      Rcpp::stop("data must be scaled before calling ldshrink");
    }
    const auto &scaled_data = a->data;
    const auto &data_vars = a->data_vars;
    const auto &indices = a->indices;
    const int num_ind = a->N;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, p_od),
        [&](const tbb::blocked_range<size_t> &r) {
          for (size_t ir = r.begin(); ir != r.end(); ir++) {
            auto idx = row_col_t(ir, p);
            const size_t i = idx.first;
            const size_t j = idx.second;
            if (i == j) {
              tripletList.push_back(Trip(indices[i], indices[i], 1.0));
            } else {
              T cov =
                  (scaled_data.col(j).dot(scaled_data.col(i)) / (num_ind - 1)) /
                  std::sqrt(data_vars[i] * data_vars[j]);
              if (((cov * cov) >= r2cutoff)) {
                tripletList.push_back(Trip(indices[i], indices[j], cov));
              }
            }
          }
          pr.increment(r.end() - r.begin());
        });
  }

  void ldshrink(const ldshrink_data<T, RN> *a, const ldshrink_data<T, RN> *b) {
    const size_t p_a = a->p;
    const size_t p_b = b->p;

    const size_t p_od = p_a * p_b;
    Progress pr(p_od, progress);
    if ((!a->is_vars) || (!b->is_vars)) {
      Rcpp::stop("data must be scaled before calling ldshrink");
    }
    const auto &map_a = a->map;
    const auto &scaled_data_a = a->data;
    const auto &data_vars_a = a->data_vars;
    const auto &index_a = a->indices;

    const auto &map_b = b->map;
    const auto &scaled_data_b = b->data;
    const auto &data_vars_b = b->data_vars;
    const auto &index_b = b->indices;
    const int num_ind = a->N;
    const bool isGeno = a->isGeno;
    const T GenoMult = isGeno ? 0.5 : 1;
    const T theta = calc_theta(m);
    const T pre_mult =
        GenoMult * (1 - theta) * (1 - theta) + 0.5 * theta * (1 - 0.5 * theta);

    tbb::parallel_for(
        tbb::blocked_range2d<size_t>(0, p_a, 0, p_b),
        [&](const tbb::blocked_range2d<size_t> &r) {
          for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
            for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
              T map_dist = map_b[j] - map_a[i];
              T rho = 4 * Ne * (map_dist) / 100;
              rho = -rho / (2 * m);
              rho = std::exp(rho);
              if (rho >= cutoff) {
                T cov = (scaled_data_a.col(i).dot(scaled_data_b.col(j)) /
                         (num_ind - 1)) *
                        rho * GenoMult * (1 - theta) * (1 - theta);
                cov = cov /
                      (pre_mult * std::sqrt(data_vars_a[i] * data_vars_b[j]));
                if (((cov * cov) >= r2cutoff) && (i != j)) {
                  tripletList.push_back(Trip(index_a[i], index_b[j], cov));
                  //			   tripletList.push_back(Trip(i,j,cov));
                }
              }
            }
          }
          size_t tot_prog = (r.rows().end() - r.rows().begin()) *
                            (r.cols().end() - r.cols().begin());
          pr.increment(tot_prog);
        });
  }

  void cor(const ldshrink_data<T, RN> *a, const ldshrink_data<T, RN> *b) {
    const size_t p_a = a->p;
    const size_t p_b = b->p;
    const size_t p_od = p_a * p_b;
    Progress pr(p_od, progress);
    if ((!a->is_vars) || (!b->is_vars)) {
      Rcpp::stop("data must be scaled before calling ldshrink");
    }
    const auto &scaled_data_a = a->data;
    const auto &data_vars_a = a->data_vars;
    const auto &index_a = a->indices;
    const auto &scaled_data_b = b->data;
    const auto &data_vars_b = b->data_vars;
    const auto &index_b = b->indices;
    const int num_ind = a->N;
    tbb::parallel_for(
        tbb::blocked_range2d<size_t>(0, p_a, 0, p_b),
        [&](const tbb::blocked_range2d<size_t> &r) {
          for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
            for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
              T cov = (scaled_data_a.col(i).dot(scaled_data_b.col(j)) /
                       (num_ind - 1)) /
                      std::sqrt(data_vars_a[i] * data_vars_b[j]);
              if (((cov * cov) >= r2cutoff) && (i != j)) {
                tripletList.push_back(Trip(i, j, cov));
              }
            }
          }

          size_t tot_prog = (r.rows().end() - r.rows().begin()) *
                            (r.cols().end() - r.cols().begin());
          pr.increment(tot_prog);
        });
  }
  Eigen::SparseMatrix<T> sparseMatrix(const size_t p_a,
                                      const size_t p_b) const {
    const size_t num_elem = tripletList.size();
    Eigen::SparseMatrix<T> object(p_a, p_b);
    Rcpp::Rcerr << "About to create sparse matrix of size: ";
    Rcpp::Rcerr << num_elem << std::endl;
    object.setFromTriplets(tripletList.begin(), tripletList.end());
    Rcpp::Rcerr << "Sparse Matrix created" << std::endl;
    return (object);
  }
  Rcpp::S4 dsCMatrix(const size_t p_a, const size_t p_b) const {
    auto object = this->sparseMatrix(p_a, p_b);
    using namespace Rcpp;
    const int nnz = object.nonZeros();
    Rcpp::Rcerr << "Number of entries: " << nnz << std::endl;
    S4 ans("dsCMatrix");
    ans.slot("Dim") = Dimension(object.rows(), object.cols());
    ans.slot("i") =
        IntegerVector(object.innerIndexPtr(), object.innerIndexPtr() + nnz);
    ans.slot("p") =
        IntegerVector(object.outerIndexPtr(),
                      object.outerIndexPtr() + object.outerSize() + 1);
    ans.slot("x") =
        Rcpp::NumericVector(object.valuePtr(), object.valuePtr() + nnz);
    return (ans);
  }
  template <int RTYPE>
  Rcpp::DataFrame ld2df(const Rcpp::Vector<RTYPE> rsid_a,
                        const Rcpp::Vector<RTYPE> rsid_b) const {
    const size_t num_trip = tripletList.size();
    Rcpp::Vector<RTYPE> rowid(num_trip);
    Rcpp::Vector<RTYPE> colid(num_trip);
    Rcpp::NumericVector rvec(num_trip);
    for (size_t i = 0; i < num_trip; i++) {
      auto trip = tripletList[i];
      rowid(i) = rsid_a(trip.row());
      colid(i) = rsid_b(trip.col());
      rvec(i) = trip.value();
    }
    using namespace Rcpp;
    const bool stringsAsFactors = false;
    return (Rcpp::DataFrame::create(_["rowsnp"] = rowid, _["colsnp"] = colid,
                                    _["r"] = rvec,
                                    _["stringsAsFactors"] = stringsAsFactors));
  }
  Rcpp::DataFrame ld2df(const size_t p_a, const size_t p_b) const {
    Rcpp::IntegerVector indices_a(p_a);
    std::iota(indices_a.begin(), indices_a.end(), 1);
    Rcpp::IntegerVector indices_b(p_b);
    std::iota(indices_b.begin(), indices_b.end(), 1);
    return (ld2df(indices_a, indices_b));
  }
};

template <typename T> class ldshrinker {
  std::vector<T> Sdat;

public:
  const T m;
  const T theta;
  const T cutoff;
  const T Ne;
  Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> S;
  Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> mapd;
  //  const size_t N;
  ldshrinker(const T m_, const T Ne_, const T cutoff_)
      : m(m_), theta(calc_theta(m)), cutoff(cutoff_), Ne(Ne_), mapd(nullptr, 1),
        S(nullptr, 1, 1) {}

  ldshrinker(Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> &S_,
             Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> &map_,
             const T m_, const T Ne_, const T cutoff_,
             Rcpp::NumericMatrix genotype_data, const bool isGeno = true)
      : m(m_), theta(calc_theta(m)), cutoff(cutoff_), Ne(Ne_),
        mapd(map_.data(), map_.size()), S(S_.data(), S_.rows(), S_.cols()) {

    const size_t p = mapd.size();
    const size_t Srows = S.rows();
    if (S.cols() != p || Srows != p) {
      Rcpp::stop("Covariance Matrix have column(row) number equal to the size "
                 "of genetic map");
    }
    const size_t genrows = genotype_data.rows();
    const bool SNPfirst = (genrows == p);
    const size_t gencols = genotype_data.cols();
    if (SNPfirst && (gencols == p)) {
      Rcpp::Rcerr << "Warning: genoype_data matrix is of dimension " << p << "x"
                  << p << ", assuming data is stored SNPxSample!";
    }
    const size_t N = SNPfirst ? genotype_data.cols() : genotype_data.rows();
    if (!SNPfirst) {

      Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> geno_map(
          &genotype_data(0, 0), N, p);
      calc_cov(geno_map, isGeno);
    } else {

      Eigen::Map<
          Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
          geno_map(&genotype_data(0, 0), N, p);
      calc_cov(geno_map, isGeno);
    }
    if (!is_sorted(mapd.data(), mapd.data() + mapd.size(), std::less<T>())) {
      Rcpp::stop("Recombination map must be non-decreasing\n");
    }
  }

  ldshrinker(Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> &S_,
             Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> &map_,
             const T m_, const T Ne_, const T cutoff_)
      : m(m_), theta(calc_theta(m)), cutoff(cutoff_), Ne(Ne_),
        mapd(map_.data(), map_.size()), S(S_.data(), S_.rows(), S_.cols()) {
    if (!is_sorted(mapd.data(), mapd.data() + mapd.size(), std::less<T>())) {
      Rcpp::stop("Recombination map must be non-decreasing\n");
    }
    const size_t p = mapd.size();
    if (S.cols() != p || S.rows() != p) {
      Rcpp::stop("Covariance Matrix have column(row) number equal to the size "
                 "of genetic map");
    }
  }

  T shrinkp(const T map_dist) {
    T rho = 4 * Ne * (map_dist) / 100;
    rho = -rho / (2 * m);
    rho = std::exp(rho);
    return (rho < cutoff ? 0 : rho);
  }
  // template<typename Derived,int Flags=Eigen::internal::traits<Derived>::Flags
  // & Eigen::RowMajorBit ? Eigen::RowMajor : Eigen::ColMajor>
  void cov_2_cor() {
    // typedef typename Derived::Scalar Scalar;
    Eigen::Array<T, Eigen::Dynamic, 1> rowvar = 1 / S.diagonal().array().sqrt();
    S.array().colwise() *= rowvar;
    S.array().rowwise() *= rowvar.transpose();
    S.diagonal().setOnes();
  }
  template <int RN>
  void
  calc_cov(Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, RN>> mat,
           const bool isGeno) {
    const T GenoMult = isGeno ? 0.5 : 1;
    // hmata = hmata.rowwise()-hmata.colwise().mean();
    mat = mat.rowwise() - mat.colwise().mean();
    S = (((mat.adjoint() * mat) / (mat.rows() - 1))) * GenoMult;
  }
  void Shrink() {
    const int numSNP = mapd.size();
    T ti = 0;
    T tshrinkage;
    for (int i = 0; i < numSNP; i++) {
      ti = mapd[i];
      for (int j = i + 1; j < numSNP; j++) {
        S(i, j) = S(i, j) * shrinkp(mapd[j] - ti);
      }
    }
    S.template triangularView<Eigen::StrictlyLower>() = S.transpose();
    S = ((1 - theta) * (1 - theta)) * S.array();
    for (int i = 0; i < numSNP; i++) {
      S(i, i) += 0.5 * theta * (1 - 0.5 * theta);
    }
  }
};
