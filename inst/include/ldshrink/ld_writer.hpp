#pragma once
#include <RcppParallel.h>
#include <RcppEigen.h>

#include <atomic>
//enum { Anno, Cor };

// template<int Def>
// struct c2n{
//   static constexpr auto colname = "unknown";
// };

// template<>
// struct c2n<2>{
//   static constexpr auto colname = "anno";
// };

// template<>
// struct c2n<1>{
//   static constexpr auto colname = "r";
// };


template<typename F,int Def>
struct c2d{
  static constexpr std::array<F, Def> default_v = {-1};
};

template<>
struct c2d<double,1>{
  static constexpr std::array<double,1> default_v = {1};
};

template<>
struct c2d<float,1>{
  static constexpr std::array<double,1> default_v = {1};
};

template<>
struct c2d<float,2>{
  static constexpr std::array<float,2> default_v = {2,0};
};

template<>
struct c2d<double,2>{
  static constexpr std::array<double,2> default_v = {2,0};
};



template<typename T,typename F,int N,int AS> class MatTup{
public:
  std::array<T,N> m_idx;
  std::array<F,AS> m_value;
public:
  MatTup(const MatTup<T, F, N, AS> &other)
    : m_idx(other.m_idx), m_value(other.m_value) {};
  MatTup(const std::array<T,N> idx):m_idx(idx),m_value(c2d<F,AS>::default_v){
  };
  MatTup(const std::array<T, N> idx, const std::array<F, AS> value_)
      : m_idx(idx), m_value(value_) {}
  MatTup operator=(const MatTup<int, double, 1, 2> &other) {
    return (MatTup(other));
  }
  T row() const {
    static_assert((AS > 1),"Can't get row from a MatTup that stores only cols");
    return m_idx[0]; }
  T col() const { return m_idx[1]; }
  F value() const { return m_value[0] ;}
  template <int I> T rc() const { return m_idx[I] ;}
  template <int I> F val() const { return m_value[I];}
  //  const std::pair<T,T> idx() const {return std::make_pair(m_row, m_col);}
};


template<int AN>
class Skyline_data_store{
  using T=int;
  using F=double;
protected:
  std::vector<tbb::concurrent_vector<MatTup<T,F,1,AN>> > data_store;
  const size_t p_a;
  const size_t p_b;
  std::atomic<int> tot_elem;
  const bool isSymmetric;
public:
  Skyline_data_store(const size_t p_);
  Skyline_data_store(const size_t p_a, const size_t p_b);
  Skyline_data_store(const std::array<Rcpp::NumericMatrix,AN> mat_dat,const bool isSymmetric = true);
  void write_symm(const std::array<T,2> idx);
  void write_nsymm(const std::array<T, 2> idx, const std::array<F,2> res);
  Eigen::SparseMatrix<F> toSparseMatrix() const;
  Rcpp::StringVector fetch_colnames(const Rcpp::StringVector rownames) const;
  void finalize();
  void finalize_sort();
  SEXP toType(const std::string output_type,Rcpp::StringVector rownames = Rcpp::StringVector::create(),
	      Rcpp::StringVector colnames = Rcpp::StringVector::create())const;
  Rcpp::S4 todsCMatrix(Rcpp::StringVector rownames = Rcpp::StringVector::create(),
		   Rcpp::StringVector colnames = Rcpp::StringVector::create()) const;

  Rcpp::List toDataFrame(
				    Rcpp::StringVector rownames = Rcpp::StringVector::create(),
				    Rcpp::StringVector colnames = Rcpp::StringVector::create()) const ;
  Rcpp::NumericMatrix toMatrix(
				    Rcpp::StringVector rownames = Rcpp::StringVector::create(),
				    Rcpp::StringVector colnames = Rcpp::StringVector::create()) const ;


  Rcpp::List subDataFrame(const int i,const Rcpp::StringVector colnames)const ;

  template<int IN>
  class const_iterator{
    using data_t=typename std::vector<tbb::concurrent_vector<MatTup<T,F,1,AN>> >;

    using	row_it_t = typename data_t::const_iterator;
    const data_t &tdata;
    using	col_it_t=typename tbb::concurrent_vector<MatTup<int,double,1,AN>>::const_iterator;
    // using T=int;
    // using F=double;
  protected:
    const size_t tot_num;
    row_it_t row_it;
    col_it_t col_it;
    T cur_row;
    T cur_col;
    F cur_value;
    Eigen::Triplet<F> cur_el;
  public:
    const_iterator(const const_iterator &other)
      : tdata(other.tdata),tot_num(other.tot_num),row_it(other.row_it),col_it(other.col_it),cur_row(other.cur_row),cur_col(other.cur_col),cur_value(other.cur_value),cur_el(other.cur_el){
      // if (other.col_it != row_it->end()) {
      // 	cur_col = other.cur_col;
      // 	cur_value = other.cur_value;
      // 	cur_el = other.cur_el;
      // }else{
      // 	cur_col=-1;
      // }
    }
    const_iterator(const data_t &data_,const size_t tot_num_)
      : tdata(data_),
	tot_num(tot_num_),
	row_it(data_.begin()),
	col_it(row_it->begin()),
	cur_row(0){
      if (col_it != row_it->end()) {
	cur_col=col_it->template rc<0>();
	cur_value=col_it->template val<IN>();
      } else {
        cur_col = -1;
        cur_value = -1;
      }
      cur_el = Eigen::Triplet<F>(cur_row, cur_col, cur_value);
    }
    const_iterator operator++(){
      col_it++;
      if (col_it == row_it->end()) {
	cur_row++;
	row_it++;
        if (row_it == tdata.end()) {
          return *this;
        }
	col_it=row_it->begin();
      }
      if (col_it == row_it->end()) {
        return ++(*this);
      }
      cur_col = col_it->template rc<0>();
      cur_value = col_it->template val<IN>();
      cur_el = Eigen::Triplet<F>(cur_row, cur_col, cur_value);
      return *this;
    }
    const Eigen::Triplet<F>& operator *(){
      // if (cur_num>=tot_num) {
      //   Rcpp::Rcerr << "In element " << cur_num << std::endl;
      //   Rcpp::stop("can't derefernce");
      // }
      return (cur_el);
    }
    const Eigen::Triplet<F>* operator ->(){
      // if (cur_num >= tot_num) {
      //   Rcpp::Rcerr << "In element " << cur_num << std::endl;
      //   Rcpp::stop("can't derefernce");
      // }
      return(&cur_el);
    }
    bool operator==(const const_iterator &rhs)const {
      return col_it == rhs.col_it;
    }
    bool operator!=(const const_iterator &rhs) const {
      return col_it != rhs.col_it;
    }
    friend class Skyline_data_store;
  };

  template<int IN>
  const_iterator<IN> begin() const{
    return(const_iterator<IN>(data_store,tot_elem));
  }
  template <int IN> const_iterator<IN> end() const {
    const size_t cur_tot_elem = tot_elem;
    const_iterator<IN> ret(data_store, cur_tot_elem);
    std::advance(ret.row_it, data_store.size() - 1);
    ret.col_it = ret.row_it->end();
   // ret.cur_num = cur_tot_elem;
    return ret;
  }
};



  template <typename T, typename F, int AN> class Pvec_data_store {
    tbb::concurrent_vector<MatTup<T, F, 2, AN>> data_store;
    const size_t p_a;
    const size_t p_b;

  public:
    Pvec_data_store(const size_t p) : p_a(p), p_b(p) {
      data_store.reserve(((p * p - p) / 2 + p) / 5);
    }
    Pvec_data_store(const size_t p_aa, const size_t p_bb)
        : p_a(p_aa), p_b(p_bb) {
      const size_t tot_p = p_a * p_b;
      data_store.reserve(tot_p / 10);
    }
    void write_symm(const std::array<T, 2> idx) {
      data_store.push_back(MatTup<T, F, 2, AN>(idx));
    }
    void write_nsymm(const std::array<T, 2> idx, const std::array<F, 2> res) {
      std::array<F, AN> tres;
      std::copy_n(res.begin(), AN, tres.begin());
      data_store.push_back(MatTup<T, F, 2, AN>(idx, tres));
    }
    void finalize() {}

    Eigen::SparseMatrix<double> toSparseMatrix() const {
      const size_t num_elem = data_store.size();
      Eigen::SparseMatrix<double> object(p_a, p_b);
      // Rcpp::Rcerr << "About to create sparse matrix of size: ";
      // Rcpp::Rcerr << num_elem << std::endl;
      object.setFromTriplets(data_store.begin(), data_store.end());
      //    Rcpp::Rcerr << "Sparse Matrix created" << std::endl;
      return (object);
    }

    Rcpp::NumericMatrix
    toMatrix(Rcpp::StringVector rownames = Rcpp::StringVector::create(),
             Rcpp::StringVector colnames = Rcpp::StringVector::create()) const {
      const size_t num_elem = data_store.size();
      Rcpp::NumericMatrix retmat(p_a, p_b);
      for (auto &tr : data_store) {
        retmat(tr.row(), tr.col()) = tr.value();
      }
      retmat.attr("dimnames") = Rcpp::List::create(rownames, colnames);
      return (retmat);
    }

    SEXP toDataFrame(
        Rcpp::StringVector rownames = Rcpp::StringVector::create(),
        Rcpp::StringVector colnames = Rcpp::StringVector::create()) const {
      const size_t num_elem = data_store.size();
      using namespace Rcpp;
      const bool use_rownames = colnames.size() == 0;
      if (use_rownames) {
        Rcpp::StringVector colname(num_elem);
        Rcpp::StringVector rowname(num_elem);
        std::array<Rcpp::NumericVector, AN> value_vec{};
        for (int k = 0; k < AN; k++) {
          value_vec[k] = Rcpp::NumericVector(num_elem);
        }
        for (int j = 0; j < num_elem; j++) {
          const auto &tref = data_store[j];
          const auto &trv = tref.m_value;
          for (int k = 0; k < AN; k++) {
            value_vec[k][j] = trv[k];
          }
          colname[j] = colnames[tref.template rc<1>()];
          rowname[j] = colnames[tref.template rc<0>()];
        }

        if (AN == 2) {
          return (List::create(_["rowsnp"] = rowname, _["colsnp "] = colname,
                               _["anno"] = value_vec[1],
                               _["r"] = value_vec[0]));
        } else {
          return (List::create(_["rowsnp"] = rowname, _["colsnp"] = colname,
                               _["r"] = value_vec[0]));
        }
      }
      Rcpp::IntegerVector colname(num_elem);
      Rcpp::IntegerVector rowname(num_elem);
      std::array<Rcpp::NumericVector, AN> value_vec{};
	for(int k=0; k<AN; k++){
	  value_vec[k]=Rcpp::NumericVector(num_elem);
	}
	for  ( int  j  = 0; j  < num_elem; j++) {
	  const auto &tref = data_store[j];
	  const auto &trv=tref.m_value;
	  for (int k = 0; k < AN; k++) {
	    value_vec[k][j] = trv[k];
	  }
	  colname[j] = tref.template rc<1>();
	  rowname[j] = tref.template rc<0>();
	}
    if(AN==2){
      return(List::create(_["rowsnp"]=rowname,_["colsnp"]=colname,_["anno"]=value_vec[1],_["r"]=value_vec[0]));
    } else {
      return (List::create(_["rowsnp"]=rowname,_["colsnp"] = colname, _["r"] = value_vec[0]));
    }
  }
};



template<int AN>
inline Skyline_data_store<AN>::Skyline_data_store(const std::array<Rcpp::NumericMatrix,AN> mat_dat,const bool is_symmetric): data_store(mat_dat[0].nrow()), p_a(mat_dat[0].nrow()), p_b(mat_dat[0].ncol()), tot_elem(0),isSymmetric(is_symmetric) {
  std::array<double,2>	telem;
  if(!isSymmetric){
    for (int i = 0; i < p_a; i++) {
      data_store[i].reserve(p_b);
      for(int j=0; j<p_b; j++){
	for(int k=0; k<AN; k++){
	  telem[k]=mat_dat[k](i,j);
	}
        write_nsymm({i, j}, telem);
      }
    }
  } else {
    for (int i = 0; i < p_a; i++) {
      data_store[i].reserve(p_b);
      for(int j=i; j<p_b; j++){
	for(int k=0; k<AN; k++){
	  telem[k]=mat_dat[k](i,j);
	}
        write_nsymm({i, j}, telem);
      }
    }
  }
}

template <int AN>
inline Skyline_data_store<AN>::Skyline_data_store(const size_t p_)
  : data_store(p_), p_a(p_), p_b(p_), tot_elem(0),isSymmetric(true) {
  for (size_t i = 0; i < p_a; i++) {
    data_store[i].reserve(p_a  / 5);
  }
}

template<int AN>
inline Skyline_data_store<AN>::Skyline_data_store(const size_t p_aa, const size_t p_bb):data_store(p_a),p_a(p_aa),p_b(p_bb),tot_elem(0),isSymmetric(false){
  for (size_t i = 0; i < p_a; i++) {
    data_store[i].reserve(p_b / 10);
  }
  const size_t tot_p = p_a * p_b;
  data_store.reserve(tot_p / 10);
}

template<int AN>
inline void Skyline_data_store<AN>::write_symm(const std::array<T,2> idx){
  data_store[idx[0]].push_back(MatTup<T, F, 1, AN>({idx[1]}));
  tot_elem.fetch_add(1, std::memory_order_relaxed);
}

template<int AN>
inline void Skyline_data_store<AN>::write_nsymm(const std::array<T, 2> idx, const std::array<F,2> res) {
  std::array<F,AN> tres;
  std::copy_n(res.begin(),AN,tres.begin());
    data_store[idx[0]].push_back(MatTup<T, F, 1, AN>({idx[1]},tres));
  tot_elem.fetch_add(1, std::memory_order_relaxed);
}

template <int AN>
inline void Skyline_data_store<AN>::finalize(){
  this->finalize_sort();
}

template<int AN>
inline void Skyline_data_store<AN>::finalize_sort(){
  using elt=MatTup<T,F,1,AN>;
  tbb::parallel_for(tbb::blocked_range<size_t>(0, p_a),
                    [&](const tbb::blocked_range<size_t> &r) {
                      for (size_t i = r.begin(); i != r.end(); i++) {
			std::sort(
                            data_store[i].begin(), data_store[i].end(),
                            [](const elt &a, const elt &b) {
                              return (a.m_idx < b.m_idx);
                            });
		      }
		    });
}

template <int AN>
inline Rcpp::NumericMatrix
Skyline_data_store<AN>::toMatrix(Rcpp::StringVector rownames,
                                      Rcpp::StringVector colnames) const {
  Rcpp::NumericMatrix ret(p_a,p_b);
  if(isSymmetric){
    for (int i = 0; i < p_a; i++) {
      const auto & dsr = data_store[i];
      const size_t tp = dsr.size();
      for (int j = 0; j < tp; j++) {
        int tcol = dsr[j].template rc<0>();
        ret(i, tcol) = dsr[j].template val<0>();
        ret(tcol, i) = ret(i, tcol);
      }
    }
  } else {
    for (int i = 0; i < p_a; i++) {
      const auto & dsr = data_store[i];
      const size_t tp = dsr.size();
      for (int j = 0; j < tp; j++) {
        int tcol = dsr[j].template rc<0>();
        ret(i, tcol) = dsr[j].template val<0>();
      }
    }
  }
  if (rownames.size() ==p_a && colnames.size()==p_b) {
    ret.attr("dimnames") = Rcpp::List::create(rownames, colnames);
  }
  return (ret);
}

template<int AN>
inline SEXP Skyline_data_store<AN>::toType(const std::string output_type,Rcpp::StringVector rownames,Rcpp::StringVector colnames)const{
  if (output_type == "dsCMatrix") {
    return (this->todsCMatrix(rownames, colnames));
  }
  if (output_type == "matrix") {
    return (this->toMatrix(rownames, colnames));
  }
  if (output_type == "data.frame" || output_type == "data_frame") {
    return (this->toDataFrame(rownames, colnames));
  }
  Rcpp::Rcerr << "Unknown type " << output_type << std::endl;
  Rcpp::stop("output type must be one of dsCMatrix,matrix,data.frame");
}

template<int AN>
inline Rcpp::S4 Skyline_data_store<AN>::todsCMatrix(Rcpp::StringVector rownames,
					 Rcpp::StringVector colnames) const {
  if(!isSymmetric){
    Rcpp::stop("can't store as dsCMatrix if data is not symmetric");
  }
  auto object = this->toSparseMatrix();
  using namespace Rcpp;
  const int nnz = object.nonZeros();
  //  Rcpp::Rcerr << "Number of entries: " << nnz << std::endl;
  S4 ans("dsCMatrix");
  ans.slot("Dim") = Dimension(object.rows(), object.cols());
  ans.slot("i") =
      IntegerVector(object.innerIndexPtr(), object.innerIndexPtr() + nnz);
  //  ans.slot("uplo")="L";
  ans.slot("p") = IntegerVector(object.outerIndexPtr(), object.outerIndexPtr() + object.outerSize() + 1);
  ans.slot("x") = NumericVector(object.valuePtr(), object.valuePtr() + nnz);
    if (rownames.size() ==p_a && colnames.size()==p_b) {
    ans.slot("Dimnames") = List::create(rownames, colnames);
  }
  return (ans);
}

template<int AN>
Eigen::SparseMatrix<double> Skyline_data_store<AN>::toSparseMatrix() const{
  Eigen::ArrayXi OuterStarts(p_b+1);
  const size_t cur_tot_elem=tot_elem;
  // Eigen::ArrayXd mat_values(cur_tot_elem);
  // Eigen::ArrayXi InnerIndices(cur_tot_elem);
  // Eigen::ArrayXi InnerNNZ(p_a);

  //  Rcpp::NumericVector	retvec(tot_elem);
  //  auto retvec_b=retvec.begin();
  Eigen::SparseMatrix<double> tr(p_a,p_b);
  tr.setFromTriplets(this->begin<0>(),this->end<0>());
  return(tr);
}

template<int AN>
inline Rcpp::List Skyline_data_store<AN>::toDataFrame(Rcpp::StringVector rownames,
							       Rcpp::StringVector colnames) const {
  using namespace Rcpp;
  CharacterVector classes = CharacterVector::create("tbl_df", "data.frame");
  const bool use_rownames=rownames.size()>0;
  List DataList = no_init(p_a);
  for (int i = 0; i < p_a; i++) {
    DataList[i] = subDataFrame(i, colnames);
  }
  Rcpp::IntegerVector row_id(p_a);
  std::iota(row_id.begin(), row_id.end(), 1);
  if (!use_rownames) {
    List retl = List::create(_["rowsnp"] = row_id, _["data"] = DataList);
    retl.attr("class") = classes;
    retl.attr("row.names")=row_id;
    return (retl);
  }
  List retl = List::create(_["rowsnp"] = rownames, _["data"] = DataList);
  retl.attr("row.names")=row_id;
  retl.attr("class") = classes;

  return (retl);
}

template<int AN>
inline Rcpp::List Skyline_data_store<AN>::subDataFrame(const int i,const Rcpp::StringVector colnames)const {
  using namespace Rcpp;
  const bool use_rownames=colnames.size()>0;
  const size_t tp=data_store[i].size();
  CharacterVector classes = CharacterVector::create("tbl_df", "data.frame");
  Rcpp::IntegerVector row_id(tp);
  std::iota(row_id.begin(), row_id.end(), 1);
  if(use_rownames){
    Rcpp::StringVector colname(tp);
    std::array<Rcpp::NumericVector, AN> value_vec{};
    for(int k=0; k<AN; k++){
      value_vec[k]=Rcpp::NumericVector(tp);
    }
    for (int j = 0; j < tp; j++) {
      const auto &tref = data_store[i][j];
      const auto &trv=tref.m_value;
      for (int k = 0; k < AN; k++) {
        value_vec[k][j] = trv[k];
      }
      colname[j] = colnames[tref.template rc<0>()];
    }

    if(AN==2){
      List retl=List::create(_["colsnp"]=colname,_["anno"]=value_vec[1],_["r"]=value_vec[0]);
      retl.attr("class")=classes;
      retl.attr("row.names")=row_id;
      return(retl);
    } else {

      List retl = List::create(_["colsnp"] = colname, _["r"] = value_vec[0]);
      retl.attr("class") = classes;
      retl.attr("row.names")=row_id;
      return (retl);
    }
  }
  Rcpp::IntegerVector colname(tp);
  std::array<Rcpp::NumericVector, AN> value_vec{};
  for(int k=0; k<AN; k++){
    value_vec[k]=Rcpp::NumericVector(tp);
  }
  for (int j = 0; j < tp; j++) {
    const auto &tref = data_store[i][j];
    const auto &trv=tref.m_value;
    for (int k = 0; k < AN; k++) {
      value_vec[k][j] = trv[k];
    }
    colname[j] = tref.template rc<0>()+1;
  }
 if(AN==2){
      List retl=List::create(_["colsnp"]=colname,_["anno"]=value_vec[1],_["r"]=value_vec[0]);
      retl.attr("class")=classes;
      retl.attr("row.names")=row_id;
      return(retl);
    } else {
      List retl = List::create(_["colsnp"] = colname, _["r"] = value_vec[0] );
      retl.attr("class") = classes;
      retl.attr("row.names")=row_id;
      return (retl);
    }
}

// template<typename S>

// class LDshrinkWr  i ter{
// protected:
//   using T=int;
//   tbb::concurrent_vector<MatTup<T> > data_store;
//   tbb::concurrent_vector<MatTup<T> > anno_data_store;
//   using	name_vec=Rcpp::StringVector;
//   const name_vec snpnames_row;
//   const name_vec snpnames_col;
//   const size_t p_row;
//   const size_t p_col;
//   const bool write_annotations;
// public:
//   LDshrinkWriter(const std::vector<T> &inp,const name_vec names,const bool write_annotations=false);
//   LDshrinkWriter(const std::pair<std::vector<T>,std::vector<T>> &inpconst, const std::pair<name_vec,name_vec> names,const bool write_annotations=false);
//   void write_symm(const std::pair<T,T> idx);
//   void write_nsymm(const std::pair<T,T> idx,const double res,const double annom);
//   Eigen::SparseMatrix<double> sparseMatrix() const;
//   void finalize();
//   SEXP dsCMatrix() const;
//   Rcpp::DataFrame toDataFrame(const bool use_rownames=true)const;
//   Rcpp::DataFrame toAnnotationDataFrame(const bool use_rownames)const;
//   SEXP toType(const std::string	output_type)const;
// };
