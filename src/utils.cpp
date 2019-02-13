#include <ldshrink.hpp>
#include <Rcpp.h>
#include <Rinternals.h>

#if __cplusplus > 201402L

namespace Rcpp{
  namespace traits{



    SEXP wrap(const std::vector<std::string_view> & obj){
      const size_t p=obj.size();
      SEXP out = PROTECT(Rf_allocVector(STRSXP, p));
      for(size_t i=0; i<p;i++){
	SET_STRING_ELT(out, i, Rf_mkChar(obj[i].data()));
      }
      UNPROTECT(1);

      return(out);
    };
        template<> class Exporter< std::vector<std::string_view> > {
            typedef typename std::vector<std::string_view> OUT ;

            // Convert the type to a valid rtype.
           // const static int RTYPE = Rcpp::traits::r_sexptype_traits< T >::rtype ;
            SEXP vec;

            public:
            Exporter(SEXP x):vec(x) {
                if (TYPEOF(x) != STRSXP)
                    throw std::invalid_argument("Wrong R type for mapped 1D array");
            }
            OUT get() {
                const size_t p=XLENGTH(vec);
                // Need to figure out a way to perhaps do a pointer pass?
                OUT x(p);

                for(size_t i=0; i<p; i++){
                  x[i]=std::string_view(CHAR(STRING_ELT(vec, i)));
                }
                return x;
            }
        };
    }
}
#endif

inline char flip(const char t) noexcept{
  switch(t){
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  default:
    return '0';
  }
};



int allele_check(const char* query,const char* alt){
  const char q_ref=query[0];
  const char q_alt=query[2];
  const char t_ref=alt[0];
  const char t_alt=alt[2];
  // Alt and ref being equal is a no-go
  if((q_ref==q_alt) or (t_ref==t_alt)){
    return 0;
  }
  
  //2 scenarios if first letters match
  if(q_ref==t_ref){
    //matches in both places
    if(q_alt==t_alt){
      return 1;
    }
    //matches in only one place
    return 0;
  }
  //this also falls under only matches in one place
  if(q_alt==t_alt){
    return 0;
  }
  //"classic" flip
  if(q_ref==t_alt){
    if(q_alt==t_ref){
      return -1;
    }
    return 0;
  }
  //strand flipped ref matches ref
  if((flip(q_ref)==t_ref)){
    //strand flipped alt mathces too, this means flip
    if(flip(q_alt) == t_alt){
      return -1;
    }
    
  }
  //Last chance for a match
  if(flip(q_ref)==t_alt){
    if(flip(t_ref)==q_alt){
      return 1;
    }
  }
  //Neither match
  return 0;

}




//' Determine whether 2 alleles are compatible
// '@param query_ref_alt is a ref/alt pair
// '@param target_ref_alt is another ref/alt pair
// '@return returns a vector with 1 if the query matches the target, -1 if a flip is required, or 0 if they are incompatible;
//' @export
//[[Rcpp::export]]
Rcpp::StringVector strand_flip(Rcpp::StringVector ref_alt,bool reverse=false){
  Rcpp::StringVector retvec(ref_alt.size());
  
  
  if(!reverse){
    std::transform(ref_alt.begin(),ref_alt.end(),retvec.begin(),[](const char* query)  {
      if(strlen( query)!=3 ){
        Rcpp::stop("alleles must all be length 3");
      } 
      const char ret[]={flip(query[0]),',',flip(query[2])};
      return(Rcpp::String(ret));
    });
  }else{
    std::transform(ref_alt.begin(),ref_alt.end(),retvec.begin(),[](const char* query)  {
      if(strlen( query)!=3 ){
        Rcpp::stop("alleles must all be length 3");
      } 
      const char ret[]={flip(query[2]),',',flip(query[0])};
      return(Rcpp::String(ret));
    });
  }
  
  
  return(retvec);                                        
}
  

//' Determine whether 2 
// '@param query_ref_alt is a ref/alt pair
// '@param target_ref_alt is another ref/alt pair
// '@return returns a vector with 1 if the query matches the target, -1 if a flip is required, or 0 if they are incompatible;
//' @export
//[[Rcpp::export]]
Rcpp::IntegerVector flip_alleles(Rcpp::StringVector query_ref_alt,Rcpp::StringVector target_ref_alt ){
  
  if(query_ref_alt.size()!=target_ref_alt.size()){
    Rcpp::stop("query ref_alt and target ref_alt must be the same size!");
  }
  

  
  Rcpp::IntegerVector retvec(target_ref_alt.size());
  std::transform(query_ref_alt.begin(),query_ref_alt.end(),target_ref_alt.begin(),retvec.begin(),[](const char* query,const char* alt) {
    if(strlen( query)!=3 || strlen(alt)!=3 ){
      Rcpp::stop("alleles must all be length 3");
    } 
    return(allele_check(query,alt));
  });
    
  return(retvec);                                        
}
                                            
  



