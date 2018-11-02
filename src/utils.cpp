#include <ldshrink.hpp>
#include <Rcpp.h>
#include <Rinternals.h>


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
