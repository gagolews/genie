#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

class V : public std::vector<RObject> {
   public:
   V(int n):std::vector<RObject>(n) {}

   ~V() {
      while (!this->empty()) {
         RObject obj = this->back();
         this->pop_back();
         R_ReleaseObject(obj);
      }
   }
   static const char* ClassName;
};

const char* V::ClassName = "Vector";

void checkIsVectorClass(XPtr< V >& _stack)
{
   if(strcmp(_stack.attr("class"), V::ClassName))
      stop("not a Vector object");
}

// [[Rcpp::export("as.list.Vector")]]
List vector_as_list(SEXP vec) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   int n = (*_vec).size();
   List retval(n);
   for (int i=0; i < n; ++i)
      retval[i] = (*_vec)[i];
   return retval;
}

// [[Rcpp::export]]
SEXP vector_create(int n=0) {
   V* vec = new V(n);
   XPtr< V > retval =  XPtr< V >(vec, true);
   retval.attr("class") = V::ClassName;
   return retval;
}

// [[Rcpp::export]]
bool vector_empty(SEXP vec) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   return (*_vec).empty();
}

// [[Rcpp::export]]
int vector_size(SEXP vec) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   return (*_vec).size();
}

// [[Rcpp::export]]
RObject vector_front(SEXP vec) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   if ((*_vec).empty())
      stop("empty vector");
   return (*_vec).front();
}

// [[Rcpp::export]]
RObject vector_back(SEXP vec) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   if ((*_vec).empty())
      stop("empty vector");
   return (*_vec).back();
}

// [[Rcpp::export]]
void vector_push_back(SEXP vec, RObject obj) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   (*_vec).push_back(obj);
}

// [[Rcpp::export]]
void vector_pop_back(SEXP vec) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   if ((*_vec).empty())
      stop("empty vector");
   (*_vec).pop_back();
}

// [[Rcpp::export]]
RObject& vector_at(SEXP vec, int i) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   return (*_vec).at(i-1);
}

// [[Rcpp::export]]
void vector_set_at(SEXP vec, int i, RObject obj) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   //Rcout << vec <<" "<< i<<" " << obj << std::endl;
   (*_vec).at(i-1)=obj;
   //Rcout << vec <<" "<< i<<" " << obj << std::endl;
}

