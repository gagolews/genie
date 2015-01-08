#include <Rcpp.h>
using namespace Rcpp;
#include <stack>

class S : public std::stack<SEXP> {
   public: ~S() {
//      Rcout << "DESCRUCT!\n";
      while (!this->empty()) {
         SEXP obj = this->top();
         this->pop();
         R_ReleaseObject(obj);
      }
   }
};

// [[Rcpp::export]]
SEXP stack_create() {
   S* stack = new S();
   return XPtr< S >(stack, true);
}

// [[Rcpp::export]]
bool stack_empty(SEXP stack) {
   XPtr< S > _stack = Rcpp::as< XPtr< S > > (stack);
   return (*_stack).empty();
}

// [[Rcpp::export]]
void stack_push(SEXP stack, SEXP obj) {
   XPtr< S > _stack = Rcpp::as< XPtr< S > > (stack);
   R_PreserveObject(obj);
   (*_stack).push(obj);
}

// [[Rcpp::export]]
SEXP stack_pop(SEXP stack) {
   XPtr< S > _stack = Rcpp::as< XPtr< S > > (stack);
   if ((*_stack).empty())
      stop("empty stack");
   SEXP obj = (*_stack).top();
   (*_stack).pop();
   R_ReleaseObject(obj);
   return obj;
}


/***
q <- stack_create()
x1 <- 1:10
x2 <- c("gagas", "gdrher")
stack_push(q, x1)
stack_push(q, x2)
rm(x1)
rm(x2)
#rm(q)
gc(TRUE)
stack_pop(q)
stack_pop(q)
*/
