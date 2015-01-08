#include <Rcpp.h>
using namespace Rcpp;
#include <queue>

class Q : public std::queue<SEXP> {
   public: ~Q() {
//      Rcout << "DESCRUCT!\n";
      while (!this->empty()) {
         SEXP obj = this->front();
         this->pop();
         R_ReleaseObject(obj);
      }
   }
};

// [[Rcpp::export]]
SEXP queue_create() {
   Q* queue = new Q();
   return XPtr< Q >(queue, true);
}

// [[Rcpp::export]]
bool queue_empty(SEXP queue) {
   XPtr< Q > _queue = Rcpp::as< XPtr< Q > > (queue);
   return (*_queue).empty();
}

// [[Rcpp::export]]
void queue_push(SEXP queue, SEXP obj) {
   XPtr< Q > _queue = Rcpp::as< XPtr< Q > > (queue);
   R_PreserveObject(obj);
   (*_queue).push(obj);
}

// [[Rcpp::export]]
SEXP queue_pop(SEXP queue) {
   XPtr< Q > _queue = Rcpp::as< XPtr< Q > > (queue);
   if ((*_queue).empty())
      stop("empty queue");
   SEXP obj = (*_queue).front();
   (*_queue).pop();
   R_ReleaseObject(obj);
   return obj;
}


/***
q <- queue_create()
x1 <- 1:10
x2 <- c("gagas", "gdrher")
queue_push(q, x1)
queue_push(q, x2)
rm(x1)
rm(x2)
#rm(q)
gc(TRUE)
queue_pop(q)
queue_pop(q)
*/
