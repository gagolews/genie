#include <Rcpp.h>
using namespace Rcpp;
#include <deque>

class Q : public std::deque<SEXP> {
   public: ~Q() {
//      Rcout << "DESCRUCT!\n";
      while (!this->empty()) {
         SEXP obj = this->front();
         this->pop_front();
         R_ReleaseObject(obj);
      }
   }
};

void checkIsQueueClass(XPtr< Q >& _queue)
{
   if(strcmp(_queue.attr("class"), "Queue"))
      stop("not a Queue object");
}

// [[Rcpp::export("as.list.Queue")]]
List queue_as_list(SEXP queue) {
   XPtr< Q > _queue = Rcpp::as< XPtr< Q > > (queue);
   checkIsQueueClass(_queue);
   int n = (*_queue).size();
   List retval(n);
   std::deque<SEXP>::iterator it = (*_queue).begin();
   for (int i=0; i < n; ++i)
      retval[i] = *(it++);
   return retval;
}


// [[Rcpp::export]]
SEXP queue_create() {
   Q* queue = new Q();
   XPtr< Q > retval =  XPtr< Q >(queue, true);
   retval.attr("class") = "Queue";
   return retval;
}

// [[Rcpp::export]]
bool queue_empty(SEXP queue) {
   XPtr< Q > _queue = Rcpp::as< XPtr< Q > > (queue);
   checkIsQueueClass(_queue);
   return (*_queue).empty();
}

// [[Rcpp::export]]
void queue_push(SEXP queue, SEXP obj) {
   XPtr< Q > _queue = Rcpp::as< XPtr< Q > > (queue);
   checkIsQueueClass(_queue);
   R_PreserveObject(obj);
   (*_queue).push_back(obj);
}

// [[Rcpp::export]]
SEXP queue_pop(SEXP queue) {
   XPtr< Q > _queue = Rcpp::as< XPtr< Q > > (queue);
   checkIsQueueClass(_queue);
   if ((*_queue).empty())
      stop("empty queue");
   SEXP obj = (*_queue).front();
   (*_queue).pop_front();
   R_ReleaseObject(obj);
   return obj;
}
