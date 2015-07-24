/* ************************************************************************* *
 *   This file is part of the `DataStructures` package.                      *
 *                                                                           *
 *   Copyright 2015 Maciej Bartoszuk, Anna Cena, Marek Gagolewski,           *
 *                                                                           *
 *   Parts of the code are taken from the 'CITAN' R package by M. Gagolewski *
 *                                                                           *
 *   'DataStructures' is free software: you can redistribute it and/or       *
 *   modify it under the terms of the GNU Lesser General Public License      *
 *   as published by the Free Software Foundation, either version 3          *
 *   of the License, or (at your option) any later version.                  *
 *                                                                           *
 *   'DataStructures' is distributed in the hope that it will be useful,     *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU Lesser General Public License for more details.                     *
 *                                                                           *
 *   You should have received a copy of the GNU Lesser General Public        *
 *   License along with 'DataStructures'.                                    *
 *   If not, see <http://www.gnu.org/licenses/>.                             *
 * ************************************************************************* */


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
   static const char* ClassName;
};

const char* Q::ClassName = "Queue";

void checkIsQueueClass(XPtr< Q >& _queue)
{
   if(strcmp(_queue.attr("class"), Q::ClassName))
      stop("not a Queue object");
}

//' @rdname queue
//' @details
//' \code{as.list.Queue} converts a given queue object to an R list.
//'
//' @return
//' \code{as.list.Queue} returns an R list object.
//' @param queue a queue object
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

//' @rdname queue
//' @details
//' \code{queue_create} creates a new queue.
//'
//' @return
//' \code{queue_create} returns a new, empty queue.
// [[Rcpp::export]]
SEXP queue_create() {
   Q* queue = new Q();
   XPtr< Q > retval =  XPtr< Q >(queue, true);
   retval.attr("class") = Q::ClassName;
   return retval;
}

//' @rdname queue
//' @details
//' \code{queue_empty} determines if a given queue is empty or not.
//'
//' @return
//' \code{queue_empty} returns a single logical value.
// [[Rcpp::export]]
bool queue_empty(SEXP queue) {
   XPtr< Q > _queue = Rcpp::as< XPtr< Q > > (queue);
   checkIsQueueClass(_queue);
   return (*_queue).empty();
}

//' @rdname queue
//' @details
//' \code{queue_push} pushes a given object to the back of the queue.
//'
//' @return
//' \code{queue_push} does not return anything interesting.
//'
//' @param obj an R object
// [[Rcpp::export]]
void queue_push(SEXP queue, SEXP obj) {
   XPtr< Q > _queue = Rcpp::as< XPtr< Q > > (queue);
   checkIsQueueClass(_queue);
   R_PreserveObject(obj);
   (*_queue).push_back(obj);
}

//' @rdname queue
//' @details
//' \code{queue_pop} pops an object from the top of the queue.
//'
//' @return
//' \code{queue_pop} returns an object at the top of te queue or
//' throws an error if the stack is empty.
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
