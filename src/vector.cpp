/* ************************************************************************* *
 *   This file is part of the `DataStructures` package.                      *
 *                                                                           *
 *   Copyright 2015 Maciej Bartoszuk, Anna Cena, Marek Gagolewski,           *
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

//' @rdname vector
//' @details
//' \code{as.list.Vector} converts a given vector object to an R list.
//'
//' @return
//' \code{as.list.Vector} returns an R list object.
//' @param vec a vector object
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

//' @rdname vector
//' @details
//' \code{vector_create} creates a new vector.
//'
//' @return
//' \code{vector_create} returns a new, empty vector.
// [[Rcpp::export]]
SEXP vector_create(int n=0) {
   V* vec = new V(n);
   XPtr< V > retval =  XPtr< V >(vec, true);
   retval.attr("class") = V::ClassName;
   return retval;
}

//' @rdname vector
//' @details
//' \code{vector_empty} determines if a given stack is empty or not.
//'
//' @return
//' \code{vector_empty} returns a single logical value.
// [[Rcpp::export]]
bool vector_empty(SEXP vec) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   return (*_vec).empty();
}

//' @rdname vector
//' @details
//' \code{vector_size} returns a size of vector.
//'
//' @return
//' \code{vector_size} returns a single integer value.
// [[Rcpp::export]]
int vector_size(SEXP vec) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   return (*_vec).size();
}

//' @rdname vector
//' @details
//' \code{vector_front} returns the first object of a given vector.
//'
//' @return
//' \code{vector_front} returns a single RObject value or
//' throws an error if a vector is empty.
// [[Rcpp::export]]
RObject& vector_front(SEXP vec) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   if ((*_vec).empty())
      stop("empty vector");
   return (*_vec).front();
}

//' @rdname vector
//' @details
//' \code{vector_back} returns the last object of a given vector.
//'
//' @return
//' \code{vector_back} returns a single RObject value or
//' throws an error if a vector is empty.
// [[Rcpp::export]]
RObject& vector_back(SEXP vec) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   if ((*_vec).empty())
      stop("empty vector");
   return (*_vec).back();
}

//' @rdname stack
//' @details
//' \code{vector_push_back} pushes a given object at the end of the vector.
//'
//' @return
//' \code{vector_push_back} does not return anything interesting.
//'
//' @param obj an R object
// [[Rcpp::export]]
void vector_push_back(SEXP vec, RObject obj) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   R_PreserveObject(obj);
   (*_vec).push_back(obj);
}

//' @rdname stack
//' @details
//' \code{vector_pop_back} removes the last element in the vector, effectively
//' reducing the container size by one.
//'
//' @return
//' \code{vector_pop_back} does not return anything interesting. Throws an
//' error if a vector is empty.
// [[Rcpp::export]]
void vector_pop_back(SEXP vec) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   if ((*_vec).empty())
      stop("empty vector");
   RObject  obj = (*_vec).back();
   (*_vec).pop_back();
   R_ReleaseObject(obj);
}

//' @rdname vector
//' @details
//' \code{vector_at} returns an object at position i in the vector.
//'
//' @return
//' \code{vector_at} returns a single RObject value.
//'
//' @param i position i in the vector
// [[Rcpp::export]]
RObject& vector_at(SEXP vec, int i) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   return (*_vec).at(i-1);
}

//' @rdname vector
//' @details
//' \code{vector_set_at} sets an object obj at position i in the vector.
//'
//' @return
//' \code{vector_set_at} does not return anything interesting.
//'
//' @param i position i in the vector
//' @param obj an R object
// [[Rcpp::export]]
void vector_set_at(SEXP vec, int i, RObject obj) {
   XPtr< V > _vec = Rcpp::as< XPtr< V > > (vec);
   checkIsVectorClass(_vec);
   R_PreserveObject(obj);
   R_ReleaseObject((*_vec).at(i-1));
   (*_vec).at(i-1)=obj;
}

