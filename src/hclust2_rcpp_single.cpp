/* ************************************************************************* *
 *   This file is part of the `grup` package.                                *
 *                                                                           *
 *   Copyright 2015 Maciej Bartoszuk, Anna Cena, Marek Gagolewski,           *
 *                                                                           *
 *   'grup' is free software: you can redistribute it and/or                 *
 *   modify it under the terms of the GNU Lesser General Public License      *
 *   as published by the Free Software Foundation, either version 3          *
 *   of the License, or (at your option) any later version.                  *
 *                                                                           *
 *   'grup' is distributed in the hope that it will be useful,               *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU Lesser General Public License for more details.                     *
 *                                                                           *
 *   You should have received a copy of the GNU Lesser General Public        *
 *   License along with 'grup'.                                              *
 *   If not, see <http://www.gnu.org/licenses/>.                             *
 * ************************************************************************* */


#include "hclust2_vptree_single.h"
#include "hclust2_naive_single.h"


// [[Rcpp::export(".hclust2_single")]]
RObject hclust2_single(RObject distance, RObject objects, RObject control=R_NilValue) {
   MESSAGE_2("[%010.3f] starting timer\n", clock()/(double)CLOCKS_PER_SEC);
   RObject result(R_NilValue);
   grup::Distance* dist = grup::Distance::createDistance(distance, objects, control);

   try { /* Rcpp::checkUserInterrupt(); may throw an exception */

      grup::HClustVpTreeSingle hclust(dist, control);
      // grup::HClustNaiveSingle hclust(dist, control);

      grup::HClustResult result2 = hclust.compute();
      result = Rcpp::as<RObject>(
         result2.toR(hclust.getStats(), hclust.getOptions(), dist->getStats())
      );

      // hclust.print();
   }
   catch(...) {
      // do nothing yet
   }

#if VERBOSE > 0
   dist->getStats().print();
#endif
   if (dist) delete dist;
   MESSAGE_2("[%010.3f] done\n", clock()/(double)CLOCKS_PER_SEC);
   if (Rf_isNull(result)) stop("stopping on error or explicit user interrupt");
   return result;
}


// // [[Rcpp::export]]
// RObject nntest(RObject distance, RObject objects, RObject control=R_NilValue) {
// #if VERBOSE > 5
//    Rprintf("[%010.3f] starting timer\n", clock()/(double)CLOCKS_PER_SEC);
// #endif
//    RObject result(R_NilValue);
//    grup::Distance* dist = grup::Distance::createDistance(distance, objects);
//
//    try {
//       /* Rcpp::checkUserInterrupt(); may throw an exception */
//       grup::HClustVpTreeSingle hclust(dist, control);
//       hclust.getNearestNeighbor(0);
//       hclust.getNearestNeighbor(25);
//       hclust.getNearestNeighbor(4653);
//       hclust.getNearestNeighbor(352);
//       result = Rcpp::as<RObject>(List::create(
//          _["merge"]  = R_NilValue,
//          _["height"] = R_NilValue,
//          _["order"]  = R_NilValue,
//          _["labels"] = R_NilValue,
//          _["call"]   = R_NilValue,
//          _["method"] = "single",
//          _["dist.method"] = R_NilValue,
//          _["stats"] = List::create(
//             _["vptree"] = hclust.getStats().toR(),
//             _["distance"] = dist->getStats().toR()
//          ),
//          _["control"] = List::create(
//             _["vptree"] = hclust.getOptions().toR()
//          )
//       ));
//    }
//    catch(...) {
//
//    }
//
//    if (dist) delete dist;
// #if VERBOSE > 5
//    Rprintf("[%010.3f] done\n", clock()/(double)CLOCKS_PER_SEC);
// #endif
//    if (Rf_isNull(result)) stop("stopping on error or explicit user interrupt");
//    return result;
// }
