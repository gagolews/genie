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


#include "hclust2_vptree_single.h"



// [[Rcpp::export(".hclust2_single")]]
RObject hclust2_single(RObject distance, RObject objects, RObject control=R_NilValue) {
   MESSAGE_2("[%010.3f] starting timer\n", clock()/(double)CLOCKS_PER_SEC);
   RObject result(R_NilValue);
   DataStructures::Distance* dist = DataStructures::Distance::createDistance(distance, objects);

   try {
      /* Rcpp::checkUserInterrupt(); may throw an exception */
      DataStructures::HClustVpTreeSingle hclust(dist, control);
      RObject merge = hclust.compute();
      result = Rcpp::as<RObject>(List::create(
         _["merge"]  = merge,
         _["height"] = R_NilValue,
         _["order"]  = R_NilValue,
         _["labels"] = R_NilValue,
         _["call"]   = R_NilValue,
         _["method"] = "single",
         _["dist.method"] = R_NilValue,
         _["stats"] = List::create(
            _["vptree"] = hclust.getStats().toR(),
            _["distance"] = dist->getStats().toR()
         ),
         _["control"] = List::create(
            _["vptree"] = hclust.getOptions().toR()
         )
      ));
      result.attr("class") = "hclust";
      //hclust.print();
   }
   catch(...) {

   }

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
//    DataStructures::Distance* dist = DataStructures::Distance::createDistance(distance, objects);
//
//    try {
//       /* Rcpp::checkUserInterrupt(); may throw an exception */
//       DataStructures::HClustVpTreeSingle hclust(dist, control);
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
