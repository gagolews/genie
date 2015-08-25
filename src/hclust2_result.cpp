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

#include "hclust2_result.h"
using namespace DataStructures;
using namespace Rcpp;


HClustResult::HClustResult(size_t n, Distance* dist) :
      i(0),
      n(n),
      links(n-1, 2),         // this may be meaningless for some methods, do not return
      merge(n-1, 2),         // TO DO: merge determination during link()
      height(n-1),
      order(n, NA_REAL),     // TO DO: how to implement that?
      labels(n, NA_STRING),  // TO DO: set by a Distance object
      dist_method(NA_STRING) // TO DO: set by a Distance object
{

}


void HClustResult::link(size_t i1, size_t i2, double d12) {
   STOPIFNOT(i < n-1);
   links(i,0) = (double)i1;
   links(i,1) = (double)i2;
   height(i) = d12;
   ++i;
}


List HClustResult::toR(
      const HClustStats& hclustStats,
      const HClustOptions& hclustOptions,
      const DistanceStats& distStats)
{
   MESSAGE_2("[%010.3f] generating output matrix\n", clock()/(float)CLOCKS_PER_SEC);
   MergeMatrixGenerator mmg(n-1);
   merge = mmg.generateMergeMatrix(links);

   List result = List::create(
      _["merge"]  = merge,
      _["height"] = height,
      _["order"]  = order,
      _["labels"] = labels,
      _["call"]   = R_NilValue,
      _["method"] = R_NilValue,
      _["dist.method"] = dist_method,
      // _["links"]  = links,
      _["stats"] = List::create(
         _["method"] = hclustStats.toR(),
         _["distance"] = distStats.toR()
      ),
      _["control"] = List::create(
         _["method"] = hclustOptions.toR()
      )
   );

   result.attr("class") = "hclust";
   return result;
}
