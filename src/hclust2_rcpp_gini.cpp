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


#include "hclust2_vptree_gini.h"


// [[Rcpp::export(".hclust2_gini")]]
RObject hclust2_gini(RObject distance, RObject objects, RObject control=R_NilValue) {
   MESSAGE_2("[%010.3f] starting timer\n", clock()/(double)CLOCKS_PER_SEC);
   RObject result(R_NilValue);
   grup::Distance* dist = grup::Distance::createDistance(distance, objects, control);

   try { /* Rcpp::checkUserInterrupt(); may throw an exception */

      grup::HClustVpTreeGini hclust(dist, control);

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
