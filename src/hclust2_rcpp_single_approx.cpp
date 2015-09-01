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

#include "hclust2_nnbased_single_approx.h"

// [[Rcpp::export(".hclust2_single_approx")]]
RObject hclust2_single_approx(RObject distance, RObject objects, RObject control=R_NilValue) {
   MESSAGE_2("[%010.3f] starting timer\n", clock()/(double)CLOCKS_PER_SEC);
   RObject result(R_NilValue);
   DataStructures::Distance* dist = DataStructures::Distance::createDistance(distance, objects);

   try { /* Rcpp::checkUserInterrupt(); may throw an exception */

      DataStructures::HClustNNbasedSingleApprox hclust(dist, control);
      // DataStructures::HClustNaiveSingle hclust(dist, control);

      DataStructures::HClustResult result2 = hclust.compute();
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
