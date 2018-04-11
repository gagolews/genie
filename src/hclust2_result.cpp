/* ************************************************************************* *
 *   This file is part of the `genie` package for R.                         *
 *                                                                           *
 *   Copyright 2015-2018 Marek Gagolewski, Maciej Bartoszuk, Anna Cena       *
 *                                                                           *
 *   'genie' is free software: you can redistribute it and/or                *
 *   modify it under the terms of the GNU General Public License             *
 *   as published by the Free Software Foundation, either version 3          *
 *   of the License, or (at your option) any later version.                  *
 *                                                                           *
 *   'genie' is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with 'genie'. If not, see <http://www.gnu.org/licenses/>.         *
 * ************************************************************************* */

#include "hclust2_result.h"
using namespace grup;
using namespace Rcpp;


HClustResult::HClustResult(size_t n, Distance* dist, bool lite) :
      curiter(0),
      n(n),
      links(n-1, 2),         // this may be meaningless for some methods, do not return
      merge(n-1, 2),
      height(n-1),
      order(n, NA_REAL),
      labels(dist->getLabels()),
      dist_method(dist->getDistMethod()),
      lite(lite)
{
   // no-op
}


void HClustResult::link(size_t i1, size_t i2, double d12)
{
   STOPIFNOT(curiter < n-1);
   links(curiter, 0) = (double)i1;
   links(curiter, 1) = (double)i2;
   height(curiter) = d12;
   ++curiter;

   if (curiter == n-1 && !lite) {
      generateMergeMatrix();
      generateOrderVector();
   }
}


void HClustResult::generateOrderVector()
{
   std::vector< std::list<double> > relord(n+1);
   size_t clusterNumber = 1;
   for (size_t k=0; k<n-1; ++k, ++clusterNumber) {
      double i = merge(k, 0);
      if (i < 0)
         relord[clusterNumber].push_back(-i);
      else
         relord[clusterNumber].splice(relord[clusterNumber].end(), relord[(size_t)i]);

      double j = merge(k, 1);
      if (j < 0)
         relord[clusterNumber].push_back(-j);
      else
         relord[clusterNumber].splice(relord[clusterNumber].end(), relord[(size_t)j]);
   }

   STOPIFNOT(relord[n-1].size() == n);
   size_t k = 0;
   for (std::list<double>::iterator it = relord[n-1].begin();
         it != relord[n-1].end(); ++it) {
      order[k++] = (*it);
   }
}


void HClustResult::generateMergeMatrix()
{
   STOPIFNOT(curiter == n-1);

   std::vector<size_t> elements(n+1, 0);
   std::vector<size_t> parents(n+1, 0);

   size_t clusterNumber = 1;
   for (size_t k=0; k<n-1; ++k, ++clusterNumber) {
      size_t i = (size_t)links(k, 0) + 1;
      size_t j = (size_t)links(k, 1) + 1;
      size_t si = elements[i];
      size_t sj = elements[j];
      elements[i] = clusterNumber;
      elements[j] = clusterNumber;

      if (si == 0)
         merge(k, 0) = -(double)i;
      else {
         while (parents[si] != 0) {
            size_t sinew = parents[si];
            parents[si] = clusterNumber;
            si = sinew;
         }
         if (si != 0) parents[si] = clusterNumber;
         merge(k,0) = (double)si;
      }

      if (sj == 0)
         merge(k, 1) = -(double)j;
      else {
         while (parents[sj] != 0) {
            size_t sjnew = parents[sj];
            parents[sj] = clusterNumber;
            sj = sjnew;
         }
         if (sj != 0) parents[sj] = clusterNumber;
         merge(k,1) = (double)sj;
      }

      if (merge(k, 0) < 0) {
         if (merge(k, 1) < 0 && merge(k, 0) < merge(k, 1)) std::swap(merge(k, 0), merge(k, 1));
      }
      else {
         if (merge(k, 0) > merge(k, 1)) std::swap(merge(k, 0), merge(k, 1));
      }
   }
}


List HClustResult::toR(
      const HClustStats& hclustStats,
      const HClustOptions& hclustOptions,
      const DistanceStats& distStats)
{
   MESSAGE_2("[%010.3f] generating output matrix\n", clock()/(float)CLOCKS_PER_SEC);

   List result = List::create(
      _["merge"]  = merge,
      _["height"] = height,
      _["order"]  = order,
      _["labels"] = labels,
      _["call"]   = R_NilValue,
      _["method"] = R_NilValue,
      _["dist.method"] = dist_method,
      _["links"]  = links,
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
