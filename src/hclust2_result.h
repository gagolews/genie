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

#ifndef __HCLUST2_RESULT_H
#define __HCLUST2_RESULT_H

#include "defs.h"
#include "disjoint_sets.h"
#include "hclust2_common.h"
#include "hclust2_distance.h"
#include <Rcpp.h>



namespace grup
{

class HClustResult
{
private:
   size_t curiter;
   size_t n;

   Rcpp::NumericMatrix links;
   Rcpp::NumericMatrix merge;
   Rcpp::NumericVector height;
   Rcpp::NumericVector order;
   // call is set by R
   // method is set by R
   Rcpp::RObject labels;
   Rcpp::RObject dist_method;
   bool lite;

   void generateMergeMatrix();
   void generateOrderVector();


public:
   HClustResult(size_t n, Distance* dist, bool lite=false);

   Rcpp::NumericMatrix getLinks() { return links; }
   Rcpp::NumericVector getHeight() { return height; }

   void link(size_t i1, size_t i2, double d12);

   Rcpp::List toR(
         const grup::HClustStats& hclustStats,
         const grup::HClustOptions& hclustOptions,
         const grup::DistanceStats& distStats
   );

}; // struct HClustResult



} // namespace grup

#endif
