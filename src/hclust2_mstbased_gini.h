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

#ifndef __HCLUST2_MSTBASED_GINI_H
#define __HCLUST2_MSTBASED_GINI_H



// ************************************************************************


#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <deque>
#include "hclust2_common.h"
#include "disjoint_sets.h"
#include "hclust2_result.h"

namespace grup
{



class HClustMSTbasedGini
{
protected:

   HClustOptions* opts;
   size_t n;
   HClustStats stats;
   Distance* distance;

   HclustPriorityQueue getMST();
   void linkAndRecomputeGini(PhatDisjointSets& ds, double& lastGini, size_t s1, size_t s2);

public:

   HClustMSTbasedGini(Distance* dist, HClustOptions* opts);
   virtual ~HClustMSTbasedGini();

   HClustResult compute();

   inline const HClustStats& getStats()     { return stats; }
   inline const HClustOptions& getOptions() { return *opts; }

}; // class

} // namespace grup


#endif
