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

#ifndef __HCLUST2_NNBASED_SINGLE_H
#define __HCLUST2_NNBASED_SINGLE_H



// ************************************************************************


#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
// #include <fstream>
// #include <deque>
// #include <exception>
// #include <string>
// #include <boost/property_map/property_map.hpp>
// #include <boost/tuple/tuple_comparison.hpp>
// #include <algorithm>

#include "hclust2_common.h"
#include "hclust2_merge.h"
#include "disjoint_sets.h"


namespace DataStructures
{



class HClustNNbasedSingle
{
protected:

   HClustTreeOptions opts;
   size_t _n;
   Distance* _distance;
   std::vector<size_t> _indices;
   // std::vector<size_t> _indicesinv;

   std::vector<size_t> neighborsCount;
   std::vector<double> minRadiuses;
   // std::vector<double> maxRadiuses;
   std::vector<bool> shouldFind;
   std::vector< deque<HeapNeighborItem> > nearestNeighbors;

   HClustTreeStats stats;

   DisjointSets ds;
   bool prefetch;

   virtual HeapNeighborItem getNearestNeighbor(size_t index, double distMax=INFINITY) = 0;

public:

   HClustNNbasedSingle(Distance* dist, RObject control);
   virtual ~HClustNNbasedSingle();

   virtual void print() { Rcout << "this print method is a stub" << std::endl; }

   NumericMatrix compute();

   inline const HClustTreeStats& getStats() { return stats; }
   inline const HClustTreeOptions& getOptions() { return opts; }

}; // class

} // namespace DataStructures


#endif
