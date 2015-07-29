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

#ifndef __HCLUST2_COMMON_H
#define __HCLUST2_COMMON_H

// ---------------------------------------------------------------------------

#define DEFAULT_MAX_LEAVES_ELEMS 4
#define DEFAULT_MAX_NN_PREFETCH 2
#define DEFAULT_VP_SELECT_SCHEME 3
#define DEFAULT_VP_SELECT_CAND 5
#define DEFAULT_VP_SELECT_TEST 12

// ---------------------------------------------------------------------------

#include "hclust2_distance.h"


namespace DataStructures
{

struct HeapNeighborItem
{
   size_t index;
   double dist;

   HeapNeighborItem(size_t index, double dist) :
      index(index), dist(dist) {}

   HeapNeighborItem() :
      index(SIZE_MAX), dist(-INFINITY) {}

   inline bool operator<( const HeapNeighborItem& o ) const {
      return dist < o.dist;
   }
};


struct HeapHierarchicalItem
{
   size_t index1;
   size_t index2;
   double dist;

   HeapHierarchicalItem(size_t index1, size_t index2, double dist) :
      index1(index1), index2(index2), dist(dist) {}

   inline bool operator<( const HeapHierarchicalItem& o ) const {
      return dist >= o.dist;
   }
};


struct DistanceComparator
{
   size_t index;
   Distance* distance;

   DistanceComparator(size_t index, Distance* distance)
      : index(index), distance(distance) {}

   inline bool operator()(size_t a, size_t b) {
      return (*distance)( index, a ) < (*distance)( index, b );
   }
};


struct IndexComparator
{
   size_t index;

   IndexComparator(size_t index)
      : index(index) {}

   inline bool operator()(size_t a) {
      return (a <= index);
   }
};


struct HClustBiVpTreeOptions
{
   size_t maxLeavesElems;
   size_t maxNNPrefetch;
   size_t vpSelectScheme;
   size_t vpSelectCand;   // for vpSelectScheme == 1
   size_t vpSelectTest;   // for vpSelectScheme == 1

   HClustBiVpTreeOptions(Rcpp::RObject control);
   Rcpp::NumericVector toR() const;
};


struct HClustBiVpTreeStats
{
   size_t nodeCount; // how many nodes are there in the tree
   size_t leafCount; // how many leaves
   size_t nodeVisit; // how many nodes were visited during NN search
   size_t nnCals;    // how many times NN search job was launched
   size_t nnCount;   // how many NNs were obtained in overall

   HClustBiVpTreeStats();
   ~HClustBiVpTreeStats();
   Rcpp::NumericVector toR() const;
};


struct SortedPoint
{
   size_t i;
   size_t j;

   SortedPoint()
      :i(0),j(0) {}

   SortedPoint(size_t _i, size_t _j) {
      if(_j < _i) {
         i = _j;
         j = _i;
      }
      else {
         i = _i;
         j = _j;
      }
   }

   inline bool operator==(const SortedPoint &other) const {
      return (i == other.i && j == other.j);
   }
};

} // namespace DataStructures


namespace std
{

template <> struct hash<DataStructures::SortedPoint>
{
   std::size_t operator()(const DataStructures::SortedPoint& k) const {
     std::size_t seed = 0;
     boost::hash_combine(seed, k.i);
     boost::hash_combine(seed, k.j);
     return seed;
   }
};

} // namespace std



#endif
