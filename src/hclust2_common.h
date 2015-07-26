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


#include "hclust2_distance.h"


namespace DataStructures{

struct HeapNeighborItem {
   size_t index;
   double dist;

   HeapNeighborItem(size_t index, double dist) :
      index(index), dist(dist) {}

   HeapNeighborItem() :
      index(SIZE_MAX), dist(-INFINITY) {}

   bool operator<( const HeapNeighborItem& o ) const {
      return dist < o.dist;
   }
};


struct HeapHierarchicalItem {
   size_t index1;
   size_t index2;
   double dist;

   HeapHierarchicalItem(size_t index1, size_t index2, double dist) :
      index1(index1), index2(index2), dist(dist) {}

   bool operator<( const HeapHierarchicalItem& o ) const {
      return dist >= o.dist;
   }
};


struct DistanceComparator
{
   size_t index;
   Distance* distance;

   DistanceComparator(size_t index, Distance* distance )
      : index(index), distance(distance) {}

   bool operator()(size_t a, size_t b) {
      return (*distance)( index, a ) < (*distance)( index, b );
   }
};


struct IndexComparator
{
   size_t index;

   IndexComparator(size_t index)
      : index(index) {}

   bool operator()(size_t a) {
      return a <= index;
   }
};


struct HClustBiVpTreeNode
{
   size_t vpindex;
   size_t left;
   size_t right;
   double radius;
   bool sameCluster;
   HClustBiVpTreeNode *ll, *lr, *rl, *rr;

   HClustBiVpTreeNode() :
      vpindex(SIZE_MAX), left(SIZE_MAX), right(SIZE_MAX), radius(-INFINITY),
      sameCluster(false), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

   HClustBiVpTreeNode(size_t left, size_t right) :
      vpindex(SIZE_MAX), left(left), right(right), radius(-INFINITY),
      sameCluster(false), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

   HClustBiVpTreeNode(size_t vpindex, double radius) :
      vpindex(vpindex), left(SIZE_MAX), right(SIZE_MAX), radius(radius),
      sameCluster(false), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

   ~HClustBiVpTreeNode() {
      if(ll) delete ll;
      if(lr) delete lr;
      if(rl) delete rl;
      if(rr) delete rr;
   }
};


struct HClustBiVpTreeStats {
   size_t nodeCount; // now many nodes are there in the tree
   size_t nodeVisit; // now many nodes were visited during NN search
   size_t nnCals;    // how many times NN search job was launched
   size_t nnCount;   // how many NNs were obtained in overall

   HClustBiVpTreeStats() :
      nodeCount(0), nodeVisit(0), nnCals(0), nnCount(0) {}

   ~HClustBiVpTreeStats() {
      #if VERBOSE > 0
      Rprintf("             vp-tree: nodeCount=%.0f, nodeVisit=%.0f\n",
         (double)nodeCount, (double)nodeVisit);
      Rprintf("             vp-tree: nnCals=%.0f, nnCount=%.0f\n",
         (double)nnCals, (double)nnCount);
      #endif
   }
};


struct SortedPoint
{
   size_t i;
   size_t j;

   SortedPoint()
      :i(0),j(0) {}

   SortedPoint(size_t _i, size_t _j)
   {
      if(_j < _i)
      {
         i = _j;
         j = _i;
      }
      else
      {
         i = _i;
         j = _j;
      }
   }

   bool operator==(const SortedPoint &other) const
   {
      return (i == other.i && j == other.j);
   }
};

} // namespace DataStructures


namespace std {

   template <>
      struct hash<DataStructures::SortedPoint>
   {
      std::size_t operator()(const DataStructures::SortedPoint& k) const
      {
        std::size_t seed = 0;
        boost::hash_combine(seed, k.i);
        boost::hash_combine(seed, k.j);
        return seed;
      }
   };
} // namespace std



#endif
