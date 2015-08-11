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

#ifndef __HCLUST2_VPTREE_MEDOID_H
#define __HCLUST2_VPTREE_MEDOID_H

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

struct HClustBiVpTreeMedoidNode
{
   size_t vpindex;
   size_t left;
   size_t right;
   double radius;
   bool sameCluster;
   size_t maxindex;
   HClustBiVpTreeMedoidNode* childL;
   HClustBiVpTreeMedoidNode* childR;

   HClustBiVpTreeMedoidNode() :
         vpindex(SIZE_MAX), left(SIZE_MAX), right(SIZE_MAX), radius(-INFINITY),
         sameCluster(false), maxindex(SIZE_MAX), childL(NULL), childR(NULL)  { }

   HClustBiVpTreeMedoidNode(size_t left, size_t right) :
         vpindex(SIZE_MAX), left(left), right(right), radius(-INFINITY),
         sameCluster(false), maxindex(SIZE_MAX), childL(NULL), childR(NULL)  { }

   HClustBiVpTreeMedoidNode(size_t vpindex, size_t left, size_t right, double radius) :
         vpindex(vpindex), left(left), right(right), radius(radius),
         sameCluster(false), maxindex(SIZE_MAX), childL(NULL), childR(NULL)  { }

   ~HClustBiVpTreeMedoidNode() {
      if (childL) delete childL;
      if (childR) delete childR;
   }
};


class HClustBiVpTreeMedoid
{
protected:
   struct HeapHierarchicalItemMedoid {
         size_t index1;
         size_t index2;
         double dist;
         size_t iter;

         HeapHierarchicalItemMedoid(size_t index1, size_t index2, double dist, size_t iter) :
            index1(index1), index2(index2), dist(dist), iter(iter) {}

         bool operator<( const HeapHierarchicalItemMedoid& o ) const {
            return dist >= o.dist;
         }
      };

   struct KKItem {
            double dist;
            size_t iter;

            KKItem(double dist, size_t iter) :
               dist(dist), iter(iter) {}
            KKItem() :
               dist(INFINITY), iter(0) {}
         };

   HClustTreeOptions opts;

   HClustBiVpTreeMedoidNode* _root;
   size_t _n;
   Distance* _distance;
   std::vector<size_t> _indices;
   //std::vector<size_t> _indicesinv;

   std::vector<size_t> neighborsCount;
   std::vector<double> minRadiuses;
   // std::vector<double> maxRadiuses;
   std::vector<bool> shouldFind;
   std::vector< deque<HeapNeighborItem> > nearestNeighbors;
   std::vector<double> distances;

   vector<size_t> medoids;
   vector<bool> medoidFound;

   HClustTreeStats stats;

   PhatDisjointSets ds;
   bool prefetch;

   size_t chooseNewVantagePoint(size_t left, size_t right);
   HClustBiVpTreeMedoidNode* buildFromPoints(size_t left, size_t right);

   void getNearestNeighborsFromMinRadiusRecursive(HClustBiVpTreeMedoidNode* node,
      size_t index, size_t clusterIndex, double minR, double& maxR, NNHeap& nnheap);

   void print(HClustBiVpTreeMedoidNode* n);
   //HeapHierarchicalItemMedoid calculateCluster2ClusterMedoidDistance(size_t item1, size_t item2, size_t iter);
   size_t medoidForCluster(size_t s); 

public:

   HClustBiVpTreeMedoid(Distance* dist, RObject control);
   ~HClustBiVpTreeMedoid();

   void print();
   NumericMatrix compute();

   HeapNeighborItem getNearestNeighbor(size_t index);

   inline const HClustTreeStats& getStats() { return stats; }
   inline const HClustTreeOptions& getOptions() { return opts; }

}; // class

} // namespace DataStructures


#endif
