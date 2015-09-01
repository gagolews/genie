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

#ifndef __HCLUST2_VPTREE_SINGLE_APPROX_H
#define __HCLUST2_VPTREE_SINGLE_APPROX_H

#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
// #include <fstream>
#include <deque>
// #include <exception>
// #include <string>
// #include <boost/property_map/property_map.hpp>
// #include <boost/tuple/tuple_comparison.hpp>
#include <algorithm>

#include "hclust2_common.h"
#include "disjoint_sets.h"
#include "hclust2_result.h"

using namespace std;
using namespace Rcpp;

// ************************************************************************

namespace DataStructures
{

struct HClustVpTreeSingleNodeApprox
{
   size_t vpindex;
   size_t left;
   size_t right;
   double radius;
   bool sameCluster;
   size_t maxindex;
   HClustVpTreeSingleNodeApprox* childL;
   HClustVpTreeSingleNodeApprox* childR;

   HClustVpTreeSingleNodeApprox() :
         vpindex(SIZE_MAX), left(SIZE_MAX), right(SIZE_MAX), radius(-INFINITY),
         sameCluster(false), maxindex(SIZE_MAX), childL(NULL), childR(NULL)  { }

   HClustVpTreeSingleNodeApprox(size_t left, size_t right) :
         vpindex(SIZE_MAX), left(left), right(right), radius(-INFINITY),
         sameCluster(false), maxindex(SIZE_MAX), childL(NULL), childR(NULL)  { }

   HClustVpTreeSingleNodeApprox(size_t vpindex, size_t left, size_t right, double radius) :
         vpindex(vpindex), left(left), right(right), radius(radius),
         sameCluster(false), maxindex(SIZE_MAX), childL(NULL), childR(NULL)  { }

   ~HClustVpTreeSingleNodeApprox() {
      if (childL) delete childL;
      if (childR) delete childR;
   }
};


class HClustVpTreeSingleApprox
{
protected:
   HClustOptions& opts;
   HClustStats& stats;
   bool& prefetch;
   DisjointSets& ds;
   size_t& n;
   Distance* distance;
   std::vector<size_t>& indices;

   HClustVpTreeSingleNodeApprox* root;
   // bool visitAll; // for testing only

   size_t chooseNewVantagePoint(size_t left, size_t right);
   HClustVpTreeSingleNodeApprox* buildFromPoints(size_t left, size_t right, std::vector<double>& distances);

   void getNearestNeighborsFromMinRadiusRecursive(HClustVpTreeSingleNodeApprox* node,
      size_t index, size_t clusterIndex, double minR, std::priority_queue<double>& bestR, double& maxR, NNHeap& nnheap);
   void getNearestNeighborsFromMinRadiusRecursiveLeaf(HClustVpTreeSingleNodeApprox* node,
      size_t index, size_t clusterIndex, double minR, std::priority_queue<double>& bestR, double& maxR, NNHeap& nnheap);
   void getNearestNeighborsFromMinRadiusRecursiveNonLeaf(HClustVpTreeSingleNodeApprox* node,
      size_t index, size_t clusterIndex, double minR, std::priority_queue<double>& bestR, double& maxR, NNHeap& nnheap);

   void updateSameClusterFlag(HClustVpTreeSingleNodeApprox* node);

   void print(HClustVpTreeSingleNodeApprox* node);

public:

   HClustVpTreeSingleApprox(HClustOptions& opts, HClustStats& stats, bool& prefetch, DisjointSets& ds, size_t& n, Distance* dist, std::vector<size_t>& indices);
   ~HClustVpTreeSingleApprox();

   void build();
   void getNearestNeighborsFromMinRadius(size_t index, size_t clusterIndex, double minR, NNHeap& nnheap);
   void print();

}; // class

} // namespace DataStructures


#endif
