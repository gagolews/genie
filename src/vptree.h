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

#ifndef __VPTREE_H
#define __VPTREE_H



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

struct VpTreeNode
{
   size_t vpindex;
   size_t left;
   size_t right;
   double radius;
   VpTreeNode* childL;
   VpTreeNode* childR;

   VpTreeNode() :
         vpindex(SIZE_MAX), left(SIZE_MAX), right(SIZE_MAX), radius(-INFINITY),
         childL(NULL), childR(NULL)  { }

   VpTreeNode(size_t left, size_t right) :
         vpindex(SIZE_MAX), left(left), right(right), radius(-INFINITY),
         childL(NULL), childR(NULL)  { }

   VpTreeNode(size_t vpindex, size_t left, size_t right, double radius) :
         vpindex(vpindex), left(left), right(right), radius(radius),
         childL(NULL), childR(NULL)  { }

   ~VpTreeNode() {
      if (childL) delete childL;
      if (childR) delete childR;
   }
};


class VpTree
{
protected:

   HClustOptions opts;

   VpTreeNode* _root;
   size_t _n;
   Distance* _distance;
   std::vector<size_t> _indices;
   std::vector<size_t> _indicesinv;

   std::vector<double> distances;

   HClustStats stats;

   size_t chooseNewVantagePoint(size_t left, size_t right);
   VpTreeNode* buildFromPoints(size_t left, size_t right);

   void getNearestNeighborsFromMinRadiusRecursive(
      VpTreeNode* node, size_t index, double minR, double& maxR,
      std::priority_queue<HeapNeighborItem>& nnheap, size_t maxNNPrefetch);

public:

   VpTree(Distance* dist, RObject control);
   ~VpTree();

   vector<HeapNeighborItem> getNearestNeighbors(size_t index, int k, double minR=-INFINITY, double maxR=INFINITY);

   inline const HClustStats& getStats() { return stats; }
   inline const HClustOptions& getOptions() { return opts; }

}; // class

} // namespace DataStructures


#endif
