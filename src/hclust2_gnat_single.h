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

#ifndef __HCLUST2_GNAT_SINGLE_H
#define __HCLUST2_GNAT_SINGLE_H

// #define GNAT_DEBUG


// ************************************************************************


// ************************************************************************

#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <queue>
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

struct HClustGnatRange
{
   double min;
   double max;

   HClustGnatRange()
      : min(-INFINITY), max(INFINITY)   {   }

   HClustGnatRange(double min, double max)
      : min(min), max(max)   {   }
};

struct HClustGnatSingleNode
{
   vector<size_t> splitPoints;
   //size_t splitPointIndex;
   size_t left;
   size_t right;
   size_t degree; // splitPoints.size; ==0 => leaf
   size_t maxindex;
   bool sameCluster;
   vector<HClustGnatSingleNode *> children;
   Matrix<HClustGnatRange> splitPointsRanges;

   HClustGnatSingleNode() :
      left(SIZE_MAX), right(SIZE_MAX), degree(SIZE_MAX), maxindex(SIZE_MAX),
      sameCluster(false), splitPointsRanges()  {
      }

   HClustGnatSingleNode(size_t left, size_t right) :
      left(left), right(right), degree(SIZE_MAX), maxindex(SIZE_MAX),
      sameCluster(false), splitPointsRanges()  {
      }

   HClustGnatSingleNode(size_t vpindex) :
      left(SIZE_MAX), right(SIZE_MAX), degree(SIZE_MAX), maxindex(SIZE_MAX),
      sameCluster(false), splitPointsRanges()  {
      }

   ~HClustGnatSingleNode() {
      for (size_t i = 0; i<children.size(); ++i)
      {
         if (children[i]) //it should always be true
            delete children[i];
      }
   }
};

class HClustGnatSingle
{
protected:

   HClustTreeOptions opts;

   HClustGnatSingleNode* _root;
   size_t _n;
   Distance* _distance;
   std::vector<size_t> _indices;

   std::vector<size_t> neighborsCount;
   std::vector<double> minRadiuses;
   // std::vector<double> maxRadiuses;
   std::vector<bool> shouldFind;
   std::vector< deque<HeapNeighborItem> > nearestNeighbors;

   std::map<size_t,size_t> rank;
   std::map<size_t,size_t> parent;

   HClustTreeStats stats;
   PhatDisjointSets ds;
   bool prefetch;

   // unordered_map<Point, HClustGnatRange> splitPointsRanges; //to oznacza, ze dla Point(i,j) dostajemy range(p_i, D_pj), szczegoly w artykule, niesymetryczne!

   void chooseNewSplitPoints(HClustGnatSingleNode *node, size_t degree, size_t left, size_t right);
   vector<size_t> groupPointsToSplitPoints(HClustGnatSingleNode *node, size_t left, size_t right);
   HClustGnatSingleNode* buildFromPoints(size_t degree, size_t optdegree,  size_t left, size_t right);
   HClustGnatSingleNode* createNonLeafNode(size_t degree, size_t optdegree,  size_t left, size_t right);
   vector<size_t> chooseDegrees(size_t degree, size_t optdegree,  size_t left, size_t allPointsCount, const vector<size_t>& boundaries);

   void printIndices();
   void excludeRegions(HClustGnatSingleNode* node, vector<bool>& shouldFind, vector<double>& distances, double minR, double& maxR);
   size_t degreeDecreaser(size_t optdegree);

   void getNearestNeighborsFromMinRadiusRecursive(HClustGnatSingleNode* node,
      size_t index, size_t clusterIndex, double minR, double& maxR,
      std::priority_queue<HeapNeighborItem>& heap);

public:

   void print(HClustGnatSingleNode* n);
   void print();
   void FindNeighborTest(size_t index, double R);

   HClustGnatSingle(Distance* dist, RObject control);
   ~HClustGnatSingle();

   NumericMatrix compute();

   HeapNeighborItem getNearestNeighbor(size_t index);

   inline const HClustTreeStats& getStats() { return stats; }
   inline const HClustTreeOptions& getOptions() { return opts; }

}; // class

} // namespace DataStructures


#endif
