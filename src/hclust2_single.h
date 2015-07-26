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

#ifndef __HCLUST2_SINGLE_H
#define __HCLUST2_SINGLE_H



// ************************************************************************

// #define MB_IMPROVEMENT
// #define USE_BOOST_DISJOINT_SETS


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
#ifdef USE_BOOST_DISJOINT_SETS
#include <boost/pending/disjoint_sets.hpp>
#else
#include "disjoint_sets.h"
#endif

#ifdef MB_IMPROVEMENT
#include <unordered_set>
#endif // MB_IMPROVEMENT



namespace DataStructures{

class HClustBiVpTreeSingle
{
protected:

   HClustBiVpTreeOptions opts;

   HClustBiVpTreeNode* _root;
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

   HClustBiVpTreeStats stats;

#ifdef USE_BOOST_DISJOINT_SETS
   boost::disjoint_sets<
     associative_property_map< std::map<size_t,size_t> >,
     associative_property_map< std::map<size_t,size_t> > > ds;
#else
   PhatDisjointSets ds;
#endif

#ifdef MB_IMPROVEMENT
   std::unordered_map<SortedPoint, double> distClust;
   std::unordered_set<size_t> clusters;
   bool mbimprovement = false;
#endif  // MB_IMPROVEMENT


   int chooseNewVantagePoint(size_t left, size_t right);
   HClustBiVpTreeNode* buildFromPoints(size_t left, size_t right);

   void getNearestNeighborsFromMinRadiusRecursive(HClustBiVpTreeNode* node,
      size_t index, size_t clusterIndex, double minR, double& maxR,
      std::priority_queue<HeapNeighborItem>& heap);

   void print(HClustBiVpTreeNode* n);


public:

   HClustBiVpTreeSingle(Distance* dist, RObject control);
   ~HClustBiVpTreeSingle();

   void print();
   NumericMatrix compute();

   HeapNeighborItem getNearestNeighbor(size_t index);

   inline const HClustBiVpTreeStats& getStats() { return stats; }
   inline const HClustBiVpTreeOptions& getOptions() { return opts; }

}; // class

} // namespace DataStructures


#endif
