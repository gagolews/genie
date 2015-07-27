#ifndef __HCLUST2_COMPLETE_H
#define __HCLUST2_COMPLETE_H
// ************************************************************************

// #define MB_IMPROVEMENT
//#define USE_BOOST_DISJOINT_SETS


// ************************************************************************
#define VANTAGE_POINT_SELECT_SCHEME 2


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
#include <boost/property_map/property_map.hpp>
#else
#include "disjoint_sets.h"
#endif

#ifdef MB_IMPROVEMENT
#include <unordered_set>
#endif // MB_IMPROVEMENT
namespace DataStructures
{


class HClustBiVpTreeComplete
{
protected:

   struct HeapHierarchicalItemMax {
         size_t index1;
         size_t index2;
         double dist;
         size_t iter;

         HeapHierarchicalItemMax(size_t index1, size_t index2, double dist, size_t iter) :
            index1(index1), index2(index2), dist(dist), iter(iter) {}

         bool operator<( const HeapHierarchicalItemMax& o ) const {
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

   size_t maxNumberOfElementsInLeaves; // set in the constructor
   const size_t maxNearestNeighborPrefetch = 1;

   HClustBiVpTreeNode* _root;
   size_t _n;
   Distance* _distance;
   std::vector<size_t> _indices;

   std::vector<size_t> neighborsCount;
   std::vector<double> minRadiuses;
   std::vector<double> maxRadiuses;
   std::vector<bool> shouldFind;
   std::vector< deque<HeapNeighborItem> > nearestNeighbors;

   std::map<size_t,size_t> rank;
   std::map<size_t,size_t> parent;

#ifdef USE_BOOST_DISJOINT_SETS
   boost::disjoint_sets<
     associative_property_map< std::map<size_t,size_t> >,
     associative_property_map< std::map<size_t,size_t> > > ds;
#else
   PhatDisjointSets ds;
#endif
#ifdef MB_IMPROVEMENT
   unordered_map<SortedPoint, double> distClust;
   unordered_set<size_t> clusters;
   bool mbimprovement = false;
#endif


   int chooseNewVantagePoint(size_t left, size_t right);


   HClustBiVpTreeNode* buildFromPoints(size_t left, size_t right);

   void getNearestNeighborsFromMinRadiusRecursive( HClustBiVpTreeNode* node, size_t index,
      size_t clusterIndex, double minR, double& maxR,
      std::priority_queue<HeapNeighborItem>& heap );

   void print(HClustBiVpTreeNode* n);

public:
   HeapNeighborItem getNearestNeighbor(size_t index);

#ifdef DEBUG
   void printCounters();
#endif

public:

   // constructor (OK, we all know what this is, but I label it for faster in-code search)
   HClustBiVpTreeComplete(Distance* dist, size_t maxNumberOfElementsInLeaves);
   virtual ~HClustBiVpTreeComplete();
   void print();
   size_t clusterCount(size_t cluster);
   HeapHierarchicalItemMax calculateCluster2ClusterMaxDistance(size_t item1, size_t item2, size_t iter);
   NumericMatrix compute();


}; // class
}//namespace
#endif
