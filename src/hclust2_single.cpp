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


#ifndef HCLUST2_SINGLE_H_
#define HCLUST2_SINGLE_H_

#define VANTAGE_POINT_SELECT_SCHEME 3
// #define MB_IMPROVEMENT
// #define USE_BOOST_DISJOINT_SETS




#include <Rcpp.h>
#define USE_RINTERNALS
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>
#include <stdio.h>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <fstream>
#include <deque>
#include <exception>
#include <string>
#include <boost/property_map/property_map.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>


#include "hclust2_distance.h"
#include "hclust2_merge.h"
#ifdef USE_BOOST_DISJOINT_SETS
#include <boost/pending/disjoint_sets.hpp>
#else
#include "disjoint_sets.h"
#endif


using namespace Rcpp;
using namespace std;
using namespace boost;


namespace DataStructures{


class HClustSingleBiVpTree
{
protected:

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

   struct Node
   {
      size_t vpindex;
      size_t left;
      size_t right;
      double radius;
      bool sameCluster;
      Node *ll, *lr, *rl, *rr;

      Node() :
         vpindex(SIZE_MAX), left(SIZE_MAX), right(SIZE_MAX), radius(-INFINITY),
         sameCluster(false), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

      Node(size_t left, size_t right) :
         vpindex(SIZE_MAX), left(left), right(right), radius(-INFINITY),
         sameCluster(false), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

      Node(size_t vpindex, double radius) :
         vpindex(vpindex), left(SIZE_MAX), right(SIZE_MAX), radius(radius),
         sameCluster(false), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

      ~Node() {
         if(ll) delete ll;
         if(lr) delete lr;
         if(rl) delete rl;
         if(rr) delete rr;
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

   size_t maxNumberOfElementsInLeaves; // set in the constructor
   const size_t maxNearestNeighborPrefetch = 1;

   Node* _root;
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
   DisjointSets ds;
#endif

#ifdef MB_IMPROVEMENT
   unordered_map<SortedPoint, double> distClust;
   unordered_set<size_t> clusters;
   bool mbimprovement = false;
#endif


   int chooseNewVantagePoint(size_t left, size_t right)
   {
#if VANTAGE_POINT_SELECT_SCHEME == 1
      // idea by A. Fu et al., "Dynamic VP-tree indexing for n-nearest neighbor
      //    search given pair-wise distances"
      size_t numCandidates = 5;
      size_t numTest = 12;

      if (left + numCandidates + numTest > right )
         return left;

      // maximize variance
      size_t bestIndex = -1;
      double bestSigma = -INFINITY;
      for(size_t i=left; i<left+numCandidates; i++) {
         accumulators::accumulator_set< double,
            accumulators::features<accumulators::tag::variance> > acc;
         for (size_t j = left+numCandidates; j < left+numCandidates+numTest; ++j)
            acc( (*_distance)( _indices[i], _indices[j] ) );
         double curSigma = accumulators::variance(acc);
         if (curSigma > bestSigma) {
            bestSigma = curSigma;
            bestIndex = i;
         }
      }

      return bestIndex;
#elif VANTAGE_POINT_SELECT_SCHEME == 2
      // idea by T. Bozkaya and M. Ozsoyoglu, "Indexing large metric spaces
      //      for similarity search queries"

      // which one maximizes dist to _indices[left]?
      size_t bestIndex = left;
      double bestDist  = 0.0;
      for (size_t i=left+1; i<right; ++i) {
         double curDist = (*_distance)(_indices[left], _indices[i]);
         if (curDist > bestDist) {
            bestDist = curDist;
            bestIndex = i;
         }
      }
//       for (size_t i=left+2; i<right; ++i) {
//          double curDist = (*_distance)(_indices[left+1], _indices[i]);
//          if (curDist > bestDist) {
//             bestDist = curDist;
//             bestIndex = i;
//          }
//       }
//       for (size_t i=left+3; i<right; ++i) {
//          double curDist = (*_distance)(_indices[left+2], _indices[i]);
//          if (curDist > bestDist) {
//             bestDist = curDist;
//             bestIndex = i;
//          }
//       }
      return bestIndex;
#else
      // return random index == left one (sample is randomized already)
      return left;
#endif
   }


   Node* buildFromPoints(size_t left, size_t right)
   {
      if(right - left <= maxNumberOfElementsInLeaves)
      {
         for (size_t i=left; i<right; ++i) {
            size_t j = _indices[(i+1 < right)?(i+1):left];
            if (_indices[i] < j)
               maxRadiuses[ _indices[i] ] = (*_distance)(_indices[i], j);
         }

         return new Node(left, right);
      }

      size_t vpi_idx = chooseNewVantagePoint(left, right);
      std::swap(_indices[left], _indices[vpi_idx]);
      size_t vpi = _indices[left];

// #if VERBOSE > 7
//       ((EuclideanDistance*)_distance)->isVP[vpi] = true;
// #endif

      size_t median = ( right + left - 1 ) / 2;
      std::nth_element(_indices.begin() + left + 1, _indices.begin() + median,  _indices.begin() + right,
                       DistanceComparator(vpi, _distance));
      // std::sort(_indices.begin() + left+1, _indices.begin() + right,
                       // DistanceComparator(vpi, &_distance ));
      // printf("(%d,%d,%d)\n", left, median, right);
      // for (int i=left; i<right; ++i) printf("%d, ", _indices[i]+1);
      // printf("\n");
      Node* node = new Node(vpi, (*_distance)(vpi, _indices[median]));


      size_t middle1 = std::partition(_indices.begin() + left,  _indices.begin() + median + 1,  IndexComparator(vpi)) - _indices.begin();
      size_t middle2 = std::partition(_indices.begin() + median + 1,  _indices.begin() + right, IndexComparator(vpi)) - _indices.begin();
      // printf("(%d,%d,%d,%d,%d)\n", left, middle1, median, middle2, right);
      // for (int i=left; i<right; ++i) printf("%d, ", _indices[i]+1);
      // printf("\n");


      if (middle1 - left > 0)     node->ll = buildFromPoints(left, middle1);
      if (median+1 - middle1 > 0) node->lr = buildFromPoints(middle1, median + 1);
      if (middle2 - median-1 > 0) node->rl = buildFromPoints(median + 1, middle2);
      if (right-middle2 > 0)      node->rr = buildFromPoints(middle2, right);

      return node;
   }


   void getNearestNeighborsFromMinRadiusRecursive( Node* node, size_t index,
      size_t clusterIndex, double minR, double& maxR,
      std::priority_queue<HeapNeighborItem>& heap )
   {
      // search within (minR, maxR]
      if (node == NULL) return;

      if (node->sameCluster) {
         if (node->vpindex == SIZE_MAX) {
            if (ds.find_set(_indices[node->left]) == clusterIndex) return;
         } else {
            if (ds.find_set(node->vpindex) == clusterIndex) return;
         }
      }

      if(node->vpindex == SIZE_MAX) // leaf
      {
         if (node->sameCluster) {
#ifdef MB_IMPROVEMENT
         size_t s = SIZE_MAX;
         if(mbimprovement)
            s = ds.find_set(_indices[node->left]);
         std::unordered_map<SortedPoint,double>::const_iterator distToClusterIterator;
         if(mbimprovement)
            distToClusterIterator = distClust.find(SortedPoint(s, clusterIndex));
         double distToCluster = INFINITY;
         if(mbimprovement && distToClusterIterator != distClust.end())
            distToCluster = distToClusterIterator->second;
#endif
            for(size_t i=node->left; i<node->right; i++)
            {
               if(index >= _indices[i]) continue;
               double dist2 = (*_distance)(index, _indices[i]);
               if (dist2 > maxR || dist2 <= minR) continue;

#ifdef MB_IMPROVEMENT
               if(mbimprovement && dist2 > distToCluster) {
                  //Rcout << "odrzucam!" << endl;
                  continue;
               }
#endif

               if (heap.size() >= maxNearestNeighborPrefetch) {
                  if (dist2 < maxR) {
                     while (!heap.empty() && heap.top().dist == maxR) {
                        heap.pop();
                     }
                  }
               }
               heap.push( HeapNeighborItem(_indices[i], dist2) );
               maxR = heap.top().dist;
#ifdef MB_IMPROVEMENT
               if(mbimprovement)
                  distClust.emplace(SortedPoint(s, clusterIndex), dist2);
#endif
            }
         }
         else {
            size_t commonCluster = ds.find_set(_indices[node->left]);
            for(size_t i=node->left; i<node->right; i++)
            {
               size_t currentCluster = ds.find_set(_indices[i]);
#ifdef MB_IMPROVEMENT
               std::unordered_map<SortedPoint,double>::const_iterator distToClusterIterator;
               if(mbimprovement)
                  distToClusterIterator = distClust.find(SortedPoint(currentCluster, clusterIndex));
               double distToCluster = INFINITY;
               if(mbimprovement && distToClusterIterator != distClust.end())
                  distToCluster = distToClusterIterator->second;
#endif
               if (currentCluster != commonCluster) commonCluster = SIZE_MAX;
               if (currentCluster == clusterIndex) continue;

               if (index >= _indices[i]) continue;

               double dist2 = (*_distance)(index, _indices[i]);
               if (dist2 > maxR || dist2 <= minR) continue;
#ifdef MB_IMPROVEMENT
               if(mbimprovement && dist2 > distToCluster) continue;
#endif

               if (heap.size() >= maxNearestNeighborPrefetch) {
                  if (dist2 < maxR) {
                     while (!heap.empty() && heap.top().dist == maxR) {
                        heap.pop();
                     }
                  }
               }
               heap.push( HeapNeighborItem(_indices[i], dist2) );
               maxR = heap.top().dist;
#ifdef MB_IMPROVEMENT
               if(mbimprovement)
                  distClust.emplace(SortedPoint(currentCluster, clusterIndex), dist2);
#endif
            }
            if (commonCluster != SIZE_MAX) node->sameCluster = true;
         }
         return;
      }
      // else // not a leaf
      double dist = (*_distance)(node->vpindex, index);

// this MB's improvement is not well-tested:
//       if (dist < maxR && dist > minR && index < node->vpindex && ds.find_set(node->vpindex) != clusterIndex) {
//
//          heap.push( HeapNeighborItem(node->vpindex, dist) );
//          maxR = heap.top().dist;
//       }


      if ( dist < node->radius ) {
         if ( dist - maxR <= node->radius && dist + node->radius > minR ) {

            if(node->ll && index <= node->vpindex)
               getNearestNeighborsFromMinRadiusRecursive( node->ll, index, clusterIndex, minR, maxR, heap );
            if(node->lr)
               getNearestNeighborsFromMinRadiusRecursive( node->lr, index, clusterIndex, minR, maxR, heap );
         }

         if ( dist + maxR >= node->radius ) {
            if(node->rl && index <= node->vpindex)
               getNearestNeighborsFromMinRadiusRecursive( node->rl, index, clusterIndex, minR, maxR, heap );
            if(node->rr)
               getNearestNeighborsFromMinRadiusRecursive( node->rr, index, clusterIndex, minR, maxR, heap );
         }

      } else /* ( dist >= node->radius ) */ {
         if ( dist + maxR >= node->radius ) {
            if(node->rl && index <= node->vpindex)
               getNearestNeighborsFromMinRadiusRecursive( node->rl, index, clusterIndex, minR, maxR, heap );
            if(node->rr)
               getNearestNeighborsFromMinRadiusRecursive( node->rr, index, clusterIndex, minR, maxR, heap );
         }

         if ( dist - maxR <= node->radius && dist + node->radius > minR ) {
            if(node->ll && index <= node->vpindex)
               getNearestNeighborsFromMinRadiusRecursive( node->ll, index, clusterIndex, minR, maxR, heap );
            if(node->lr)
               getNearestNeighborsFromMinRadiusRecursive( node->lr, index, clusterIndex, minR, maxR, heap );
         }
      }

      if (   !node->sameCluster
         && (!node->ll || node->ll->sameCluster)
         && (!node->lr || node->lr->sameCluster)
         && (!node->rl || node->rl->sameCluster)
         && (!node->rr || node->rr->sameCluster)  )
      {
         size_t commonCluster = SIZE_MAX;
         if (node->ll) {
            size_t currentCluster = ds.find_set((node->ll->vpindex == SIZE_MAX)?_indices[node->ll->left]:node->ll->vpindex);
            if (commonCluster == SIZE_MAX) commonCluster = currentCluster;
            else if (currentCluster != commonCluster) return;
         }
         if (node->lr) {
            size_t currentCluster = ds.find_set((node->lr->vpindex == SIZE_MAX)?_indices[node->lr->left]:node->lr->vpindex);
            if (commonCluster == SIZE_MAX) commonCluster = currentCluster;
            else if (currentCluster != commonCluster) return;
         }
         if (node->rl) {
            size_t currentCluster = ds.find_set((node->rl->vpindex == SIZE_MAX)?_indices[node->rl->left]:node->rl->vpindex);
            if (commonCluster == SIZE_MAX) commonCluster = currentCluster;
            else if (currentCluster != commonCluster) return;
         }
         if (node->rr) {
            size_t currentCluster = ds.find_set((node->rr->vpindex == SIZE_MAX)?_indices[node->rr->left]:node->rr->vpindex);
            if (commonCluster == SIZE_MAX) commonCluster = currentCluster;
            else if (currentCluster != commonCluster) return;
         }
         node->sameCluster = true;
     }
   }

   void print(Node* n) {
      if (n->ll) {
         Rprintf("\"%llx\" -> \"%llx\" [label=\"LL\"];\n", (unsigned long long)n, (unsigned long long)(n->ll));
         print(n->ll);
      }
      if (n->lr) {
         Rprintf("\"%llx\" -> \"%llx\" [label=\"LR\"];\n", (unsigned long long)n, (unsigned long long)(n->lr));
         print(n->lr);
      }
      if (n->rl) {
         Rprintf("\"%llx\" -> \"%llx\" [label=\"RL\"];\n", (unsigned long long)n, (unsigned long long)(n->rl));
         print(n->rl);
      }
      if (n->rr) {
         Rprintf("\"%llx\" -> \"%llx\" [label=\"RR\"];\n", (unsigned long long)n, (unsigned long long)(n->rr));
         print(n->rr);
      }
      if (n->vpindex == SIZE_MAX) {
         for (size_t i=n->left; i<n->right; ++i)
            Rprintf("\"%llx\" -> \"%llu\" [arrowhead = diamond];\n", (unsigned long long)n, (unsigned long long)_indices[i]+1);
      }
      else {
         Rprintf("\"%llx\" [label=\"(%llu, %g)\"];\n", (unsigned long long)n, (unsigned long long)n->vpindex+1, n->radius);
      }
   }


public:

   HeapNeighborItem getNearestNeighbor(size_t index)
   {
#if VERBOSE > 5
      // Rprintf(".");
#endif
      if(shouldFind[index] && nearestNeighbors[index].empty())
      {
         std::priority_queue<HeapNeighborItem> heap;
         size_t clusterIndex = ds.find_set(index);

         double _tau = maxRadiuses[index];
//       THIS IS SLOWER:
//          size_t test = (size_t)(index+unif_rand()*(_n-index));
//          if (ds.find_set(test) != clusterIndex)
//             _tau = (*_distance)(index, test);

         getNearestNeighborsFromMinRadiusRecursive( _root, index, clusterIndex, minRadiuses[index], _tau, heap );
         while( !heap.empty() ) {
            nearestNeighbors[index].push_front(heap.top());
            heap.pop();
         }
         maxRadiuses[index] = INFINITY;
         size_t newNeighborsCount = nearestNeighbors[index].size();

         neighborsCount[index] += newNeighborsCount;
         if(neighborsCount[index] > _n - index || newNeighborsCount == 0)
            shouldFind[index] = false;

         if(newNeighborsCount > 0)
            minRadiuses[index] = nearestNeighbors[index].back().dist;
      }

      if(!nearestNeighbors[index].empty())
      {
         auto res = nearestNeighbors[index].front();
         nearestNeighbors[index].pop_front();
         return res;
      }
      else
      {
         return HeapNeighborItem(SIZE_MAX,-INFINITY);
         //stop("nie ma sasiadow!");
      }
   }

/*
   void searchRadiusKnownIndex(int index, double tau, std::vector<RObject>* results,
                               std::vector<double>* distances, bool findItself = true)
   {
      if(index < 0 || index >= _items.size()) stop("Index out of bounds.");
      std::priority_queue<HeapNeighborItem> heap;

      _tau = tau;
      search( _root, index, false, -1, heap, findItself );

      results->clear(); distances->clear();

      while( !heap.empty() ) {
         results->push_back( _items[heap.top().index] );
         distances->push_back( heap.top().dist );
         heap.pop();
      }

      std::reverse( results->begin(), results->end() );
      std::reverse( distances->begin(), distances->end() );
   }
   */

#ifdef DEBUG
   void printCounters()
   {
      _distance.printCounters();
   }
#endif

public:

   // constructor (OK, we all know what this is, but I label it for faster in-code search)
   HClustSingleBiVpTree(Distance* dist, size_t maxNumberOfElementsInLeaves) :
      maxNumberOfElementsInLeaves(maxNumberOfElementsInLeaves),
      _root(NULL), _n(dist->getObjectCount()), _distance(dist),
      _indices(dist->getObjectCount()),
      neighborsCount(vector<size_t>(dist->getObjectCount(), 0)),
      minRadiuses(vector<double>(dist->getObjectCount(), -INFINITY)),
      maxRadiuses(vector<double>(dist->getObjectCount(), INFINITY)),
      shouldFind(vector<bool>(dist->getObjectCount(), true)),
      nearestNeighbors(vector< deque<HeapNeighborItem> >(dist->getObjectCount())),
#ifdef USE_BOOST_DISJOINT_SETS
      ds(make_assoc_property_map(rank), make_assoc_property_map(parent))
#else
      ds(dist->getObjectCount())
#endif
   {
#if VERBOSE > 5
      Rprintf("[%010.3f] building vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
#endif
      // maxNumberOfElementsInLeaves = 2; //(size_t)log2(_n);

      // starting _indices: random permutation of {0,1,...,_n-1}
      for(size_t i=0;i<_n;i++) _indices[i] = i;
      for(size_t i=_n-1; i>= 1; i--)
         swap(_indices[i], _indices[(size_t)(unif_rand()*(i+1))]);

      for(size_t i=0; i<_n; i++)
      {
#ifdef USE_BOOST_DISJOINT_SETS
        ds.make_set(i);
#endif
#ifdef MB_IMPROVEMENT
        clusters.emplace(i);
#endif
      }

      _root = buildFromPoints(0, _n);
   }


   virtual ~HClustSingleBiVpTree() {
#if VERBOSE > 5
      Rprintf("[%010.3f] destroying vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
#endif
      if(_root) delete _root;
   }


   void print() {
      Rprintf("digraph vptree {\n");
      Rprintf("size=\"6,6\";\n");
	   Rprintf("node [color=lightblue2, style=filled];");
      print(_root);
      Rprintf("}\n");
   }


   NumericMatrix compute()
   {
      NumericMatrix ret(_n-1, 2);
      priority_queue<HeapHierarchicalItem> pq;

      // INIT: Pre-fetch a few nearest neighbors for each point
#if VERBOSE > 5
      Rprintf("[%010.3f] prefetching NNs\n", clock()/(float)CLOCKS_PER_SEC);
#endif
#if VERBOSE > 3
      int misses = 0;
#endif
      for(size_t i=0;i<_n;i++)
      {
         if (i % 1024 == 0) Rcpp::checkUserInterrupt(); // may throw an exception
#if VERBOSE > 7
         if (i % 1024 == 0) Rprintf("\r             prefetch NN: %d/%d", i, _n-1);
#endif
         HeapNeighborItem hi=getNearestNeighbor(i);

         if(hi.index != SIZE_MAX)
         {
            //Rcout <<"dla " << i << "najblizszym jest " << hi->index << endl;
            pq.push(HeapHierarchicalItem(i, hi.index, hi.dist));
         }
      }
#if VERBOSE > 7
      Rprintf("\r             prefetch NN: %d/%d", _n-1, _n-1);
      Rprintf("\n");
#endif
#if VERBOSE > 5
      Rprintf("[%010.3f] merging clusters\n", clock()/(float)CLOCKS_PER_SEC);
#endif

#ifdef MB_IMPROVEMENT
      int nsqrt = (int)sqrt((double)_n);
#endif

      size_t i = 0;
      while(i < _n - 1)
      {
         //Rcout << "iteracja " << i << endl;
         //Rcout << "pq size = " << pq.size()<< endl;
         HeapHierarchicalItem hhi = pq.top();
         pq.pop();

         size_t s1 = ds.find_set(hhi.index1);
         size_t s2 = ds.find_set(hhi.index2);
         if(s1 != s2)
         {
            ret(i,0)=(double)hhi.index1;
            ret(i,1)=(double)hhi.index2;
            ++i;
            ds.link(s1, s2);
#ifdef MB_IMPROVEMENT
            size_t s_new = ds.find_set(s1);
            size_t s_old;
            //Rcout << "przed usuwaniem z clusters" << endl;
            if(s1==s_new)
            {
            	clusters.erase(s2);
            	s_old = s2;
            }
            else
            {
            	clusters.erase(s1);
            	s_old = s1;
            }
            //Rcout << "clusters size = " << clusters.size() << endl;
            if(i >= _n - nsqrt)
            {
               //Rcout << "po sqrt" << endl;
               mbimprovement = true;
               //Rcout << "aktualizuje clusters dist" << endl;
               for ( auto it = clusters.begin(); it != clusters.end(); ++it )
               {
                  if(*it == s_new || *it == s_old) continue;
                  SortedPoint spold = SortedPoint(*it, s_old);
                  SortedPoint spnew = SortedPoint(*it, s_new);
                  auto dold = distClust.find(spold);
                  auto dnew = distClust.find(spnew);
                  if(dold != distClust.end() && dnew != distClust.end())
                  {
                     //Rcout << "znaleziono obie" << endl;
                     dnew->second = min(dnew->second, dold->second);
                     //Rcout << "po aktu" << endl;
                  }
                  else if(dold != distClust.end())
                  {
                     //Rcout << "znaleziono tylko stara" << endl;
                     distClust.emplace(spnew, dold->second);
                     //Rcout << "po aktu2" << endl;
                  }

                  if(dold != distClust.end())
                  {
                     //Rcout << "znaleziono stara" << endl;
                     distClust.erase(spold);
                     //Rcout << "po aktu3" << endl;
                  }
               }
               //Rcout << "po aktualizacji clusters dist" << endl;
            }
#endif
            if (i % 10000 == 0) Rcpp::checkUserInterrupt(); // may throw an exception
         }
#if VERBOSE > 3
         else
            ++misses;
#endif
#if VERBOSE > 7
         if (i % 1024 == 0) Rprintf("\r             %d / %d / %d ", i+1, _n, misses);
#endif

         // ASSERT: hhi.index1 < hhi.index2
         HeapNeighborItem hi=getNearestNeighbor(hhi.index1);
         if(hi.index != SIZE_MAX)
            pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));
      }
#if VERBOSE > 7
      Rprintf("\r             %d / %d / %d ", _n, _n, misses);
#endif
#if VERBOSE > 3
      Rprintf("Total ignored NNs: %d\n", misses);
#endif
#if VERBOSE > 5
      Rprintf("[%010.3f] generating output matrix\n", clock()/(float)CLOCKS_PER_SEC);
#endif
      Rcpp::checkUserInterrupt();

      MergeMatrixGenerator mmg(ret.nrow());
      return mmg.generateMergeMatrix(ret);
   }

}; // class

} // namespace DataStructures


// [[Rcpp::export(".hclust2_single")]]
RObject hclust2_single(RObject distance, RObject objects, int maxNumberOfElementsInLeaves=2) {
#if VERBOSE > 5
   Rprintf("[%010.3f] starting timer\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   RObject result(R_NilValue);
   DataStructures::Distance* dist = DataStructures::Distance::createDistance(distance, objects);

   try {
      /* Rcpp::checkUserInterrupt(); may throw an exception */
      DataStructures::HClustSingleBiVpTree hclust(dist, (int)maxNumberOfElementsInLeaves);
      RObject merge = hclust.compute();
      result = Rcpp::as<RObject>(List::create(
         _["merge"]  = merge,
         _["height"] = R_NilValue,
         _["order"]  = R_NilValue,
         _["labels"] = R_NilValue,
         _["call"]   = R_NilValue,
         _["method"] = "single",
         _["dist.method"] = R_NilValue
      ));
      result.attr("class") = "hclust";
      //hclust.print();
   }
   catch(...) {

   }

   if (dist) delete dist;
#if VERBOSE > 5
   Rprintf("[%010.3f] done\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   if (Rf_isNull(result)) stop("stopping on error or explicit user interrupt");
   return result;
}


#endif /* VPTREEBINHIERARCHICAL_H_ */
