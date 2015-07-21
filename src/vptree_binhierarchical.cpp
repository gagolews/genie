#ifndef VPTREEBINHIERARCHICAL_H_
#define VPTREEBINHIERARCHICAL_H_
#include <iostream>
#include <Rcpp.h>
#define USE_RINTERNALS
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <queue>
#include <limits>
#include <numeric>
#include <fstream>
#include <deque>
#include <exception>
#include <unordered_map>
#include <string>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/functional/hash.hpp>

using namespace Rcpp;
using namespace std;
using namespace boost;

#define HASHMAP_DISABLE
// #define HASHMAP_COUNTERS
#define VERBOSE 0

namespace DataStructures{
namespace HClustSingleBiVpTree{

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

} // namespace HClustSingleBiVpTree
} // namespace DataStructures




namespace std {

   template <>
      struct hash<DataStructures::HClustSingleBiVpTree::SortedPoint>
   {
      std::size_t operator()(const DataStructures::HClustSingleBiVpTree::SortedPoint& k) const
      {
        std::size_t seed = 0;
        boost::hash_combine(seed, k.i);
        boost::hash_combine(seed, k.j);
        return seed;
      }
   };
}


namespace DataStructures{
namespace HClustSingleBiVpTree{


struct Distance {
   size_t n;

   Distance(size_t n) :
      n(n) { }

   virtual ~Distance() { }

   inline size_t getObjectCount() { return n; }

   virtual double operator()(size_t v1, size_t v2) = 0;
};

struct EuclideanDistance : public Distance
{
   SEXP robj1;
#ifndef HASHMAP_DISABLE
   vector< unordered_map<size_t, double> > hashmap;
#endif
   const double* items;
   size_t m;
#ifdef HASHMAP_COUNTERS
   size_t hashmapHit;
   size_t hashmapMiss;
#endif

   EuclideanDistance(const NumericMatrix& points) :
      Distance(points.nrow()), robj1(points),
#ifndef HASHMAP_DISABLE
      hashmap(vector< unordered_map<size_t, double> >(points.nrow())),
#endif
      items(REAL((SEXP)points)), m(points.ncol())
   {
      R_PreserveObject(robj1);

#ifdef HASHMAP_COUNTERS
      hashmapHit=0;
      hashmapMiss=0;
#endif
   }

#ifdef HASHMAP_COUNTERS
   virtual ~EuclideanDistance()
   {
#if VERBOSE > 5
      Rprintf("[%010.3f] destroying distance object\n", clock()/(float)CLOCKS_PER_SEC);
#endif
#if VERBOSE > 3
      Rprintf("Distance Hashmap: #hits=%d, #miss=%d, est.mem.used>=%.1fMB (vs %.1fMB)\n",
         hashmapHit, hashmapMiss, 8.0f*hashmapMiss/1000.0f/1000.0f,
         8.0f*(n-1)*(n-1)*0.5f/1000.0f/1000.0f);
#endif
      R_ReleaseObject(robj1);
   }
#endif

   virtual double operator()(size_t v1, size_t v2)
   {
      if (v1 == v2) return 0.0;
#ifndef HASHMAP_DISABLE
      if (v1 > v2) std::swap(v1, v2);

      // std::unordered_map<SortedPoint,double>::iterator got = hashmap.find(p);
      auto got = hashmap[v1].find(v2);
      if ( got == hashmap[v1].end() )
      {
#endif
         double d = 0.0;
         for (size_t i=0; i<m; ++i)
            d += (items[v1+i*n]-items[v2+i*n])*(items[v1+i*n]-items[v2+i*n]);
         d = sqrt(d);
#ifdef HASHMAP_COUNTERS
         ++hashmapMiss;
#endif
#ifndef HASHMAP_DISABLE
         hashmap[v1].emplace(v2, d);
#endif
         return d;
#ifndef HASHMAP_DISABLE
      }
      else
      {
#ifdef HASHMAP_COUNTERS
         ++hashmapHit;
#endif
         return got->second;
      }
#endif
   }
};

struct GenericRDistance : public Distance
{
   Function distfun;
#ifndef HASHMAP_DISABLE
   vector< unordered_map<size_t, double> > hashmap;
#endif
   vector<RObject> items;
#ifdef HASHMAP_COUNTERS
   size_t hashmapHit;
   size_t hashmapMiss;
#endif

   GenericRDistance(const Function& distfun, const vector<RObject>& items) :
      Distance(items.size()), distfun(distfun),
#ifndef HASHMAP_DISABLE
      hashmap(vector< unordered_map<size_t, double> >(items.size())),
#endif
      items(items)
   {
      R_PreserveObject(distfun);
#ifdef HASHMAP_COUNTERS
      hashmapHit=0;
      hashmapMiss=0;
#endif
   }

#ifdef HASHMAP_COUNTERS
   virtual ~GenericRDistance()
   {
#if VERBOSE > 5
      Rprintf("[%010.3f] destroying distance object\n", clock()/(float)CLOCKS_PER_SEC);
#endif
#if VERBOSE > 3
      Rprintf("Distance Hashmap: #hits=%d, #miss=%d, est.mem.used>=%.1fMB (vs %.1fMB)\n",
         hashmapHit, hashmapMiss, 8.0f*hashmapMiss/1000.0f/1000.0f,
         8.0f*(items.size()-1)*(items.size()-1)*0.5f/1000.0f/1000.0f);
#endif
      R_ReleaseObject(distfun);
   }
#endif

   virtual double operator()(size_t v1, size_t v2)
   {
      if (v1 == v2) return 0;
#ifndef HASHMAP_DISABLE
      SortedPoint p(v1,v2);
      // std::unordered_map<SortedPoint,double>::iterator got = hashmap.find(p);
      auto got = hashmap[p.i].find(p.j);
      if ( got == hashmap[p.i].end() )
      {
#endif
#ifdef HASHMAP_COUNTERS
         ++hashmapMiss;
#endif
         double d = ((NumericVector)distfun(items[v1],items[v2]))[0];
#ifndef HASHMAP_DISABLE
         hashmap[p.i].emplace(p.j, d);
#endif
         return d;
#ifndef HASHMAP_DISABLE
      }
      else
      {
#ifdef HASHMAP_COUNTERS
         ++hashmapHit;
#endif
         return got->second;
      }
#endif
   }
};


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

   boost::disjoint_sets<
     associative_property_map< std::map<size_t,size_t> >,
     associative_property_map< std::map<size_t,size_t> > > ds;

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

      size_t vpi = _indices[left];//(int)((double)rand() / RAND_MAX * (upper - lower - 1) ) + lower;

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



   /*

   size_t calculateNodeSize(Node* node)
   {
      return node->radiuses.size()*sizeof(double)
         + node->points.size()*sizeof(int)
         + node->children.size()*sizeof(Node*)
         + sizeof(Node);
   }

   size_t treeSize_rec(Node* node)
   {
      size_t size = calculateNodeSize(node);
      for(int i=0;i<node->childCount;i++)
      {
         size += treeSize_rec(node->children[i]);
      }
      return size;
   }

   int treeHeight_rec(Node* node)
   {
      int maxH = 0;
      for(int i=0;i<node->childCount;i++)
      {
         maxH = max(treeHeight_rec(node->children[i]), maxH);
      }
      return maxH+1;
   }
*/

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
            for(size_t i=node->left; i<node->right; i++)
            {
               if(index >= _indices[i]) continue;
               double dist2 = (*_distance)(index, _indices[i]);
               if (dist2 > maxR || dist2 <= minR) continue;
               if (heap.size() >= maxNearestNeighborPrefetch) {
                  if (dist2 < maxR) {
                     while (!heap.empty() && heap.top().dist == maxR) {
                        heap.pop();
                     }
                  }
               }
               heap.push( HeapNeighborItem(_indices[i], dist2) );
               maxR = heap.top().dist;
            }
         }
         else {
            size_t commonCluster = ds.find_set(_indices[node->left]);
            for(size_t i=node->left; i<node->right; i++)
            {
               size_t currentCluster = ds.find_set(_indices[i]);
               if (currentCluster != commonCluster) commonCluster = SIZE_MAX;
               if (currentCluster == clusterIndex) continue;

               if (index >= _indices[i]) continue;

               double dist2 = (*_distance)(index, _indices[i]);
               if (dist2 > maxR || dist2 <= minR) continue;

               if (heap.size() >= maxNearestNeighborPrefetch) {
                  if (dist2 < maxR) {
                     while (!heap.empty() && heap.top().dist == maxR) {
                        heap.pop();
                     }
                  }
               }
               heap.push( HeapNeighborItem(_indices[i], dist2) );
               maxR = heap.top().dist;
            }
            if (commonCluster != SIZE_MAX) node->sameCluster = true;
         }
         return;
      }
      // else // not a leaf
      double dist = (*_distance)(node->vpindex, index);
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
    static NumericMatrix generateMergeMatrix(const NumericMatrix& x) {
      // x has 0-based indices
      size_t n = x.nrow();
      if (x.ncol() != 2) stop("x should have 2 columns");
      NumericMatrix y(n, 2);

      /* -------------------------------------------------------------- */
      /* TO DO: new method, O(n)                                        */




      /* -------------------------------------------------------------- */

      std::vector< std::unordered_set<size_t> > curclust(n);
      std::vector< bool > alreadyInSomeCluster(n+1, false);

      for (size_t k=0; k<n; ++k) {
         if (k % 1000 == 0) Rcpp::checkUserInterrupt(); // may throw an exception
         size_t i = (size_t)x(k,0)+1;
         size_t j = (size_t)x(k,1)+1;
         size_t si = (k > 0) ? k-1 : SIZE_MAX;
         size_t sj = (k > 0) ? k-1 : SIZE_MAX;
         if (alreadyInSomeCluster[i])
            while (si != SIZE_MAX && curclust[si].find(i) == curclust[si].end())
               si = (si>0) ? si-1 : SIZE_MAX;
         else si = SIZE_MAX;

         if (alreadyInSomeCluster[j])
            while (sj != SIZE_MAX && curclust[sj].find(j) == curclust[sj].end())
               sj = (sj>0) ? sj-1 : SIZE_MAX;
         else sj = SIZE_MAX;

         if (si == SIZE_MAX && sj == SIZE_MAX) {
            curclust[k].insert(i);
            curclust[k].insert(j);
            y(k,0) = -(double)i;
            y(k,1) = -(double)j;
            alreadyInSomeCluster[i] = true;
            alreadyInSomeCluster[j] = true;
         }
         else if (si == SIZE_MAX && sj != SIZE_MAX) {
            curclust[k].insert(curclust[sj].begin(), curclust[sj].end());
            curclust[k].insert(i);
            curclust[sj].clear(); // no longer needed
            y(k,0) = -(double)i;
            y(k,1) = (double)sj+1;
            alreadyInSomeCluster[i] = true;
         }
         else if (si != SIZE_MAX && sj == SIZE_MAX) {
            curclust[k].insert(curclust[si].begin(), curclust[si].end());
            curclust[k].insert(j);
            curclust[si].clear(); // no longer needed
            y(k,0) = (double)si+1;
            y(k,1) = -(double)j;
            alreadyInSomeCluster[j] = true;
         }
         else {
            curclust[k].insert(curclust[si].begin(), curclust[si].end());
            curclust[k].insert(curclust[sj].begin(), curclust[sj].end());
            curclust[si].clear(); // no longer needed
            curclust[sj].clear(); // no longer needed
            y(k,0) = (double)si+1;
            y(k,1) = (double)sj+1;
         }
      }
      return y;
   }

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
   HClustSingleBiVpTree(Distance* dist) :
      _root(NULL), _n(dist->getObjectCount()), _distance(dist),
      _indices(dist->getObjectCount()),
      neighborsCount(vector<size_t>(dist->getObjectCount(), 0)),
      minRadiuses(vector<double>(dist->getObjectCount(), -INFINITY)),
      maxRadiuses(vector<double>(dist->getObjectCount(), INFINITY)),
      shouldFind(vector<bool>(dist->getObjectCount(), true)),
      nearestNeighbors(vector< deque<HeapNeighborItem> >(dist->getObjectCount())),
      ds(make_assoc_property_map(rank), make_assoc_property_map(parent))
   {
#if VERBOSE > 5
      Rprintf("[%010.3f] building vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
#endif
      maxNumberOfElementsInLeaves = 12; //(size_t)log2(_n);

      // starting _indices: random permutation of {0,1,...,_n-1}
      for(size_t i=0;i<_n;i++) _indices[i] = i;
      for(size_t i=_n-1; i>= 1; i--)
         swap(_indices[i], _indices[(size_t)(unif_rand()*(i+1))]);

      for(size_t i=0; i<_n; i++)
        ds.make_set(i);

      _root = buildFromPoints(0, _n);
   }


   virtual ~HClustSingleBiVpTree() {
#if VERBOSE > 5
      Rprintf("[%010.3f] destroying vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
#endif
      if(_root) delete _root;
   }


   /*size_t treeSize()
   {
      if(_root==NULL) return sizeof(VpTree);
      return sizeof(VpTree) + treeSize_rec(_root);
   }

   int treeHeight()
   {
      if(_root==NULL) return 0;
      return treeHeight_rec(_root);
   }*/

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
         if (i % 1000 == 0) Rcpp::checkUserInterrupt(); // may throw an exception
#if VERBOSE > 7
      Rprintf("\r             prefetch NN: %d/%d", i, _n-1);
#endif
         HeapNeighborItem hi=getNearestNeighbor(i);

         if(hi.index != SIZE_MAX)
         {
            //Rcout <<"dla " << i << "najblizszym jest " << hi->index << endl;
            pq.push(HeapHierarchicalItem(i, hi.index, hi.dist));
         }
      }
#if VERBOSE > 7
      Rprintf("\n");
#endif
#if VERBOSE > 5
      Rprintf("[%010.3f] merging clusters\n", clock()/(float)CLOCKS_PER_SEC);
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
            if (i % 10000 == 0) Rcpp::checkUserInterrupt(); // may throw an exception
         }
#if VERBOSE > 3
         else
            ++misses;
#endif
#if VERBOSE > 7
         Rprintf("\r             %d / %d / %d ", i+1, _n, misses);
#endif

         // ASSERT: hhi.index1 < hhi.index2
         HeapNeighborItem hi=getNearestNeighbor(hhi.index1);
         if(hi.index != SIZE_MAX)
            pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));
      }
#if VERBOSE > 3
      Rprintf("Total ignored NNs: %d\n", misses);
#endif
#if VERBOSE > 5
   Rprintf("[%010.3f] generating output matrix\n", clock()/(float)CLOCKS_PER_SEC);
#endif
      Rcpp::checkUserInterrupt();
      return generateMergeMatrix(ret);
   }

}; // class






/*
template <>
   void vptree<RObject>::findIndex(const RObject& target) // specialize only one member
   {
      Rcout << "specialized" << endl;
      for(int i = 0; i<_items.size(); i++)
      {
         if(Rcpp::all(_items[i] == target))
            //if(_items[i] == target)
            return i;
      }
      stop("There is no such element in the tree.");
   }
*/


} // namespace HClustSingleBiVpTree
} // namespace DataStructures


NumericMatrix transpose(const NumericMatrix& matrix)
{
   size_t width = matrix.ncol();
   size_t height = matrix.nrow();
   NumericMatrix transposed(width, height);

   for (size_t i = 0; i < width; i++)
   {
      for (size_t j = 0; j < height; j++)
      {
         transposed(i,j) = matrix(j,i);
      }
   }

   return transposed;
}


// [[Rcpp::export]]
SEXP hclust2(RObject objects, RObject distance=R_NilValue) {
#if VERBOSE > 5
   Rprintf("[%010.3f] starting timer\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   SEXP result;
   DataStructures::HClustSingleBiVpTree::Distance* dist = NULL;
   try {
      /* Rcpp::checkUserInterrupt(); may throw an exception */
      if (Rf_isVectorList(objects) && Rf_isFunction(distance)) {
         Function distance2(distance);
         List objects2(objects);
         dist = (DataStructures::HClustSingleBiVpTree::Distance*)
            new DataStructures::HClustSingleBiVpTree::GenericRDistance(
               distance2,
               vector<RObject>(objects2.begin(), objects2.end())
            );
      }
      else if (Rf_isMatrix(objects)) {
         NumericMatrix objects2(objects);
         dist = (DataStructures::HClustSingleBiVpTree::Distance*)
            new DataStructures::HClustSingleBiVpTree::EuclideanDistance(
               objects2
            );
      }

      if (!dist) stop("Incorrect input data");
      DataStructures::HClustSingleBiVpTree::HClustSingleBiVpTree hclust(dist);
      result = (SEXP)hclust.compute();
      //hclust.print();
   }
   catch(...) {
      result = R_NilValue;
   }
   if (dist) delete dist;
#if VERBOSE > 5
   Rprintf("[%010.3f] done\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   if (Rf_isNull(result))
      stop("error or stop");
   return result;
}

// [[Rcpp::export]]
NumericMatrix generateMergeMatrix(NumericMatrix x) {
   return DataStructures::HClustSingleBiVpTree::HClustSingleBiVpTree::generateMergeMatrix(x);
}


#endif /* VPTREEBINHIERARCHICAL_H_ */
