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
#include <exception>
#include <unordered_map>
#include <string>

#include <boost/functional/hash.hpp>

using namespace Rcpp;
using namespace std;

struct RFunction2 {
   RFunction2(const Function& _f) : f(_f) {
      R_PreserveObject(f);
   }

   ~RFunction2() {
      R_ReleaseObject(f);
   }

   Function f;
};

class SortedPoint
{
   public:
   int i;
   int j;
   SortedPoint():i(0),j(0){}
   SortedPoint(int _i, int _j)
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


namespace std {

   template <>
      struct hash<SortedPoint>
   {
      std::size_t operator()(const SortedPoint& k) const
      {
        std::size_t seed = 0;
        boost::hash_combine(seed, k.i);
        boost::hash_combine(seed, k.j);
        return seed;
      }
   };
}

struct distClass2
{
   RFunction2* distance;
   vector<RObject> *items;

   unordered_map<SortedPoint, double> hashmap;

   distClass2(RFunction2* distance, vector<RObject> *items) : distance(distance), items(items)
   {
   #ifdef DEBUG
    metricCalculated=0;
    hashmapHit=0;
   #endif
   }

   #ifdef DEBUG
   int metricCalculated;
   int hashmapHit;
   #endif

   double operator()(int v1, int v2)
   {
#ifdef DEBUG
      metricCalculated++;
#endif
      if(v1==v2) return 0;
      SortedPoint p(v1,v2);
      std::unordered_map<SortedPoint,double>::iterator got = hashmap.find(p);
      if ( got == hashmap.end() )
      {
         NumericVector res = distance->f((*items)[v1],(*items)[v2]);
         double d = res[0];
         hashmap.emplace(p, d);
         return d;
      }
      else
      {
#ifdef DEBUG
         got->second.incrementCounter();
         hashmapHit++;
#endif
         //Rcout << "trafilem"<<endl;
         return got->second;
      }
   }

#ifdef DEBUG
   void printCounters()
   {
      for (auto it=hashmap.begin(); it != hashmap.end(); ++it)
         Rcout << it->first << " => " << it->second.counter << endl;
      Rcout << "hashmap count = " << hashmap.size() << endl;
   }
#endif
};


class VpTreeHierarchical2
{
   public:

   VpTreeHierarchical2(vector<RObject>* items, RFunction2* rf)
    : _root(NULL), _distance(rf, items), _items(items), _indices(items->size()), _n(items->size())
   {
      for(size_t i=0;i<_n;i++)
         _indices[i] = i;

      for(int i=_n-1; i>= 1; i--)
      {
         int j = (int)(unif_rand()*(i+1));
         swap(_indices[i], _indices[j]);
      }

      _root = buildFromPoints(0, _n);
      #ifdef DEBUG
      Rcout << "metric calculated = " << _distance.metricCalculated << endl;
      Rcout << "hashmapHit = " << _distance.hashmapHit << endl;
      #endif
   }



   virtual ~VpTreeHierarchical2() {
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
/*
   virtual void searchKNNKnownIndex(int index, int k, std::vector<RObject>* results,
                            std::vector<double>* distances, bool findItself = true)
   {
      if(index < 0 || index >= _items->size()) stop("Index out of bounds.");
      #ifdef DEBUG
      _distance.metricCalculated = 0;
      _distance.hashmapHit = 0;
      #endif
      std::priority_queue<HeapItem> heap;

      _tau = std::numeric_limits<double>::max();
      search( _root, index, true, k, heap, findItself );

      results->clear(); distances->clear();

      while( !heap.empty() ) {
         results->push_back( _items[heap.top().index] );
         distances->push_back( heap.top().dist );
         heap.pop();
      }

      std::reverse( results->begin(), results->end() );
      std::reverse( distances->begin(), distances->end() );
      #ifdef DEBUG
      Rcout << "metric calculated = " << _distance.metricCalculated << endl;
      Rcout << "hashmapHit = " << _distance.hashmapHit << endl;
      #endif
   }

   void searchRadiusKnownIndex(int index, double tau, std::vector<RObject>* results,
                               std::vector<double>* distances, bool findItself = true)
   {
      if(index < 0 || index >= _items.size()) stop("Index out of bounds.");
      std::priority_queue<HeapItem> heap;

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

   protected:
   distClass2 _distance;
   std::vector<RObject>* _items;
   double _tau;
   std::vector<int> _indices;
   size_t _n;

   struct Node
   {
      int vpindex;
      double radius;
      int left;
      int right;
      Node *ll, *lr, *rl, *rr;

      Node() :
      vpindex(-1), left(-1), right(-1), radius(-1), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

      Node(int left, int right) :
      vpindex(-1), left(left), right(right), radius(-1), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

      Node(int vpindex, double radius) :
      vpindex(vpindex), left(-1), right(-1), radius(radius), ll(NULL), lr(NULL), rl(NULL), rr(NULL) {}

      ~Node() {
         if(ll) delete ll;
         if(lr) delete lr;
         if(rl) delete rl;
         if(rr) delete rr;
      }
   }* _root;



   struct HeapItem {
      HeapItem( int index, double dist) :
      index(index), dist(dist) {}
      int index;
      double dist;
      bool operator<( const HeapItem& o ) const {
         return dist < o.dist;
      }
   };

   struct DistanceComparator
   {
      int index;
      distClass2* distance;

      DistanceComparator(int index, distClass2* distance ) : index(index), distance(distance) {}
      bool operator()(int a, int b) {
         return (*distance)( index, a ) < (*distance)( index, b );
      }
   };

   struct IndexComparator
   {
      int index;

      IndexComparator(int index ) : index(index) {}
      bool operator()(int a) {
         return a < index;
      }
   };

   Node* buildFromPoints(int left, int right)
   {
      const int maxNumberOfElementInLeaf = 4;
      if(right - left < maxNumberOfElementInLeaf)
      {
         return new Node(left, right);
      }

      int vpi = left;//(int)((double)rand() / RAND_MAX * (upper - lower - 1) ) + lower;

      int median = ( right + left ) / 2;
      std::nth_element(_indices.begin() + left + 1, _indices.begin() + median,  _indices.begin() + right,
                       DistanceComparator(vpi, &_distance ));
      Node* node = new Node(vpi, _distance(_indices[left], _indices[median]));


      int middle1 = std::partition(_indices.begin() + left,  _indices.begin() + median + 1,  IndexComparator(vpi)) - _indices.begin();
      int middle2 = std::partition(_indices.begin() + median + 1,  _indices.begin() + right, IndexComparator(vpi)) - _indices.begin();

      node->ll = buildFromPoints(left, middle1);
      node->lr = buildFromPoints(middle1, median + 1);
      node->rl = buildFromPoints(median + 1, middle2);
      node->rr = buildFromPoints(middle2, right);

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
   /*
   void search( Node* node, const RObject& target, bool isKNN, int k,
               std::priority_queue<HeapItem>& heap, bool findItself)
   {
      if ( node == NULL ) return;


      //printf("dist=%g tau=%gn", dist, _tau );
      if(node->isLeaf)
      {
         for(size_t i=0;i<node->points.size();i++)
         {
            double dist2 = _distance( _items[node->points[i]], target );
            if(dist2 < 1e-6)
            {
                if(!findItself && R_compute_identical(_items[node->points[i]], target, 16))
                   continue;
            }
            if ( (dist2 < _tau && isKNN) || (dist2 <= _tau && !isKNN) ) {
               if ( heap.size() == (size_t)k && isKNN) heap.pop();
               heap.push( HeapItem(node->points[i], dist2) );
               if ( heap.size() ==(size_t) k && isKNN)
               {
                  _tau = heap.top().dist;
                  //Rcout << "current tau=" << _tau << endl;
               }
            }
         }
      }
      else
      {
         double dist = _distance( _items[node->vpindex], target );
         vector<bool> visited(node->childCount, false);

         for(int i=0;i<node->childCount-1;i++)
         {
            //if ( dist < node->radiuses[i] ) {
            if ( dist - _tau <= node->radiuses[i] && !visited[i]) {
               search( node->children[i], target, isKNN, k, heap, findItself );
               visited[i]=true;
            }
            if ( dist + _tau >= node->radiuses[i] && !visited[i+1]) {
               search( node->children[i+1], target, isKNN, k, heap, findItself );
               visited[i+1]=true;
            }
            //}
         }
         /*if ( dist + _tau >= node->radiuses[node->childCount-2] ) {
				search( node->children[node->childCount-1], target, k, heap );
			}
      }
   }

   virtual void search( Node* node, int index, bool isKNN, int k,
               std::priority_queue<HeapItem>& heap, bool findItself )
   {
      if ( node == NULL ) return;

      //printf("dist=%g tau=%gn", dist, _tau );
      if(node->isLeaf)
      {
         for(size_t i=0;i<node->points.size();i++)
         {
            double dist2 = _distance(node->points[i], index );
            if(dist2 < 1e-6)
            {
                if(!findItself && R_compute_identical(_items[node->points[i]], _items[index], 16))
                   continue;
            }
            if ( (dist2 < _tau && isKNN) || (dist2 <= _tau && !isKNN) ) {
               if ( heap.size() ==(size_t) k && isKNN) heap.pop();
               heap.push( HeapItem(node->points[i], dist2) );
               if ( heap.size() == (size_t) k && isKNN)
               {
                  _tau = heap.top().dist;
                  //Rcout << "current tau=" << _tau << endl;
               }
            }
         }
      }
      else
      {
         double dist = _distance(node->vpindex, index);
         vector<bool> visited(node->childCount, false);

         for(int i=0;i<node->childCount-1;i++)
         {
            //if ( dist < node->radiuses[i] ) {
            if ( dist - _tau <= node->radiuses[i] && !visited[i]) {
               search( node->children[i], index, isKNN, k, heap, findItself );
               visited[i]=true;
            }
            if ( dist + _tau >= node->radiuses[i] && !visited[i+1]) {
               search( node->children[i+1], index, isKNN, k, heap, findItself );
               visited[i+1]=true;
            }
            //}
         }

      }
   }
   */
};
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


// [[Rcpp::export]]
IntegerMatrix hclust2(Function distance, List listobj) { //https://code.google.com/p/vptree/source/browse/src/vptree/VpTreeNode.java

   RFunction2 *rf = new RFunction2(distance);
   vector<RObject> points(listobj.begin(), listobj.end());
   // prze


   VpTreeHierarchical2 _tree(&points, rf);
   //IntegerMatrix im = _tree->hierarchicalClustering();
   delete rf;
   return IntegerMatrix(0,0);
}


#endif /* VPTREEBINHIERARCHICAL_H_ */
