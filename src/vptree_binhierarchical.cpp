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
#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/functional/hash.hpp>

using namespace Rcpp;
using namespace std;
using namespace boost;

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

   struct HeapItem {
      HeapItem( int index, double dist) :
      index(index), dist(dist) {}
      HeapItem():index(-1), dist(0) {}
      int index;
      double dist;
      bool operator<( const HeapItem& o ) const {
         return dist < o.dist;
      }
   };

   struct HeapHierarchicalItem {
      HeapHierarchicalItem( int index1, int index2, double dist) :
      index1(index1), index2(index2), dist(dist) {}
      int index1;
      int index2;
      double dist;
      bool operator<( const HeapHierarchicalItem& o ) const {
         return dist >= o.dist;
      }
   };

   IntegerMatrix hierarchicalClustering()
   {
      //Rcout << "wszedl do dobrego" << endl;
      IntegerMatrix ret(_n-1, 2);
      neighborsCount = vector<int>(_n, 0);
      minRadiuses = vector<double>(_n, 0);
      shouldFind = vector<bool>(this->_items->size(), true);
      nearestNeighbors = vector<queue<HeapItem>>(_n);
      //Rcout << "dociagam sasiadow. po raz pierwszy.." << endl;
      priority_queue<HeapHierarchicalItem> pq;
      //Rcout << "_items.size() = " << this->_items.size() << endl;
      for(int i=0;i<_n;i++)
      {
         //Rcout << i << endl;
         //stop("kazdemu na poczatek znajduje sasiada");
         HeapItem hi=getNearestNeighbor(i);

         if(hi.index != -1)
         {
            //Rcout <<"dla " << i << "najblizszym jest " << hi->index << endl;
            pq.push(HeapHierarchicalItem(i, hi.index, hi.dist));
         }
         //stop("po pierwszym wstepnym");
      }
      //stop("po wstepnym zebraniu sasiadow");
      std::map<int,int> rank;
      std::map<int,int> parent;

      boost::disjoint_sets<
         associative_property_map<std::map<int,int>>,
         associative_property_map<std::map<int,int>> > ds(
            make_assoc_property_map(rank),
            make_assoc_property_map(parent));

      for(int i=0; i<_n; i++)
          ds.make_set(i);
      //stop("po stworzeniu union find");
      //Rcout << "pq size = " << pq.size() << endl;
      //stop(std::to_string(pq.size()));
      int i = 0;
      while(i < _n - 1)
      //for(int i=0;i<this->_items.size() - 1 ; i++)
      {
         //Rcout << "iteracja " << i << endl;
         //Rcout << "pq size = " << pq.size()<< endl;
         HeapHierarchicalItem hhi = pq.top();
         pq.pop();

         int s1 = ds.find_set(hhi.index1);
         int s2 = ds.find_set(hhi.index2);
         if(s1 != s2)
         {
            ret(i,0)=hhi.index1;
            ret(i,1)=hhi.index2;
            //Rcout << "el1="<<ret(i,0)<< "el2=" <<ret(i,1)<< ", i =" << i << endl;
            i++;
            ds.link(s1, s2);
            //Rcout << "el1="<<hhi.index1+1<< "el2=" <<hhi.index2 +1<< endl;

         }
         //stop("przed dociaganiem sasiadow");

         HeapItem hi=getNearestNeighbor(hhi.index1);
         if(hi.index != -1)
            pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));

         //hi=getNearestNeighbor(hhi.index2);
         //if(hi.index != -1)
         //   pq.push(HeapHierarchicalItem(hhi.index2, hi.index, hi.dist));
         //stop("po pierwszej iteracji");
      }
      return hclust_merge_matrix(ret);
   }

   IntegerMatrix hclust_merge_matrix(const IntegerMatrix x) const {
      // x has 0-based indices
      int n = x.nrow();
      if (x.ncol() != 2) stop("x should have 2 columns");

      IntegerMatrix y(n, 2);
      std::vector< std::unordered_set<int> > curclust(n);

      for (int k=0; k<n; ++k) {
         int i = x(k,0)+1;
         int j = x(k,1)+1;
         int si=k-1, sj=k-1;
         while (si >= 0 && curclust[si].find(i) == curclust[si].end()) si--;
         while (sj >= 0 && curclust[sj].find(j) == curclust[sj].end()) sj--;
         if (si < 0 && sj < 0) {
            curclust[k].insert(i);
            curclust[k].insert(j);
            y(k,0) = -i;
            y(k,1) = -j;
         }
         else if (si < 0 && sj >= 0) {
            curclust[k].insert(curclust[sj].begin(), curclust[sj].end());
            curclust[k].insert(i);
            curclust[sj].clear(); // no longer needed
            y(k,0) = -i;
            y(k,1) = sj+1;
         }
         else if (si >= 0 && sj < 0) {
            curclust[k].insert(curclust[si].begin(), curclust[si].end());
            curclust[k].insert(j);
            curclust[si].clear(); // no longer needed
            y(k,0) = si+1;
            y(k,1) = -j;
         }
         else { // if (si >= 0 && sj >= 0)
            curclust[k].insert(curclust[si].begin(), curclust[si].end());
            curclust[k].insert(curclust[sj].begin(), curclust[sj].end());
            curclust[si].clear(); // no longer needed
            curclust[sj].clear(); // no longer needed
            y(k,0) = si+1;
            y(k,1) = sj+1;
         }
      }
      return y;
   }

   HeapItem getNearestNeighbor(int index)
   {
      const int delta = 10;
      //Rcout << "nearestNeighbors[index] = " << nearestNeighbors[index].size() << endl;
      if(nearestNeighbors[index].empty() && shouldFind[index])
      {
         //Rcout << "kolejka pusta, trzeba dociagnac sasiadow..." << endl;

         priority_queue<HeapItem> heap;
         vector<HeapItem> res;
         if(shouldFind[index])
            res=searchKNNKnownIndex(index, neighborsCount[index]+delta, minRadiuses[index]);
         if(neighborsCount[index]+delta > _n - index)
            shouldFind[index] = false;

         if(res.size() > 0)
            minRadiuses[index] = res[res.size()-1].dist;

         //Rcout << "index=" << index << endl;
         for(int i=neighborsCount[index]; i < min(neighborsCount[index]+delta, (int)res.size()); i++)
         {
            //Rcout << "results[i]" << results[i] << endl;
            //Rcout << "distances[i]" << distances[i] << endl;
            nearestNeighbors[index].push(HeapItem(res[i].index, res[i].dist ));
         }
         //Rcout << "dociagnalem " << heap.size() << "sasiadow" << endl;
         //nearestNeighbors[index] = heap;
         neighborsCount[index] = min(neighborsCount[index]+delta, (int)_n);
         //Rcout << "nearestNeighbors[index] = " << nearestNeighbors[index].size() << endl;
      }

      if(!nearestNeighbors[index].empty())
      {
         HeapItem hiret = nearestNeighbors[index].front();
         nearestNeighbors[index].pop();
         return hiret;
      }
      else
      {
         return HeapItem(-1,-1);
         //stop("nie ma sasiadow!");
      }
   }

   vector<HeapItem> searchKNNKnownIndex(int index, int k, double minR)
   {
      if(index < 0 || index >= _items->size()) stop("Index out of bounds.");
      #ifdef DEBUG
      _distance.metricCalculated = 0;
      _distance.hashmapHit = 0;
      #endif

      std::priority_queue<HeapItem> heap;
      _tau = std::numeric_limits<double>::max();
      search( _root, index, true, k, minR, heap );
      vector<HeapItem> results(heap.size());
      int i=results.size() - 1;

      while( !heap.empty() ) {
         results[i] = heap.top();
         //Rcout << "index = " << results[i].index << ", dist = " << results[i].dist << endl;
         i--;
         heap.pop();
      }
      #ifdef DEBUG
      Rcout << "metric calculated = " << _distance.metricCalculated << endl;
      Rcout << "hashmapHit = " << _distance.hashmapHit << endl;
      #endif
      return results;
   }
/*
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
   std::vector<bool> shouldFind;
   std::vector<int> neighborsCount;
   std::vector<double> minRadiuses;
   std::vector<queue<HeapItem>> nearestNeighbors;

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
         return a <= index;
      }
   };

   Node* buildFromPoints(int left, int right)
   {
      const int maxNumberOfElementInLeaf = 4;
      if(right - left <= maxNumberOfElementInLeaf)
      {
//          printf("(%d,%d)\n", left, right);
//          for (int i=left; i<right; ++i) printf("%d, ", _indices[i]+1);
//          printf("\n");
         return new Node(left, right);
      }

      int vpi = _indices[left];//(int)((double)rand() / RAND_MAX * (upper - lower - 1) ) + lower;

      int median = ( right + left - 1 ) / 2;
      std::nth_element(_indices.begin() + left + 1, _indices.begin() + median,  _indices.begin() + right,
                       DistanceComparator(vpi, &_distance ));
      // std::sort(_indices.begin() + left+1, _indices.begin() + right,
                       // DistanceComparator(vpi, &_distance ));
      // printf("(%d,%d,%d)\n", left, median, right);
      // for (int i=left; i<right; ++i) printf("%d, ", _indices[i]+1);
      // printf("\n");
      Node* node = new Node(vpi, _distance(vpi, _indices[median]));


      int middle1 = std::partition(_indices.begin() + left,  _indices.begin() + median + 1,  IndexComparator(vpi)) - _indices.begin();
      int middle2 = std::partition(_indices.begin() + median + 1,  _indices.begin() + right, IndexComparator(vpi)) - _indices.begin();
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

   virtual void search( Node* node, int index, bool isKNN, int k, double minR,
               std::priority_queue<HeapItem>& heap )
   {
      if ( node == NULL ) return;

      if(node->vpindex == -1)
      {
         for(size_t i=node->left;i<node->right;i++)
         {
            if(index < _indices[i])
            {
               double dist2 = _distance(_indices[i], index );
               if ( (dist2 < _tau && isKNN && dist2 > minR) || (dist2 <= _tau && !isKNN) )
               {

                     if ( heap.size() >=(size_t) k && isKNN) heap.pop();
                     heap.push( HeapItem(_indices[i], dist2) );
                     if ( heap.size() == (size_t) k && isKNN)
                     {
                        _tau = heap.top().dist;
                        //Rcout << "current tau=" << _tau << endl;
                     }

               }
            }
         }
      }
      else
      {


         /*if ( node->ll == NULL && node->lr == NULL && node->rl == NULL && node->rr == NULL ) {
            return;
         }*/
         double dist = _distance(node->vpindex, index);
         if ( dist < node->radius ) {
            if ( dist - _tau <= node->radius && dist + node->radius >= minR ) {

               if(node->ll != NULL && index <= node->vpindex)
                  search( node->ll, index, isKNN, k, minR, heap );
               if(node->lr != NULL)
                  search( node->lr, index, isKNN, k, minR, heap );
            }

            if ( dist + _tau >= node->radius ) {
               if(node->rl && index <= node->vpindex)
                  search( node->rl, index, isKNN, k, minR, heap );
               if(node->rr != NULL)
                  search( node->rr, index, isKNN, k, minR, heap );
            }

         } else {
            if ( dist + _tau >= node->radius ) {
               if(node->rl && index <= node->vpindex)
                  search( node->rl, index, isKNN, k, minR, heap );
               if(node->rr != NULL)
                  search( node->rr, index, isKNN, k, minR, heap );
            }

            if ( dist - _tau <= node->radius && dist + node->radius >= minR) {
               if(node->ll != NULL && index <= node->vpindex)
                  search( node->ll, index, isKNN, k, minR, heap );
               if(node->lr != NULL)
                  search( node->lr, index, isKNN, k, minR, heap );
            }
         }

      }
   }



public:
   void print() {
      Rprintf("digraph vptree {\n");
      Rprintf("size=\"6,6\";\n");
	   Rprintf("node [color=lightblue2, style=filled];");
      print(_root);
      Rprintf("}\n");
   }


protected:
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
      if (n->vpindex < 0) {
         for (int i=n->left; i<n->right; ++i)
            Rprintf("\"%llx\" -> \"%d\" [arrowhead = diamond];\n", (unsigned long long)n, _indices[i]+1);
      }
      else {
         Rprintf("\"%llx\" [label=\"(%d, %g)\"];\n", (unsigned long long)n, n->vpindex+1, n->radius);
      }
   }


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
   //_tree.searchKNNKnownIndex(438, 8, 5);
   //_tree.print();
   IntegerMatrix im = _tree.hierarchicalClustering();
   delete rf;
   return im;
}


#endif /* VPTREEBINHIERARCHICAL_H_ */
