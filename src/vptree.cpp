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

// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

using namespace Rcpp;
using namespace std;

struct RFunction {
   RFunction(Function _f) : f(_f) {
      R_PreserveObject(f);
   }

   ~RFunction() {
      R_ReleaseObject(f);
   }

   Function f;
};

struct distClass
{
   RFunction* distance;
   bool isSimilarity;

   /*friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
   //ar & distance;
   }*/

   double operator()(const RObject& v1, const RObject& v2)
   {
      NumericVector res = distance->f(v1,v2);
      return isSimilarity ? 1.0-res[0] : res[0];
   }
};

template<typename T>
class VpTree
{
   public:
   static const char* ClassName;

   VpTree(int m=4, int minM=2, int max_leaf_size=25, int vantage_point_candidates=5, int test_point_count=15)
   : _m(m), _minM(minM), _max_leaf_size(max_leaf_size), _vantage_point_candidates(vantage_point_candidates),
      _test_point_count(test_point_count), _root(NULL)
   {
      if(max_leaf_size < vantage_point_candidates + test_point_count)
      {
         stop("Warning: max leaf size is too small");
      }
   }

   ~VpTree() {
      size_t n = _items.size();
      for(size_t i=0; i<n; ++i)
         R_ReleaseObject(_items[i]);
      delete _root;
   }

   vector<T> getCopyItems()
   {
      return vector<T>(_items.begin(), _items.end());
   }

   void create( const std::vector<T>& items) {
      delete _root;
      //_distance = distance;
      _items = items;
      size_t n = items.size();
      vector<int> indices(n);
      for (size_t i=0; i<n; ++i)
         indices[i] = i;
      //_root = buildFromPoints(indices, _minM, _minM, NULL, NULL);
      _root = buildFromPoints(indices, _m, _m, NULL, NULL);
   }

   void setItems(const std::vector<T>& items)
   {
      _items = items;
   }

   void setDistanceFunction(const distClass& distance_arg)
   {
      _distance = distance_arg;
   }

   void searchKNN( const T& target, int k, std::vector<T>* results,
   std::vector<double>* distances)
   {
      std::priority_queue<HeapItem> heap;

      _tau = std::numeric_limits<double>::max();
      search( _root, target, true, k, heap );

      results->clear(); distances->clear();

      while( !heap.empty() ) {
         results->push_back( _items[heap.top().index] );
         distances->push_back( heap.top().dist );
         heap.pop();
      }

      std::reverse( results->begin(), results->end() );
      std::reverse( distances->begin(), distances->end() );
   }

   void searchRadius( const T& target, double tau, std::vector<T>* results,
   std::vector<double>* distances)
   {
      std::priority_queue<HeapItem> heap;

      _tau = tau;
      search( _root, target, false, -1, heap );

      results->clear(); distances->clear();

      while( !heap.empty() ) {
         results->push_back( _items[heap.top().index] );
         distances->push_back( heap.top().dist );
         heap.pop();
      }

      std::reverse( results->begin(), results->end() );
      std::reverse( distances->begin(), distances->end() );
   }

   void insert(const T& target)
   {
      _items.push_back(target);
      insert_rec(_root, target);
   }

   private:
   distClass _distance;
   std::vector<T> _items;
   double _tau;

   // The following condition must held:
   // MAX_LEAF_SIZE >= VANTAGE_POINT_CANDIDATES + TEST_POINT_COUNT
   int _m; //max number of children
   int _minM; //min number of children
   size_t _max_leaf_size;
   int _vantage_point_candidates;
   int _test_point_count;


   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      //ar & _items;
      //ar & _tau;
      ar & _root;
      ar & _max_leaf_size;
      ar & _vantage_point_candidates;
      ar & _test_point_count;
      ar & _m;
      ar & _minM;
      //ar & distance;
   }

   struct Node
   {
      int vpindex;
      //Node* left;
      //Node* right;
      bool isLeaf;
      vector<double> radiuses;
      vector<int> points;
      vector<Node*> children;
      int childCount;
      Node* parent;

      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
         ar & vpindex;
         //for(int i=0; i<children.size(); i++)
         //	ar & children[i];
         ar & children;
         ar & parent;
         ar & points;
         ar & radiuses;
         ar & isLeaf;
         ar & childCount;
      }

      Node() :
      vpindex(0), isLeaf(false), childCount(0) {}

      ~Node() {
         for(int i=0;i<childCount;i++)
         {
            if(children[i] != NULL)
            delete children[i];
         }
      }
   };

   Node* _root;

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
      vector<T>* items;
      int index;
      distClass distance;

      //_items, indices[0], _distance
      DistanceComparator(vector<T>* items, int index, distClass distance )
         : items(items), index(index), distance(distance) {}
      bool operator()(int a, int b) {
         return distance( (*items)[index], (*items)[a] ) < distance( (*items)[index], (*items)[b] );
      }
   };



   void splitLeaf(Node* node, const T& target)
   {
      Node* parent = node->parent;
      if(parent == NULL) //strange case, when we have only one leaf in whole tree, very small number of items
      {
         //Rcout << "splitLeaf - tree grows in height, items size = " << _items.size() << endl;
         delete _root;
         size_t n = _items.size();
         vector<int> indices(n);
         for(size_t i=0; i<n; ++i)
            indices[i] = i;
         _root = buildFromPoints(indices, _minM, _minM, NULL);
         return;
      }
      //Rcout << "splitLeaf - przebudowuje" << endl;
      vector<int> allitems;
      allitems.push_back(_items.size() - 1);

      for(int i=0;i<parent->childCount;i++)
      {
         vector<int> items = parent->children[i]->points;
         allitems.insert(allitems.end(), items.begin(), items.end());
      }

      buildFromPoints(allitems, parent->childCount+1, _minM, parent->parent, parent);
   }

   void splitNonLeaf(Node* node, const T& target)
   {
      Node* parent = node->parent;
      while(parent != NULL && parent->childCount >= _m)
      parent = parent->parent;

      if(parent == NULL) // If root is full, then tree grows in height
      {
         //Rcout << "splitNonLeaf - tree grows in height, items size = " << _items.size() << endl;
         delete _root;
         size_t n = _items.size();
         vector<int> indices(n);
         for(size_t i=0; i<n; ++i)
            indices[i] = i;
         _root = buildFromPoints(indices, _minM, _minM, NULL, NULL);
         return;
      }
      //Rcout << "splitNonLeaf - przebudowuje" << endl;
      vector<int> allitems = getAllItemsFromNode(parent);
      allitems.push_back(_items.size() - 1);
      buildFromPoints(allitems, parent->childCount+1, _minM, parent->parent, parent);
   }

   vector<int> getAllItemsFromNode(Node* node)
   {
      if(node->isLeaf)
      return node->points;

      vector<int> ret;
      for(int i=0;i<node->childCount;i++)
      {
         vector<int> items = getAllItemsFromNode(node->children[i]);
         ret.insert(ret.end(), items.begin(), items.end());
      }
      return ret;
   }

   void insert_rec(Node* node, const T& target)
   {
      if ( node == NULL ) return;
      if(node->isLeaf)
      {
         if(node->points.size() < _max_leaf_size)
         {
            //Rcout << "dodaje normalnie, itemow w drzewie: " << _items.size() << endl;
            node->points.push_back(_items.size() - 1);
            //Rcout << _items[_items.size() - 1][0] << ", " << _items[_items.size() - 1][1] << endl;
            /*for(int i=0;i<node->points.size(); i++)
            Rcout << node->points[i] << " ";
            Rcout << endl;*/
         }
         else
         {
            if(node->parent->childCount < _m)
            {
               //Rcout << "splitLeaf" << endl;
               splitLeaf(node, target);
            }
            else
            {
               //Rcout << "splitNonLeaf" << endl;
               splitNonLeaf(node, target);
            }
         }
      }
      else
      {
         double dist = _distance(_items[node->vpindex], target);
         for(int i=0;i<node->childCount-1;i++)
         {
            if ( dist < node->radiuses[i] ) {
               insert_rec(node->children[i], target);
               return;
            }
         }
         insert_rec( node->children[node->childCount-1], target);
      }
   }

   Node* buildFromPoints(vector<int> indices, int childCountThis, int childCountChildren, Node* parent, Node* node=NULL )
   {
      if(node == NULL)
      node = new Node();
      else
      {
         for(int i=0;i<node->childCount;i++)
         delete node->children[i];
      }
      node->parent = parent;

      if ( indices.size() < _max_leaf_size  ) {
         node->isLeaf=true;
         node->points = indices;
         return node;
      }

      node->isLeaf=false;
      node->childCount = childCountThis;

      // choose an arbitrary point and move it to the start
      int vpi = chooseNewVantagePoint(indices); //(int)((double)rand() / RAND_MAX * (upper - lower - 1) ) + lower;
      std::swap( indices[0], indices[vpi] );
      std::sort(indices.begin() + 1, indices.end(), DistanceComparator( &_items, indices[0], _distance ));

      node->radiuses = vector<double>(childCountThis-1);
      node->children = vector<Node*>(childCountThis);
      for(int i=0;i<childCountThis;i++)
      {
         int itemsPerNode = indices.size() / childCountThis;
         if(i < childCountThis - 1)
         {
            node->children[i] = buildFromPoints(
               vector<int>(indices.begin() + i*itemsPerNode, indices.begin() + (i+1)*itemsPerNode),
               childCountChildren, childCountChildren,
               node);
               node->radiuses[i] = (_distance( _items[indices[0]], _items[indices[0 + (i+1)*itemsPerNode - 1] ])
               + _distance( _items[indices[0]], _items[indices[0+ (i+1)*itemsPerNode] ])) / 2.0;
         }
         else
         {
            node->children[i] = buildFromPoints(
               vector<int>(indices.begin() + (childCountThis - 1)*itemsPerNode, indices.end()),
               childCountChildren, childCountChildren,
               node);
         }
      }
      node->vpindex = indices[0];
      return node;
   }

   int chooseNewVantagePoint(vector<int> indices)
   {
      vector<int> candidates(0);
      vector<int> testPoints(0);

      for (int i = 0; i < 0 + _vantage_point_candidates; ++i) {
         int basePointIndex = i + (int) (rand() % (indices.size() - i));
         candidates.push_back(indices[basePointIndex]);
         std::swap(indices[i], indices[basePointIndex] );
      }
      //Rcout << "wybralem kandydatow" << endl;
      for (int i = 0 + _vantage_point_candidates; i < 0 + _vantage_point_candidates + _test_point_count; ++i) {
         int testPointIndex = i + (int) (rand() % (indices.size() - i));
         testPoints.push_back(indices[testPointIndex]);
         std::swap(indices[i], indices[testPointIndex] );
      }
      //Rcout << "wybralem testerow" << endl;
      /*if(candidates.size() != VANTAGE_POINT_CANDIDATES)
      Rcout << "candidates.size() " << candidates.size() << endl;
      if(testPoints.size() != TEST_POINT_COUNT)
      Rcout << "testPoints.size() " << testPoints.size() << endl;*/
      double bestBasePointSigma = 0;
      int bestIndex = 0;
      for(size_t i=0;i<candidates.size();i++) {
         vector<double> distances(_test_point_count);
         for (int j = 0; j < _test_point_count; ++j) {
            distances[j] = _distance( _items[candidates[i]], _items[testPoints[j]] );
         }
         //Rcout << "przed liczeniem sigmy" << endl;
         double sigma = sigmaSquare(distances);
         //Rcout << "policzylem sigme" << endl;
         if (sigma > bestBasePointSigma) {
            bestBasePointSigma = sigma;
            bestIndex = i;
         }
      }

      return bestIndex;
   }

   static double sigmaSquare(const vector<double>& values) {
      double sum = std::accumulate(values.begin(),values.end(),0);

      double avg = sum / values.size();
      double sigmaSq = 0;

      size_t n = values.size();
      for (size_t i=0;i<n;++i) {
         double value = values[i];
         double dev = value - avg;
         sigmaSq += dev * dev;
      }

      return sigmaSq;
   }

   void search( Node* node, const T& target, bool isKNN, size_t k,
   std::priority_queue<HeapItem>& heap )
   {
      if ( node == NULL ) return;

      double dist = _distance( _items[node->vpindex], target );
      //printf("dist=%g tau=%gn", dist, _tau );
      if(node->isLeaf)
      {
         for(size_t i=0;i<node->points.size();++i)
         {
            double dist2 = _distance( _items[node->points[i]], target );
            if ( (dist2 < _tau && isKNN) || (dist2 <= _tau && !isKNN) ) {
               if ( heap.size() == k && isKNN) heap.pop();
               heap.push( HeapItem(node->points[i], dist2) );
               if ( heap.size() == k && isKNN)
               {
                  _tau = heap.top().dist;
                  //Rcout << "current tau=" << _tau << endl;
               }
            }
         }
      }
      else
      {
         vector<bool> visited(node->childCount, false);

         for(int i=0;i<node->childCount-1;i++)
         {
            //if ( dist < node->radiuses[i] ) {
            if ( dist - _tau <= node->radiuses[i] && !visited[i]) {
               search( node->children[i], target, isKNN, k, heap );
               visited[i]=true;
            }
            if ( dist + _tau >= node->radiuses[i] && !visited[i+1]) {
               search( node->children[i+1], target, isKNN, k, heap );
               visited[i+1]=true;
            }
            //}
         }
         /*if ( dist + _tau >= node->radiuses[node->childCount-2] ) {
         search( node->children[node->childCount-1], target, k, heap );
         }*/
      }
   }
};

template<typename T>
const char* VpTree<T>::ClassName = "VpTree";

template<typename T>
void checkIsVpTreeClass(XPtr< VpTree<T> >& _queue)
{
   if(strcmp(_queue.attr("class"), VpTree<T>::ClassName))
   stop("not a VpTree object");
}

// [[Rcpp::export]]
SEXP vptree_create(Function distance, bool isSimilarity = false, int m=2, int minm=4, int maxLeafPointsCount = 25, int vantageCandidatesCount = 5, int testPointsCount = 15) { //https://code.google.com/p/vptree/source/browse/src/vptree/VpTreeNode.java
VpTree<RObject>* vec = new VpTree<RObject>(m, minm, maxLeafPointsCount, vantageCandidatesCount, testPointsCount);
XPtr< VpTree<RObject> > retval =  XPtr< VpTree<RObject> >(vec, true);
retval.attr("class") = VpTree<RObject>::ClassName;
distClass distci;
distci.isSimilarity = isSimilarity;
distci.distance = new RFunction(distance);
(*retval).setDistanceFunction(distci);

return retval;
}

// [[Rcpp::export]]
void vptree_insert(SEXP tree, RObject obj) {
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (tree);
   checkIsVpTreeClass(_tree);
   R_PreserveObject(obj);
   (*_tree).insert(obj);
}

// [[Rcpp::export]]
void vptree_set_distancefunction(SEXP tree, Function distance, bool isSimilarity = false) {
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (tree);
   checkIsVpTreeClass(_tree);
   distClass distci;
   distci.isSimilarity = isSimilarity;
   distci.distance = new RFunction(distance);
   (*_tree).setDistanceFunction(distci);
}

vector<RObject> createStdVectorOfRobjects(List listobj)
{
   int n = listobj.size();
   vector<RObject> points(n);
   for(int i=0;i<n;i++)
   {
      RObject obj = listobj[i];
      R_PreserveObject(obj);
      points[i] = obj;
   }
   return points;
}

// [[Rcpp::export]]
void vptree_build(SEXP tree, List listobj) {
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (tree);
   checkIsVpTreeClass(_tree);
   vector<RObject> points = createStdVectorOfRobjects(listobj);
   (*_tree).create(points);
}

// [[Rcpp::export]]
List vptree_searchKNN(SEXP tree, RObject p, int k)
{
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (tree);
   checkIsVpTreeClass(_tree);
   std::vector<RObject> results;
   std::vector<double> distances;
   (*_tree).searchKNN( p, k, &results, &distances);

   List resultList(2);
   List elements(results.size());
   List distancesRcpp(results.size());
   int index =0;
   for (std::vector<RObject>::iterator it = results.begin(); it!=results.end(); ++it)
   {
      elements[index] = *it;
      distancesRcpp[index] = distances[index];
      index++;
   }
   resultList[0] = elements;
   resultList[1] = distancesRcpp;
   return resultList;
}

// [[Rcpp::export]]
List vptree_searchRadius(SEXP tree, RObject p, double tau)
{
   XPtr< VpTree<RObject> > _tree = Rcpp::as< XPtr< VpTree<RObject> > > (tree);
   checkIsVpTreeClass(_tree);
   std::vector<RObject> results;
   std::vector<double> distances;
   (*_tree).searchRadius( p, tau, &results, &distances);

   List resultList(2);
   List elements(results.size());
   List distancesRcpp(results.size());
   int index =0;
   for (std::vector<RObject>::iterator it = results.begin(); it!=results.end(); ++it)
   {
      elements[index] = *it;
      distancesRcpp[index] = distances[index];
      index++;
   }
   resultList[0] = elements;
   resultList[1] = distancesRcpp;
   return resultList;
}



