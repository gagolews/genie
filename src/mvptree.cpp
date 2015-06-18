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

// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "boost/serialization/unordered_map.hpp"
#include "metric_trees_helpers.h"
#include "mvptree/mvptree.h"

#define MVP_BRANCHFACTOR 3
#define MVP_PATHLENGTH   5
#define MVP_LEAFCAP     25

using namespace Rcpp;
using namespace std;


class mvptree
{
   public:
   static const char* ClassName;
   MVPTree *_tree;
   distClass _distance;

   mvptree(MVPTree *tree, distClass distance) : _tree(tree), _distance(distance){}
};

mvptree* globaltree = NULL;

float distanceWrapper(MVPDP *pointA, MVPDP *pointB)
{
   if(globaltree == NULL) stop("global tree is null.");
   return globaltree->_distance(pointA, pointB);
}



const char* mvptree::ClassName = "mvptree";

void checkIsMvpTreeClass(XPtr<mvptree>& _tree)
{
   if(strcmp(_tree.attr("class"), mvptree::ClassName))
      stop("not an mvptree object");
}


//' @rdname mvptree
//' @details
//' \code{mvptree_create} creates an empty m-tree instance. Please use
//' \code{mvptree_build} for populating the tree.
//'
//' @return
//' \code{mvptree_create} returns a new, empty m-tree.
//' @param distance an R function, which calculates metric in the space
//' @param isSimilarity a logical value, does a distance function calculate
//' distance or similarity?
//' @param m single integer, maximum number of children can be in a node
//' in the tree
//' @param minm single integer, minimum number of children can be in a node
//' in the tree
//' @param maxLeafPointsCount single integer, maximum number of
//' points in a single leaf
//' @param vantageCandidatesCount single integer, number of candidates
//' considered as a vantage point in a single node. Too big number can
//' cause performance problem in building the tree, but too small can cause
//' performance problem when searching.
//' @param testPointsCount single integer, how many points are considered with
//' a single candidate for vantage point in a node to assess a variance of
//' distances
// [[Rcpp::export]]
SEXP mvptree_create(Function distance, bool isSimilarity = false, int mvp_branchfactor=2, int mvp_pathlength=5, int mvp_leafcap=25) {
   distClass distci;
   distci.isSimilarity = isSimilarity;
   distci.distance = new RFunction(distance);

   MVPTree *tree = mvptree_alloc(NULL, distanceWrapper, mvp_branchfactor, mvp_pathlength, mvp_leafcap);

   mvptree* vec = new mvptree(tree, distci);
   globaltree = vec;
   XPtr< mvptree > retval =  XPtr< mvptree >(vec, true);
   retval.attr("class") = mvptree::ClassName;
   return retval;
}

MVPDP* generate_point(RObject obj){
   /* generate one datapoint with BYTEARRAY type data of length dp_length */
   //stop("tworze zmienna statyczna...");
   static unsigned long long uid = 0;
   //stop("tworze napis...");
   char scratch[32];
   //stop("tworze wskaznik...");
   MVPDP *newpnt = NULL;


   //stop("alokuje pamiec...");
   if(sizeof(RObject*) == 4)
   {
      newpnt = dp_alloc(MVP_UINT32ARRAY);
      newpnt->data = malloc(sizeof(RObject));
      RObject* row = (RObject*)newpnt->data;
      row[0] = obj;
      R_PreserveObject(row[0]);
   }
   else if(sizeof(RObject*) == 8)
   {
      newpnt = dp_alloc(MVP_UINT64ARRAY);
      newpnt->data = malloc(sizeof(RObject));
      RObject* row = (RObject*)newpnt->data;
      row[0] = obj;
      R_PreserveObject(row[0]);
   }
   else
      stop("sizeof(RObject*) is not 4 nor 8");
   newpnt->datalen = 1;
   //stop("kopiuje id...");
   snprintf(scratch, 32, "point%llu", ++uid);
   newpnt->id = strdup(scratch);

   return newpnt;
}

//' @rdname mvptree
//' @details
//' \code{mvptree_insert} inserts a new point into a m-tree.
//' Calling this function n times has worse performance in comparison
//' to calling one \code{mvptree_build} with n points in the beginning.
//'
//' @return
//' \code{mvptree_insert} does not return anything interesting.
//' @param obj a point to insert
// [[Rcpp::export]]
void mvptree_insert(SEXP tree, RObject obj) {
   XPtr< mvptree > _tree = Rcpp::as< XPtr< mvptree > > (tree);
   checkIsMvpTreeClass(_tree);
   R_PreserveObject(obj);
   MVPDP* obj2 = generate_point(obj);
   MVPDP **pointlist = (MVPDP**)malloc(1*sizeof(MVPDP*));
   pointlist[0] = obj2;
   MVPError err = mvptree_add(_tree->_tree, pointlist, 1);
   if(err != MVP_SUCCESS)
      stop(std::to_string(err));
}

MVPDP** createMVPDPFromVector(const vector<RObject>& vec)
{
   //stop("tworze liste punktow...");
   MVPDP **pointlist = (MVPDP**)malloc(vec.size()*sizeof(MVPDP*));
   unsigned int i;
   for (i = 0; i < vec.size(); i++)
   {
      //stop("tworze punkt...");
	   pointlist[i] = generate_point(vec[i]);
      //stop("stworzylem punkt...");
   }
   return pointlist;
}

//' @rdname mvptree
//' @details
//' \code{mvptree_build} populates an empty tree with set of points. If called
//' for non-empty tree, whole structure of previous tree is deleted. This is
//' recommended way of building a tree (instead of inserts by
//' \code{mvptree_insert}).
//'
//' @return
//' \code{mvptree_build} does not return anything interesting.
//' @param listobj an R list, where every element of list is assumed as one
//' point in a  metric space
// [[Rcpp::export]]
void mvptree_build(SEXP tree, List listobj) {
   XPtr< mvptree > _tree = Rcpp::as< XPtr< mvptree > > (tree);
   checkIsMvpTreeClass(_tree);
   vector<RObject> points = createStdVectorOfRobjects(listobj);
   //stop("Tworze ten dziwny wektor...");
   //Rcout<<"Tworze ten dziwny wektor..."<<endl;
   MVPDP  **pointlist = createMVPDPFromVector(points);
   //stop("Dodaje do drzewa...");
   //Rcout<<"Dodaje do drzewa " << endl;
   globaltree = _tree;
   MVPError err = mvptree_add(_tree->_tree, pointlist, points.size());
   if(err != MVP_SUCCESS)
      stop(std::to_string(err));
}

//' @rdname mvptree
//' @details
//' \code{mvptree_searchKNN} finds k nearest neighbours for a given object
//' in the tree.
//' Please use this function for objects, which are not in a set of points
//' inserted to the tree or which were used for building the tree. For these
//' cases please use \code{mvptree_searchKNNKnown} or
//' \code{mvptree_searchKNNKnownIndex} which can be faster.
//'
//' @return
//' \code{mvptree_searchKNN} returns a list. First element of the list
//' is a list with found elements. A second element of the list is a list
//' with distances these objects from a given object.
//' @param p an R object for which neighbours are found
//' @param k a single integer, number of neighbours to find
//' @param findItself boolean value, should results contain an given object?
// [[Rcpp::export]]
List mvptree_searchKNN(SEXP tree, RObject p, int k, bool findItself = true)
{
   XPtr< mvptree > _tree = Rcpp::as< XPtr< mvptree > > (tree);
   checkIsMvpTreeClass(_tree);

   MVPDP* p1 = generate_point(p);

   unsigned long long nbcalcs = 0;
   unsigned int nbresults;
   MVPError err;
   globaltree = _tree;
   std::priority_queue<HeapItemMVPTree> results = mvptree_retrieve(_tree->_tree, p1 , k, std::numeric_limits<double>::max(), &nbresults, &err, findItself);

   size_t resultsSize = results.size();

   List resultList(2);
   List elements(resultsSize);
   List distancesRcpp(resultsSize);
   int index = 0;
   while( !results.empty() ) {
         elements[resultsSize-index-1] = *results.top().obj;
         distancesRcpp[resultsSize-index-1] = results.top().dist;
         index++;
         results.pop();
      }
   resultList[0] = elements;
   resultList[1] = distancesRcpp;
   return resultList;
}

//' @rdname mvptree
//' @details
//' \code{mvptree_searchRadius} finds neighbours for a given object
//' within a given radius in the tree.
//' Please use this function for objects, which are not in a set of points
//' inserted to the tree or which were used for building the tree. For these
//' cases please use \code{mvptree_searchRadiusKnown} or
//' \code{mvptree_searchRadiusKnownIndex} which can be faster.
//'
//' @return
//' \code{mvptree_searchRadius} returns a list. First element of the list
//' is a list with found elements. A second element of the list is a list
//' with distances these objects from a given object.
//' @param p an R object for which neighbours are found
//' @param tau a float value, a radius
//' @param findItself boolean value, should results contain an given object?
// [[Rcpp::export]]
List mvptree_searchRadius(SEXP tree, RObject p, double tau, bool findItself = true)
{XPtr< mvptree > _tree = Rcpp::as< XPtr< mvptree > > (tree);
   checkIsMvpTreeClass(_tree);

   MVPDP* p1 = generate_point(p);

   unsigned long long nbcalcs = 0;
   unsigned int nbresults;
   MVPError err;
   globaltree = _tree;
   std::priority_queue<HeapItemMVPTree> results = mvptree_retrieve(_tree->_tree, p1 , std::numeric_limits<unsigned int>::max(), tau, &nbresults, &err, findItself);

   size_t resultsSize = results.size();

   List resultList(2);
   List elements(resultsSize);
   List distancesRcpp(resultsSize);
   int index = 0;
   while( !results.empty() ) {
         elements[resultsSize-index-1] = *results.top().obj;
         distancesRcpp[resultsSize-index-1] = results.top().dist;
         index++;
         results.pop();
      }
   resultList[0] = elements;
   resultList[1] = distancesRcpp;
   return resultList;
}
