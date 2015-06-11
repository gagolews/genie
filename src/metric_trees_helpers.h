#ifndef METRIC_TREES_HELPERS_
#define METRIC_TREES_HELPERS_
#include <Rcpp.h>
#define USE_RINTERNALS
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <unordered_map>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "boost/serialization/unordered_map.hpp"

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

class Point
{
   public:
   int i;
   int j;
   Point():i(0),j(0){}
   Point(int i, int j):i(i),j(j){}

   friend class boost::serialization::access;
   template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
   {
      ar & i;
      ar & j;
   }

   static Point createValidPoint(int i, int j)
   {
      if(j < i)
      {
         swap(i,j);
      }
      return Point(i,j);
   }

   bool operator==(const Point &other) const
   {
      return (i == other.i && j == other.j);
   }
};

ostream& operator<< (ostream& os, const Point& obj);
istream& operator>> (istream& is, Point& obj);


namespace std {

   template <>
      struct hash<Point>
   {
      std::size_t operator()(const Point& k) const
      {
         return (51 + hash<int>()(k.i)) * 51 + hash<int>()(k.j);
      }
   };

}

struct distClass
{
   RFunction* distance;
   bool isSimilarity;
   vector<RObject> *items;

   unordered_map<Point, double> hashmap;

   friend class boost::serialization::access;
   template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
   {
      ar & isSimilarity;
      ar & hashmap;
   }

   double operator()(const RObject& v1, const RObject& v2) const
   {
      NumericVector res = distance->f(v1,v2);
      return isSimilarity ? 1.0-res[0] : res[0];
   }

   double operator()(int v1, int v2)
   {
      if(v1==v2) return 0;
      Point p = Point::createValidPoint(v1,v2);
      std::unordered_map<Point,double>::const_iterator got = hashmap.find(p);
      if ( got == hashmap.end() )
      {
         NumericVector res = distance->f((*items)[v1],(*items)[v2]);
         double d = isSimilarity ? 1.0-res[0] : res[0];
         hashmap.emplace(p, d);
         return d;
      }
      else
      {
         return got->second;
      }
   }
};

vector<RObject> createStdVectorOfRobjects(List listobj);

#endif /* METRIC_TREES_HELPERS_ */
