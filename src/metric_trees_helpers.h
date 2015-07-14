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
#include "mvptree/mvptree.h"

using namespace Rcpp;
using namespace std;

//Sys.setenv(PKG_CPPFLAGS="-DDEBUG")

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

class Value
{
public:
   double dist;
#ifdef DEBUG
   int counter;
#endif

   Value() : dist(0)
#ifdef DEBUG
   , counter(1)
#endif
   {}

   Value(double dist) : dist(dist)
#ifdef DEBUG
   , counter(1)
#endif
   {}

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & dist;
#ifdef DEBUG
      ar & counter;
#endif
   }

#ifdef DEBUG
   void incrementCounter()
   {
      counter++;
   }
#endif
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

   unordered_map<Point, Value> hashmap;

   #ifdef DEBUG
   int metricCalculated;
   int hashmapHit;
   #endif

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

   float operator()(MVPDP *pointA, MVPDP *pointB) const
   {
      RObject *o1 = (RObject*)pointA->data;
      RObject *o2 = (RObject*)pointB->data;

      return (*this)(*o1,*o2);
   }

   double operator()(int v1, int v2)
   {
#ifdef DEBUG
      metricCalculated++;
#endif
      if(v1==v2) return 0;
      Point p = Point::createValidPoint(v1,v2);
      std::unordered_map<Point,Value>::iterator got = hashmap.find(p);
      if ( got == hashmap.end() )
      {
         NumericVector res = distance->f((*items)[v1],(*items)[v2]);
         double d = isSimilarity ? 1.0-res[0] : res[0];
         hashmap.emplace(p, Value(d));
         return d;
      }
      else
      {
#ifdef DEBUG
         got->second.incrementCounter();
         hashmapHit++;
#endif
         //Rcout << "trafilem"<<endl;
         return got->second.dist;
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

vector<RObject> createStdVectorOfRobjects(List listobj);

#endif /* METRIC_TREES_HELPERS_ */
