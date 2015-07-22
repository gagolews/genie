#include "hclust2_distance.h"
using namespace DataStructures;


Distance::Distance(size_t n) :
#ifndef HASHMAP_DISABLE
      hashmap(std::vector< std::unordered_map<size_t, double> >(n)),
#ifdef HASHMAP_COUNTERS
      hashmapHit(0),
      hashmapMiss(0),
#endif
#endif
      n(n)
{

}


Distance::~Distance()
{
// #if VERBOSE > 5
//    Rprintf("[%010.3f] destroying distance object (base)\n", clock()/(float)CLOCKS_PER_SEC);
// #endif
#if VERBOSE > 3 && !defined(HASHMAP_DISABLE) && defined(HASHMAP_COUNTERS)
   Rprintf("Distance Hashmap: #hits=%d, #miss=%d, est.mem.used>=%.1fMB (vs %.1fMB)\n",
      hashmapHit, hashmapMiss, 8.0f*hashmapMiss/1000.0f/1000.0f,
      8.0f*(n-1)*(n-1)*0.5f/1000.0f/1000.0f);
#endif
}


#ifndef HASHMAP_DISABLE
double Distance::operator()(size_t v1, size_t v2)
{
   if (v1 == v2) return 0.0;
   if (v1 > v2) std::swap(v1, v2);

   // std::unordered_map<SortedPoint,double>::iterator got = hashmap.find(p);
   auto got = hashmap[v1].find(v2);
   if ( got == hashmap[v1].end() )
   {
#ifdef HASHMAP_COUNTERS
      ++hashmapMiss;
#endif
      double d = compute(v1, v2);
      hashmap[v1].emplace(v2, d);
      return d;
   }
   else
   {
#ifdef HASHMAP_COUNTERS
      ++hashmapHit;
#endif
      return got->second;
   }
}
#endif




Distance* Distance::createDistance(Rcpp::RObject objects, Rcpp::RObject distance)
{
   if (Rf_isVectorList(objects) && Rf_isFunction(distance))
   {
      Rcpp::Function distance2(distance);
      Rcpp::List objects2(objects);
      return (DataStructures::Distance*)
         new DataStructures::GenericRDistance(
            distance2,
            std::vector<Rcpp::RObject>(objects2.begin(), objects2.end())
         );
   }
   else if (Rf_isMatrix(objects) && (Rf_isNull(distance) || Rf_isString(distance)))
   {
      Rcpp::NumericMatrix objects2(objects);
      Rcpp::CharacterVector distance2 =
         ((Rf_isNull(distance))?Rcpp::CharacterVector("euclidean"):Rcpp::CharacterVector(distance));

      const char* distance3 = CHAR(STRING_ELT((SEXP)distance2, 0));
      if (!strcmp(distance3, "euclidean")) {
         return (DataStructures::Distance*)
            new DataStructures::EuclideanDistance(
               objects2
            );
      }
      else if (!strcmp(distance3, "manhattan")) {
         return (DataStructures::Distance*)
            new DataStructures::ManhattanDistance(
               objects2
            );
      }
      else if (!strcmp(distance3, "maximum")) {
         return (DataStructures::Distance*)
            new DataStructures::MaximumDistance(
               objects2
            );
      }
      else {
         Rcpp::stop("`distance` should be one of: \"euclidean\" (default), \"manhattan\", \"maximum\"");
      }
   }
   else {
      Rcpp::stop("incorrect input data");
   }
   return NULL;
}



double EuclideanDistance::compute(size_t v1, size_t v2) const
{
   if (v1 == v2) return 0.0;
   double d = 0.0;
   for (size_t i=0; i<m; ++i)
      d += (items[v1+i*n]-items[v2+i*n])*(items[v1+i*n]-items[v2+i*n]);
   return sqrt(d);
}


double ManhattanDistance::compute(size_t v1, size_t v2) const
{
   if (v1 == v2) return 0.0;
   double d = 0.0;
   for (size_t i=0; i<m; ++i)
      d += std::abs(items[v1+i*n]-items[v2+i*n]);
   return d;
}


double MaximumDistance::compute(size_t v1, size_t v2) const
{
   if (v1 == v2) return 0.0;
   double d = 0.0;
   for (size_t i=0; i<m; ++i) {
      double d2 = std::abs(items[v1+i*n]-items[v2+i*n]);
      if (d2 > d) d = d2;
   }
   return d;
}


double GenericRDistance::compute(size_t v1, size_t v2) const
{
   if (v1 == v2) return 0;
   return ((Rcpp::NumericVector)distfun(items[v1],items[v2]))[0];
}
