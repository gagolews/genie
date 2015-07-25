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


#include "hclust2_distance.h"
using namespace DataStructures;


Distance::Distance(size_t n) :
#ifdef HASHMAP_ENABLED
      hashmap(std::vector< std::unordered_map<size_t, double> >(n)),
#endif
#ifdef DIST_COUNTERS
#ifdef HASHMAP_ENABLED
      hashmapHit(0),
      hashmapMiss(0),
#endif
      distCallCount(0),
#endif
      n(n)
{
#ifdef HASHMAP_ENABLED
   Rcpp::Rcout << "Warning: HASHMAP_ENABLED is defined in hclust2_distance.h\n";
#endif
#ifdef DIST_COUNTERS
   Rcpp::Rcout << "Warning: DIST_COUNTERS is defined in hclust2_distance.h\n";
#endif
}


Distance::~Distance()
{
// #if VERBOSE > 5
//    Rprintf("[%010.3f] destroying distance object (base)\n", clock()/(float)CLOCKS_PER_SEC);
// #endif
#if VERBOSE > 3 && defined(HASHMAP_ENABLED) && defined(DIST_COUNTERS)
   Rprintf("             distance function hashmap #hits: %.0f, #miss: %.0f, est.mem.used: ~%.1fMB (vs %.1fMB)\n",
      (double)hashmapHit, (double)hashmapMiss,
      8.0f*hashmapMiss/1000.0f/1000.0f,
      8.0f*(n-1)*(n-1)*0.5f/1000.0f/1000.0f);
#endif
#if VERBOSE > 0 && defined(DIST_COUNTERS)
   Rprintf("             distance function total calls: %.0f (i.e., %.2f%% of %.0f)\n",
      (double)distCallCount,
      (double)distCallCount*100.0/(0.5*n*(n-1.0)),
      0.5*n*(n-1.0)
   );
#endif
}


#ifdef HASHMAP_ENABLED
double Distance::operator()(size_t v1, size_t v2)
{
   if (v1 == v2) return 0.0;
   if (v1 > v2) std::swap(v1, v2);

#ifdef DIST_COUNTERS
   ++distCallCount;
#endif

   // std::unordered_map<SortedPoint,double>::iterator got = hashmap.find(p);
   auto got = hashmap[v1].find(v2);
   if ( got == hashmap[v1].end() )
   {
#ifdef DIST_COUNTERS
      ++hashmapMiss;
#endif
      double d = compute(v1, v2);
      hashmap[v1].emplace(v2, d);
      return d;
   }
   else
   {
#ifdef DIST_COUNTERS
      ++hashmapHit;
#endif
      return got->second;
   }
}
#endif




Distance* Distance::createDistance(Rcpp::RObject distance, Rcpp::RObject objects)
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
   else if (Rf_isNumeric(distance) && Rf_isObject(distance) && !strcmp(distance.attr("class"), "dist") && Rf_isNull(objects))
   {
      return (DataStructures::Distance*)
            new DataStructures::DistObjectDistance(
               (Rcpp::NumericVector)distance
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
      d += (items[v1*m+i]-items[v2*m+i])*(items[v1*m+i]-items[v2*m+i]);
   return sqrt(d);

// this is not faster:
//    double d = sqobs[v1]+sqobs[v2]; // already multiplied by 0.5
//    // sum((x-y)^2) == 2*(sum(x^2)/2 + sum(y^2)/2 - sum(x*y))
//    const double* items1_ptr = items+v1*m;
//    const double* items2_ptr = items+v2*m;
//    for (size_t i=0; i<m; ++i)
//       d -= (*(items1_ptr++))*(*(items2_ptr++));
//    return sqrt(2.0*d);
}


double ManhattanDistance::compute(size_t v1, size_t v2) const
{
   if (v1 == v2) return 0.0;
   double d = 0.0;
   for (size_t i=0; i<m; ++i)
      d += std::abs(items[v1*m+i]-items[v2*m+i]);
   return d;
}


double MaximumDistance::compute(size_t v1, size_t v2) const
{
   if (v1 == v2) return 0.0;
   double d = 0.0;
   for (size_t i=0; i<m; ++i) {
      double d2 = std::abs(items[v1*m+i]-items[v2*m+i]);
      if (d2 > d) d = d2;
   }
   return d;
}


double GenericRDistance::compute(size_t v1, size_t v2) const
{
   if (v1 == v2) return 0.0;
   return ((Rcpp::NumericVector)distfun(items[v1], items[v2]))[0];
}


double DistObjectDistance::compute(size_t v1, size_t v2) const
{
   if (v1 == v2) return 0.0;

   size_t i;
   if (v1 < v2)
      i = n*v1-((v1+1)*(v1))/2+v2-v1-1;
   else
      i = n*v2-((v2+1)*(v2))/2+v1-v2-1;

   return items[i];
}

