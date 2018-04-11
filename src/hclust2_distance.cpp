/* ************************************************************************* *
 *   This file is part of the `genie` package for R.                         *
 *                                                                           *
 *   Copyright 2015-2018 Marek Gagolewski, Maciej Bartoszuk, Anna Cena       *
 *                                                                           *
 *   'genie' is free software: you can redistribute it and/or                *
 *   modify it under the terms of the GNU General Public License             *
 *   as published by the Free Software Foundation, either version 3          *
 *   of the License, or (at your option) any later version.                  *
 *                                                                           *
 *   'genie' is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with 'genie'. If not, see <http://www.gnu.org/licenses/>.         *
 * ************************************************************************* */



#include <algorithm>
#include "hclust2_distance.h"
using namespace grup;


// ------------------------------------------------------------------------
#ifdef MEASURE_MEM_USE
/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 *
 *
 *          http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c
 */

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif





/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t getPeakRSS( )
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
        return (size_t)0L;      /* Can't open? */
    if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
    {
        close( fd );
        return (size_t)0L;      /* Can't read? */
    }
    close( fd );
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;          /* Unsupported. */
#endif
}


/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t getCurrentRSS( )
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
        return (size_t)0L;      /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE* fp = NULL;
    if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
        return (size_t)0L;      /* Can't open? */
    if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    {
        fclose( fp );
        return (size_t)0L;      /* Can't read? */
    }
    fclose( fp );
    return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;          /* Unsupported. */
#endif
}
#endif
// ------------------------------------------------------------------------



void DistanceStats::print() const
{
#if VERBOSE > 0
// #if defined(HASHMAP_ENABLED) && defined(GENERATE_STATS)
//    Rprintf("             distance function hashmap #hits: %.0f, #miss: %.0f, est.mem.used: ~%.1fMB (vs %.1fMB)\n",
//       (double)hashmapHit, (double)hashmapMiss,
//       8.0f*hashmapMiss/1000.0f/1000.0f,
//       8.0f*distCallTheoretical/1000.0f/1000.0f);
// #endif
#if defined(GENERATE_STATS)
   Rprintf("             distance function total calls: %.0f (i.e., %.2f%% of %.0f)\n",
      (double)distCallCount,
      (double)distCallCount*100.0/(double)distCallTheoretical,
      (double)distCallTheoretical
   );
#if defined(MEASURE_MEM_USE)
   Rprintf("             currentRSS=%.0f MB, peakRSS=%.0f MB\n",
      (double)getCurrentRSS()/1000.0/1000.0,
      (double)getPeakRSS()/1000.0/1000.0
   );
#endif
#endif
#endif
}


Distance::Distance(size_t n) :
// #ifdef HASHMAP_ENABLED
//    hashmap(std::vector< std::unordered_map<size_t, double> >(n)),
// #endif
   stats(DistanceStats(n)),
   n(n)
{
// #ifdef HASHMAP_ENABLED
//    MESSAGE_1("Warning: HASHMAP_ENABLED is defined in hclust2_distance.h\n");
// #endif
#ifdef GENERATE_STATS
   MESSAGE_1("Warning: GENERATE_STATS is defined in hclust2_distance.h\n");
#endif
}


Distance::~Distance()
{
// #if VERBOSE > 5
//    Rprintf("[%010.3f] destroying distance object (base)\n", clock()/(float)CLOCKS_PER_SEC);
// #endif
}


// #ifdef HASHMAP_ENABLED
// double Distance::operator()(size_t v1, size_t v2)
// {
//    if (v1 == v2) return 0.0;
//    if (v1 > v2) std::swap(v1, v2);
//
// #ifdef GENERATE_STATS
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
//       ++stats.distCallCount;
// #endif
//
//    // this is thread unsafe, but we use it only for testing:
//    auto got = hashmap[v1].find(v2);
//    if ( got == hashmap[v1].end() )
//    {
// #ifdef GENERATE_STATS
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
//       ++stats.hashmapMiss;
// #endif
//       double d = compute(v1, v2);
//       hashmap[v1].emplace(v2, d);
//       return d;
//    }
//    else
//    {
// #ifdef GENERATE_STATS
// #ifdef _OPENMP
// #pragma omp atomic
// #endif
//       ++stats.hashmapHit;
// #endif
//       return got->second;
//    }
// }
// #endif


Distance* Distance::createDistance(Rcpp::RObject distance, Rcpp::RObject objects, Rcpp::RObject control)
{
   if (Rf_isVectorList(objects) && Rf_isFunction(distance))
   {
      Rcpp::Function distance2(distance);
      Rcpp::List objects2(objects);
      return (grup::Distance*)
         new grup::GenericRDistance(
            distance2,
            objects2
         );
   }
   else if (Rf_isNumeric(distance) && Rf_isObject(distance) && !strcmp(distance.attr("class"), "dist") && Rf_isNull(objects))
   {
      return (grup::Distance*)
            new grup::DistObjectDistance(
               (Rcpp::NumericVector)distance
            );
   }
   else if (Rf_isVectorList(objects) && (Rf_isNull(distance) || Rf_isString(distance)))
   {
      Rcpp::List objects2(objects);
      Rcpp::CharacterVector distance2 =
         ((Rf_isNull(distance))?Rcpp::CharacterVector("levenshtein"):Rcpp::CharacterVector(distance));

      const char* distance3 = CHAR(STRING_ELT((SEXP)distance2, 0));
      if (!strcmp(distance3, "levenshtein")) {
         return (grup::Distance*)new grup::LevenshteinDistanceInt(objects2);
      }
      else if (!strcmp(distance3, "dinu")) {
         return (grup::Distance*)new grup::DinuDistanceInt(objects2);
      }
      else if (!strcmp(distance3, "hamming")) {
         return (grup::Distance*)new grup::HammingDistanceInt(objects2);
      }
      else if (!strcmp(distance3, "euclinf")) {
         Rcpp::List control2(control);
         double p, r;

         if (control2.containsElementNamed("p"))
            p = (size_t)Rcpp::as<Rcpp::NumericVector>(control2["p"])[0];
         else
            Rcpp::stop("In euclinf p should be given.");

         if (control2.containsElementNamed("r"))
            r = (size_t)Rcpp::as<Rcpp::NumericVector>(control2["r"])[0];
         else
            Rcpp::stop("In euclinf r should be given.");

         return (grup::Distance*)new grup::Euclinf(objects2, p, r);
      }
      else {
         Rcpp::stop("`distance` should be one of: \"levenshtein\" (default), \"dinu\", \"hamming\", \"euclinf\"");
      }
   }
   else if (Rf_isString(objects) && (Rf_isNull(distance) || Rf_isString(distance)))
   {
      Rcpp::CharacterVector objects2(objects);
      Rcpp::CharacterVector distance2 =
         ((Rf_isNull(distance))?Rcpp::CharacterVector("levenshtein"):Rcpp::CharacterVector(distance));

      const char* distance3 = CHAR(STRING_ELT((SEXP)distance2, 0));
      if (!strcmp(distance3, "levenshtein")) {
         return (grup::Distance*)new grup::LevenshteinDistanceChar(objects2);
      }
      else if (!strcmp(distance3, "dinu")) {
         return (grup::Distance*)new grup::DinuDistanceChar(objects2);
      }
      else if (!strcmp(distance3, "hamming")) {
         return (grup::Distance*)new grup::HammingDistanceChar(objects2);
      }
      else {
         Rcpp::stop("`distance` should be one of: \"levenshtein\" (default), \"dinu\", \"hamming\"");
      }
   }
   else if (Rf_isMatrix(objects) && Rf_isNumeric(objects) && (Rf_isNull(distance) || Rf_isString(distance)))
   {
      Rcpp::NumericMatrix objects2(objects);
      Rcpp::CharacterVector distance2 =
         ((Rf_isNull(distance))?Rcpp::CharacterVector("euclidean_squared"):Rcpp::CharacterVector(distance));

      const char* distance3 = CHAR(STRING_ELT((SEXP)distance2, 0));
      if (!strcmp(distance3, "euclidean_squared")) {
         return (grup::Distance*)
            new grup::SquaredEuclideanDistance(
               objects2
            );
      }
      else if (!strcmp(distance3, "euclidean")) {
         return (grup::Distance*)
            new grup::EuclideanDistance(
               objects2
            );
      }
      else if (!strcmp(distance3, "manhattan")) {
         return (grup::Distance*)
            new grup::ManhattanDistance(
               objects2
            );
      }
      else if (!strcmp(distance3, "maximum")) {
         return (grup::Distance*)
            new grup::MaximumDistance(
               objects2
            );
      }
      else if (!strcmp(distance3, "hamming")) {
         return (grup::Distance*)
            new grup::HammingDistance(
               objects2
            );
      }
      else {
         Rcpp::stop("`distance` should be one of: \"euclidean_squared\" (default), \"euclidean\", \"manhattan\", \"maximum\", \"hamming\"");
      }
   }
   else {
      Rcpp::stop("incorrect input data");
   }
   return NULL;
}



GenericMatrixDistance::GenericMatrixDistance(const Rcpp::NumericMatrix& points) :
      Distance(points.nrow()),
      items(REAL((SEXP)points)), m(points.ncol())  {
   // act on a transposed matrix to avoid many L1/L... cache misses
   items = new double[m*n];
   const double* items2 = REAL((SEXP)points);
   double* items_ptr = items;
   for (size_t i=0; i<n; ++i) {
      for (size_t j=0; j<m; ++j) {
         if (!std::isfinite(items2[j*n+i]))
            Rcpp::stop("missing values and infinities in input objects are not allowed");
         *(items_ptr++) = items2[j*n+i];
      }
   }
}

double SquaredEuclideanDistance::compute(size_t v1, size_t v2)
{
   if (v1 == v2) return 0.0;
   double d = 0.0;
   for (size_t i=0; i<m; ++i)
      d += (items[v1*m+i]-items[v2*m+i])*(items[v1*m+i]-items[v2*m+i]);
   return d;
}


double EuclideanDistance::compute(size_t v1, size_t v2)
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


double ManhattanDistance::compute(size_t v1, size_t v2)
{
   if (v1 == v2) return 0.0;
   double d = 0.0;
   for (size_t i=0; i<m; ++i)
      d += std::abs(items[v1*m+i]-items[v2*m+i]);
   return d;
}


double MaximumDistance::compute(size_t v1, size_t v2)
{
   if (v1 == v2) return 0.0;
   double d = 0.0;
   for (size_t i=0; i<m; ++i) {
      double d2 = std::abs(items[v1*m+i]-items[v2*m+i]);
      if (d2 > d) d = d2;
   }
   return d;
}

double HammingDistance::compute(size_t v1, size_t v2)
{
   if (v1 == v2) return 0.0;
   double d = 0.0;
   for (size_t i=0; i<m; ++i) {
      if (items[v1*m+i] != items[v2*m+i]) d += 1.0;
   }
   return d;
}


double GenericRDistance::compute(size_t v1, size_t v2)
{
   if (v1 == v2) return 0.0;
   return ((Rcpp::NumericVector)distfun(items[v1], items[v2]))[0];
}




double DistObjectDistance::compute(size_t v1, size_t v2)
{
   if (v1 == v2) return 0.0;

   size_t i;
   if (v1 < v2)
      i = n*v1-((v1+1)*(v1))/2+v2-v1-1;
   else
      i = n*v2-((v2+1)*(v2))/2+v1-v2-1;

   return items[i];
}




// --------------------------------------------------------------------------------------------


StringDistanceInt::StringDistanceInt(const Rcpp::List& strings) :
   Distance(strings.size()),
      robj()
{
   R_PreserveObject(robj = (SEXP)strings);

   items = new const int*[n];
   lengths = new size_t[n];

   for (size_t i=0; i<n; ++i) {
      SEXP cur = VECTOR_ELT(robj, i);
      if (!Rf_isInteger(cur))
         Rcpp::stop("only integer vectors are allowed in the input list; check for NULLs, NAs, etc.");
      lengths[i] = LENGTH(cur);
      items[i] = INTEGER(cur);

      for (size_t j=0; j<lengths[i]; ++j)
         if (items[i][j] == NA_INTEGER)
            Rcpp::stop("missing values in input objects are not allowed");
   }
}

StringDistanceInt::~StringDistanceInt() {
   delete [] items;
   delete [] lengths;
   R_ReleaseObject(robj);
}



StringDistanceChar::StringDistanceChar(const Rcpp::CharacterVector& strings) :
   Distance(strings.size()),
      robj()
{
   R_PreserveObject(robj = (SEXP)strings);

   items = new const char*[n];
   lengths = new size_t[n];

   for (size_t i=0; i<n; ++i) {
      SEXP cur = STRING_ELT(robj, i);
      if (cur == NA_STRING)
         Rcpp::stop("missing values are not allowed");
      // if (Rf_getCharCE(cur) != CE_ANY)
         // Rcpp::stop("only ASCII strings allowed. Try with stringi::stri_enc_toutf32()");
      lengths[i] = LENGTH(cur);
      items[i] = CHAR(cur);
   }
}

StringDistanceChar::~StringDistanceChar() {
   delete [] items;
   delete [] lengths;
   R_ReleaseObject(robj);
}


StringDistanceDouble::StringDistanceDouble(const Rcpp::List& vectors) :
   Distance(vectors.size()),
   robj()
{
   R_PreserveObject(robj = (SEXP)vectors);
   items = new const double*[n];
   lengths = new size_t[n];

   for (size_t i=0; i<n; ++i) {
      SEXP cur = VECTOR_ELT(robj, i);
      if (!Rf_isReal(cur))
      Rcpp::stop("only real vectors are allowed in the input list; check for NULLs, NAs, etc.");
      lengths[i] = LENGTH(cur);
      items[i] = REAL(cur);

      for (size_t j=0; j<lengths[i]; ++j)
         if (items[i][j] == NA_REAL)
            Rcpp::stop("missing values in input objects are not allowed");
   }
}


StringDistanceDouble::~StringDistanceDouble() {
   delete [] items;
   delete [] lengths;
   R_ReleaseObject(robj);
}



// --------------------------------------------------------------------------------------------


template<class T> double distance_levenshtein(const T* s1, const T* s2, size_t n1, size_t n2) {
  if (n1 < n2) {
      std::swap(s1, s2); // pointer swap
      std::swap(n1, n2);
   }

// #ifdef _OPENMP
   // to be thread-safe, we have to allocate these 2 arrays each time...
   size_t* v_cur = new size_t[n2+1];
   size_t* v_last = new size_t[n2+1];
// #endif

   // n2 <= n1
   for (size_t j=0; j<=n2; ++j) v_cur[j] = j;

   for (size_t i=1; i<=n1; ++i) {
      std::swap(v_last, v_cur); // pointer swap
      v_cur[0] = i;
      for (size_t j=1; j<=n2; ++j) {
         if (s1[i-1] == s2[j-1])
            v_cur[j] = v_last[j-1];
         else
            v_cur[j] = std::min(std::min(
               v_last[j-1]+1,
               v_cur[j-1]+1),
               v_last[j]+1);
      }
   }

   double ret = (double) v_cur[n2];
// #ifdef _OPENMP
   delete [] v_cur;
   delete [] v_last;
// #endif
   return ret;
}


double LevenshteinDistanceInt::compute(size_t v1, size_t v2)
{
   return distance_levenshtein(items[v1], items[v2], lengths[v1], lengths[v2]);
}

double LevenshteinDistanceChar::compute(size_t v1, size_t v2)
{
   return distance_levenshtein(items[v1], items[v2], lengths[v1], lengths[v2]);
}


// --------------------------------------------------------------------------------------------



template<class T> double distance_hamming(const T* s1, const T* s2, size_t n1, size_t n2) {
   if (n1 != n2)
      Rcpp::stop("objects should be of the same dimension");

   double d = 0.0;
   for (size_t i=0; i<n1; ++i) {
      if (s1[i] != s2[i]) d += 1.0;
   }

   return d;
}


double HammingDistanceInt::compute(size_t v1, size_t v2)
{
   return distance_hamming(items[v1], items[v2], lengths[v1], lengths[v2]);
}

double HammingDistanceChar::compute(size_t v1, size_t v2)
{
   return distance_hamming(items[v1], items[v2], lengths[v1], lengths[v2]);
}

// --------------------------------------------------------------------------------------------

template<class T> double distance_dinu(const T* x, const T* y, const size_t* ox, const size_t* oy, size_t nx, size_t ny) {
   double d = 0.0;
   size_t ix = 0, iy = 0;
   while (ix < nx && iy < ny) {
      if (x[ox[ix]] == y[oy[iy]])
         d += std::abs((ox[ix++]+1.0) - (oy[iy++]+1.0));
      else if (x[ox[ix]] < y[oy[iy]])
         d += std::abs((ox[ix++]+1.0) - 0.0);
      else
         d += std::abs(0.0 - (oy[iy++]+1.0));
   }
   while (ix < nx) d += std::abs((ox[ix++]+1.0) - 0.0);
   while (iy < ny) d += std::abs(0.0 - (oy[iy++]+1.0));

   return d;
}


double DinuDistanceInt::compute(size_t v1, size_t v2)
{
   const int* x = items[v1];
   const int* y = items[v2];
   const size_t* ox = ranks[v1].data();
   const size_t* oy = ranks[v2].data();
   size_t nx = lengths[v1];
   size_t ny = lengths[v2];

   return distance_dinu(x, y, ox, oy, nx, ny);
}


double DinuDistanceChar::compute(size_t v1, size_t v2)
{
   const char* x = items[v1];
   const char* y = items[v2];
   const size_t* ox = ranks[v1].data();
   const size_t* oy = ranks[v2].data();
   size_t nx = lengths[v1];
   size_t ny = lengths[v2];

   return distance_dinu(x, y, ox, oy, nx, ny);
}


double Euclinf::compute(size_t v1, size_t v2)
{
  const double* x = items[v1];
  const double* y = items[v2];
  size_t nx = lengths[v1];
  size_t ny = lengths[v2];
  double dist = 0.0;
  std::size_t min_nx_ny = std::min(nx, ny);
  for (std::size_t i=0; i<min_nx_ny;  ++i) dist += (x[i]-y[i])*(x[i]-y[i]);
  for (std::size_t i=min_nx_ny; i<nx; ++i) dist += x[i]*x[i];
  for (std::size_t i=min_nx_ny; i<ny; ++i) dist += y[i]*y[i];
  dist += p*fabs(std::pow(nx, r)-std::pow(ny, r));
  return dist;
}


