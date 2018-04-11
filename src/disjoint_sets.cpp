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


#include "disjoint_sets.h"

using namespace grup;


DisjointSets::DisjointSets(std::size_t _n) :
   clusterParent(std::vector< std::size_t >(_n)),
   n(_n)
{
#ifdef DISJOINT_SETS_DEBUG
   MESSAGE_1("Warning: DISJOINT_SETS_DEBUG defined in disjoint_sets.h\n");
#endif
   for (std::size_t i=0; i<n; ++i)
      clusterParent[i] = i;
}


DisjointSets::~DisjointSets() {
#ifdef DISJOINT_SETS_DEBUG
   std::size_t px = find_set(0);
   for (size_t i=1; i<n; ++i) STOPIFNOT(px == find_set(i));
#endif
}


std::size_t DisjointSets::link(std::size_t x, std::size_t y, std::size_t z) {
#ifdef DISJOINT_SETS_DEBUG
   STOPIFNOT(find_set(z) == x || find_set(z) == y);
   STOPIFNOT(clusterParent[y] == y && clusterParent[x] == x);
   STOPIFNOT(x != y);
#endif

   clusterParent[z] = z;
   clusterParent[y] = z;
   clusterParent[x] = z;
   return z;
}

std::size_t DisjointSets::link(std::size_t x, std::size_t y) {
#ifdef DISJOINT_SETS_DEBUG
   STOPIFNOT(clusterParent[y] == y && clusterParent[x] == x);
   STOPIFNOT(x != y);
#endif
   return (clusterParent[y] = x);
}


std::size_t DisjointSets::union_set(std::size_t x, std::size_t y) {
   // will call overriden link() method
   return link(find_set(x), find_set(y));
}


PhatDisjointSets::PhatDisjointSets(std::size_t _n) :
   DisjointSets(_n),
   clusterSize(std::vector< std::size_t >(_n, 1)),
   clusterMembers(_n),
   clusterNext(std::vector< std::size_t >(_n)),
   clusterPrev(std::vector< std::size_t >(_n))
{
   clusterCount = _n;
   minClusterSize = 1;
   minClusterCount = _n;
   for (std::size_t i=0; i<_n; ++i) {
      clusterMembers[i] = (std::size_t*)malloc(sizeof(std::size_t)*1);
      clusterMembers[i][0] = i;
      // clusterSize[i] = 1; // already initialized
      clusterNext[i] = (i<_n-1)?(i+1):0;
      clusterPrev[i] = (i>0)?(i-1):(_n-1);
   }
}


PhatDisjointSets::~PhatDisjointSets()
{
#ifdef DISJOINT_SETS_DEBUG
   STOPIFNOT(getClusterCount() == 1);
   STOPIFNOT(getClusterSize(find_set(0)) == n);
   STOPIFNOT(getClusterNext(find_set(0)) == find_set(0));
   STOPIFNOT(getClusterPrev(find_set(0)) == find_set(0));
#endif
   for (std::size_t i=0; i<n; ++i)
      if (clusterMembers[i]) free(clusterMembers[i]);
}


std::size_t PhatDisjointSets::link(std::size_t x, std::size_t y)
{
   std::size_t sizex = getClusterSize(x);
   std::size_t sizey = getClusterSize(y);
   std::size_t z = DisjointSets::link(x, y);
#ifdef DISJOINT_SETS_DEBUG
   STOPIFNOT(z == x);
   STOPIFNOT(clusterMembers[x]);
   STOPIFNOT(clusterMembers[y]);
#endif

   if (clusterCount > 2) {
      std::size_t oldprev = clusterPrev[y];
      std::size_t oldnext = clusterNext[y];
      clusterPrev[ oldnext ] = oldprev;
      clusterNext[ oldprev ] = oldnext;
   }
   else {
      clusterPrev[ z ] = z;
      clusterNext[ z ] = z;
   }

   clusterMembers[z] = (std::size_t*)realloc(clusterMembers[z], (clusterSize[x]+clusterSize[y])*sizeof(std::size_t));
   memcpy(clusterMembers[z]+clusterSize[x], clusterMembers[y], clusterSize[y]*sizeof(std::size_t));
   free(clusterMembers[y]);
   clusterMembers[y] = NULL;

   clusterSize[z] += clusterSize[y];
   // clusterMembers[z]->splice(clusterMembers[z]->end(), *(clusterMembers[y])); // O(1)

   --clusterCount;
   if (minClusterCount > 0 && sizex == minClusterSize) minClusterCount--;
   if (minClusterCount > 0 && sizey == minClusterSize) minClusterCount--;
   if (minClusterCount == 0) recomputeMinClusterSize();

   return z;
}


std::size_t PhatDisjointSets::link(std::size_t x, std::size_t y, std::size_t z)
{
   std::size_t z2 = DisjointSets::link(x, y, z);
#ifdef DISJOINT_SETS_DEBUG
   STOPIFNOT(z == z2);
#endif

   if (clusterCount > 2) {
      std::size_t oldprev = clusterPrev[y];
      std::size_t oldnext = clusterNext[y];
      clusterPrev[ oldnext ] = oldprev;
      clusterNext[ oldprev ] = oldnext;

      oldprev = clusterPrev[x];
      oldnext = clusterNext[x];
      clusterPrev[z2] = oldprev;
      clusterNext[z2] = oldnext;
   #ifdef DISJOINT_SETS_DEBUG
      STOPIFNOT(clusterPrev[ oldnext ] == x);
      STOPIFNOT(clusterNext[ oldprev ] == x);
   #endif
      clusterPrev[ oldnext ] = z2;
      clusterNext[ oldprev ] = z2;
   }
   else {
      clusterPrev[ z2 ] = z2;
      clusterNext[ z2 ] = z2;
   }


//    clusterMembers[x]->splice(clusterMembers[x]->end(), (*clusterMembers[y])); // O(1)
//    delete clusterMembers[y];
//    clusterMembers[y] = NULL;
// #ifdef DISJOINT_SETS_DEBUG
//    if (x != z2 && clusterMembers[z2])
//       Rcpp::stop("PhatDisjointSets::link assert failed");
// #endif
//    std::swap(clusterMembers[z2], clusterMembers[x]);
//    clusterSize[z2] = clusterSize[x] + clusterSize[y];

   clusterMembers[x] = (std::size_t*)realloc(clusterMembers[x], (clusterSize[x]+clusterSize[y])*sizeof(std::size_t));
   memcpy(clusterMembers[x]+clusterSize[x], clusterMembers[y], clusterSize[y]*sizeof(std::size_t));
   free(clusterMembers[y]);
   clusterMembers[y] = NULL;
#ifdef DISJOINT_SETS_DEBUG
   STOPIFNOT(x == z2 || !clusterMembers[z2]);
#endif
   std::swap(clusterMembers[z2], clusterMembers[x]);
#ifdef DISJOINT_SETS_DEBUG
   STOPIFNOT(x == z2 || (!clusterMembers[x] && clusterMembers[z2]));
#endif

   clusterSize[z2] = clusterSize[x] + clusterSize[y];

   --clusterCount;
   return z2;
}



void PhatDisjointSets::recomputeMinClusterSize() {
   std::size_t start = find_set(0);
   minClusterSize = getClusterSize(start);
   minClusterCount = 1;
   size_t curi = getClusterNext(start);
   while (curi != start) {
      std::size_t curSize = getClusterSize(curi);
      if (curSize == minClusterSize)
         minClusterCount++;
      else if (curSize < minClusterSize) {
         minClusterSize = curSize;
         minClusterCount = 1;
      }
      curi = getClusterNext(curi);
   }
}

