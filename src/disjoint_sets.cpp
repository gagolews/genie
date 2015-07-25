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


#include "disjoint_sets.h"

using namespace DataStructures;


DisjointSets::DisjointSets(std::size_t n) :
   clusterParent(std::vector< std::size_t >(n)),
   n(n)
{
#ifdef DISJOINT_SETS_DEBUG
   Rcpp::Rcout << "Warning: DISJOINT_SETS_DEBUG defined in disjoint_sets.h\n";
#endif
   for (std::size_t i=0; i<n; ++i)
      clusterParent[i] = i;
}


DisjointSets::~DisjointSets() {
#ifdef DISJOINT_SETS_DEBUG
   std::size_t px = find_set(0);
   for (size_t i=1; i<n; ++i) {
      if (px != find_set(i)) Rcpp::stop("~DisjointSets: assert failed");
   }
#endif
}


std::size_t DisjointSets::link(std::size_t x, std::size_t y) {
#ifdef DISJOINT_SETS_DEBUG
   if (clusterParent[y] != y || clusterParent[x] != x)
      Rcpp::stop("DisjointSets::link assert failed");
#endif
   return (clusterParent[x] = y);
}


std::size_t DisjointSets::union_set(std::size_t x, std::size_t y) {
   // will call overriden link() method
   return link(find_set(x), find_set(y));
}


PhatDisjointSets::PhatDisjointSets(std::size_t n) :
   DisjointSets(n),
   clusterSize(std::vector< std::size_t >(n, 1)),
   clusterCount(n),
   clusterMembers(std::vector< std::list<std::size_t> >(n)),
   clusterNext(std::vector< std::size_t >(n)),
   clusterPrev(std::vector< std::size_t >(n))
{
   for (std::size_t i=0; i<n; ++i) {
      clusterMembers[i].push_front(i);
      clusterNext[i] = (i<n-1)?(i+1):0;
      clusterPrev[i] = (i>0)?(i-1):(n-1);
   }
}


PhatDisjointSets::~PhatDisjointSets()
{
#ifdef DISJOINT_SETS_DEBUG
   if (getClusterCount() != 1 ||
         getClusterSize(find_set(0)) != n ||
         getClusterNext(find_set(0)) != find_set(0) ||
         getClusterPrev(find_set(0)) != find_set(0)
      )
      Rcpp::stop("~PhatDisjointSets: assert failed");
#endif
}


std::size_t PhatDisjointSets::link(std::size_t x, std::size_t y)
{
   std::size_t cx = clusterSize[x];
   std::size_t z = DisjointSets::link(x, y);
#ifdef DISJOINT_SETS_DEBUG
   if (z != y)
      Rcpp::stop("DisjointSets::link assert failed");
#endif

   std::size_t oldprev = clusterPrev[x];
   std::size_t oldnext = clusterNext[x];
   clusterPrev[ oldnext ] = oldprev;
   clusterNext[ oldprev ] = oldnext;

   clusterSize[z] += cx;

   clusterMembers[z].splice(clusterMembers[z].end(), clusterMembers[x]); // O(1)

   --clusterCount;
   return z;
}
