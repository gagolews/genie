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


#ifndef __DISJOINT_SETS_H
#define __DISJOINT_SETS_H

#include "defs.h"
#include <vector>
#include <list>
#include <cstdint>
#include <Rcpp.h>

/* see defs.h */
#ifndef DISJOINT_SETS_DEBUG
#define DISJOINT_SETS_DEBUG_CONST const
#else
#define DISJOINT_SETS_DEBUG_CONST /* const */
#endif

namespace grup {

class DisjointSets {
private:
   std::vector< std::size_t > clusterParent;

protected:
   std::size_t n;

public:
   DisjointSets(std::size_t n);
   virtual ~DisjointSets();

   virtual std::size_t link(std::size_t x, std::size_t y, std::size_t z);
   virtual std::size_t link(std::size_t x, std::size_t y);
   std::size_t union_set(std::size_t x, std::size_t y);

   inline std::size_t find_set(std::size_t x) {
      if (clusterParent[x] != x)
         return clusterParent[x] = find_set(clusterParent[x]);
      else
         return clusterParent[x];
   }

};


class PhatDisjointSets : public DisjointSets {
private:
   std::vector< std::size_t > clusterSize;
   std::vector< std::size_t* > clusterMembers;
   std::vector< std::size_t > clusterNext;
   std::vector< std::size_t > clusterPrev;
   std::size_t clusterCount;
   std::size_t minClusterSize;
   std::size_t minClusterCount;

   void recomputeMinClusterSize();

public:
   PhatDisjointSets(std::size_t n);
   virtual ~PhatDisjointSets();

   virtual std::size_t link(std::size_t x, std::size_t y);
   virtual std::size_t link(std::size_t x, std::size_t y, std::size_t z);

   inline std::size_t getClusterCount() const { return clusterCount; }
   inline std::size_t getMinClusterSize() const { return minClusterSize; }

   inline const std::size_t* getClusterMembers(std::size_t x) DISJOINT_SETS_DEBUG_CONST {
      #ifdef DISJOINT_SETS_DEBUG
      STOPIFNOT(find_set(x) == x);
      STOPIFNOT(clusterMembers[x]);
      #endif
      return clusterMembers[x];
   }

   inline std::size_t getClusterSize(std::size_t x) DISJOINT_SETS_DEBUG_CONST {
      #ifdef DISJOINT_SETS_DEBUG
      STOPIFNOT(find_set(x) == x);
      STOPIFNOT(clusterSize[x] == 0 || clusterMembers[x] != NULL);
      #endif
      return clusterSize[x];
   }

   inline std::size_t getClusterPrev(std::size_t x) DISJOINT_SETS_DEBUG_CONST {
      #ifdef DISJOINT_SETS_DEBUG
      STOPIFNOT(find_set(x) == x);
      STOPIFNOT(find_set(clusterPrev[x]) == clusterPrev[x]);
      STOPIFNOT(find_set(clusterNext[x]) == clusterNext[x]);
      #endif
      return clusterPrev[x];
   }

   inline std::size_t getClusterNext(std::size_t x) DISJOINT_SETS_DEBUG_CONST {
      /*
      to iterate over all clusters starting from x, use something like:
      for (size_t nx = ds.getClusterNext(x); nx != x; nx = ds.getClusterNext(nx)) {
         // e.g.:
         for (auto it = ds.getClusterMembers(nx).cbegin(); it != ds.getClusterMembers(nx).cend(); ++it)
            // play with *it
      }
      */
      #ifdef DISJOINT_SETS_DEBUG
      STOPIFNOT(find_set(x) == x);
      STOPIFNOT(find_set(clusterPrev[x]) == clusterPrev[x]);
      STOPIFNOT(find_set(clusterNext[x]) == clusterNext[x]);
      #endif
      return clusterNext[x];
   }
};

} /* namespace grup */

#endif
