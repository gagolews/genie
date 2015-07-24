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


#ifndef __DISJOINT_SETS_H
#define __DISJOINT_SETS_H

#include <vector>
#include <cstdint>

namespace DataStructures {

class DisjointSets {
private:
   std::vector< std::size_t > clusterParent;
   std::size_t n;

public:
   DisjointSets(std::size_t n) :
      clusterParent(std::vector< std::size_t >(n)),
      n(n)
   {
      for (std::size_t i=0; i<n; ++i)
         clusterParent[i] = i;
   }

   ~DisjointSets()
   {

   }

   inline std::size_t find_set(std::size_t x) {
      if (clusterParent[x] != x)
         return clusterParent[x] = find_set(clusterParent[x]);
      else
         return clusterParent[x];
   }

   inline std::size_t link(std::size_t x, std::size_t y) {
      // we have: clusterParent[y] == y
      // we have: clusterParent[x] == x
      return (clusterParent[x] = y);
   }

   inline std::size_t union_set(std::size_t x, std::size_t y) {
      return link(find_set(x), find_set(y));
   }
};

} /* namespace DataStructures */

#endif
