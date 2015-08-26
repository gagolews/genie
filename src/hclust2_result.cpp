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

#include "hclust2_result.h"
using namespace DataStructures;
using namespace Rcpp;


#pragma message "@TODO: HClustResult::order is not yet generated"


/*
 * The algorithm used in hclust is to order the subtree so that the tighter
 *  cluster is on the left (the last, i.e., most recent, merge of the left
 *   subtree is at a lower value than the last merge of the right subtree).
 *   Single observations are the tightest clusters possible, and merges
 *   involving two observations place them in order by their
 *   observation sequence number.
 */

/*
 * If you can instead create an object like 'hcl', then pass the results to
.Fortran("hcass2",...), you should be done.

.Fortran(stats:::C_hcass2)  -- that won't work on CRAN...
 */

/*
 * http://lib.stat.cmu.edu/S/multiv
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                               C
C  Given a HIERARCHIC CLUSTERING, described as a sequence of    C
C  agglomerations, prepare the seq. of aggloms. and "horiz."    C
C  order of objects for plotting the dendrogram using S routine C
C  'plclust'.                                                   C
C                                                               C
C  Parameters:                                                  C
C                                                               C
C  IA, IB:       vectors of dimension N defining the agglomer-  C
C                 ations.                                       C
C  IIA, IIB:     used to store IA and IB values differently     C
C                (in form needed for S command `plclust`        C
C  IORDER:       "horiz." order of objects for dendrogram       C
C                                                               C
C  F. Murtagh, ESA/ESO/STECF, Garching, June 1991               C
C                                                               C
C  HISTORY                                                      C
C                                                               C
C  Adapted from routine HCASS, which additionally determines    C
C   cluster assignments at all levels, at extra comput. expense C
C                                                               C
C---------------------------------------------------------------C
      SUBROUTINE HCASS2(N,IA,IB,IORDER,IIA,IIB)
      INTEGER IA(N),IB(N),IORDER(N),IIA(N),IIB(N)
C
C     Following bit is to get seq. of merges into format acceptable to plclust
C     I coded clusters as lowest seq. no. of constituents; S's `hclust' codes
C     singletons as -ve numbers, and non-singletons with their seq. nos.
C
      DO 912 I=1,N
         IIA(I)=IA(I)
         IIB(I)=IB(I)
  912 CONTINUE
      DO 915 I=1,N-2
C        In the following, smallest (+ve or -ve) seq. no. wanted
         K=MIN(IA(I),IB(I))
         DO 913 J=I+1, N-1
            IF(IA(J).EQ.K) IIA(J)=-I
            IF(IB(J).EQ.K) IIB(J)=-I
  913    CONTINUE
  915 CONTINUE
      DO 916 I=1,N-1
         IIA(I)=-IIA(I)
         IIB(I)=-IIB(I)
  916 CONTINUE
      DO 917 I=1,N-1
         IF (IIA(I).GT.0.AND.IIB(I).LT.0) THEN
            K = IIA(I)
            IIA(I) = IIB(I)
            IIB(I) = K
         ENDIF
         IF (IIA(I).GT.0.AND.IIB(I).GT.0) THEN
            K1 = MIN(IIA(I),IIB(I))
            K2 = MAX(IIA(I),IIB(I))
            IIA(I) = K1
            IIB(I) = K2
         ENDIF
  917 CONTINUE
C
C
C     NEW PART FOR `ORDER'
C
      IORDER(1) =IIA(N-1)
      IORDER(2) =IIB(N-1)
      LOC=2
      DO 175 I=N-2,1,-1
        DO 169 J=1,LOC
          IF(IORDER(J).EQ.I) THEN
C           REPLACE IORDER(J) WITH IIA(I) AND IIB(I)
            IORDER(J)=IIA(I)
            IF (J.EQ.LOC) THEN
                LOC=LOC+1
                IORDER(LOC)=IIB(I)
                GOTO 171
            ENDIF
            LOC=LOC+1
            DO 95 K=LOC,J+2,-1
               IORDER(K)=IORDER(K-1)
  95        CONTINUE
            IORDER(J+1)=IIB(I)
            GOTO 171
          ENDIF
 169    CONTINUE
C       SHOULD NEVER REACH HERE
 171    CONTINUE
 175  CONTINUE
C
C
C
      DO 181 I=1,N
         IORDER(I) = -IORDER(I)
 181  CONTINUE
C
C
      RETURN
      END



 */

HClustResult::HClustResult(size_t n, Distance* dist) :
      curiter(0),
      n(n),
      links(n-1, 2),         // this may be meaningless for some methods, do not return
      merge(n-1, 2),         // TO DO: merge determination during link()
      height(n-1),
      order(n, NA_REAL),     // TO DO: how to implement that?
      labels(dist->getLabels()),
      dist_method(dist->getDistMethod())
{

}


void HClustResult::link(size_t i1, size_t i2, double d12) {
   STOPIFNOT(curiter < n-1);
   links(curiter, 0) = (double)i1;
   links(curiter, 1) = (double)i2;
   height(curiter) = d12;
   ++curiter;

   if (curiter == n-1) generateMergeMatrix();
}


void HClustResult::generateMergeMatrix() {
   STOPIFNOT(curiter == n-1);

   std::vector<size_t> elements(n+1, 0);
   std::vector<size_t> parents(n+1, 0);

   size_t clusterNumber = 1;
   for (size_t k=0; k<n-1; ++k, ++clusterNumber) {
      size_t i = (size_t)links(k, 0) + 1;
      size_t j = (size_t)links(k, 1) + 1;
      size_t si = elements[i];
      size_t sj = elements[j];
      elements[i] = clusterNumber;
      elements[j] = clusterNumber;

      if (si == 0)
         merge(k, 0) = -(double)i;
      else {
         while (parents[si] != 0) {
            size_t sinew = parents[si];
            parents[si] = clusterNumber;
            si = sinew;
         }
         if (si != 0) parents[si] = clusterNumber;
         merge(k,0) = (double)si;
      }

      if (sj == 0)
         merge(k, 1) = -(double)j;
      else {
         while (parents[sj] != 0) {
            size_t sjnew = parents[sj];
            parents[sj] = clusterNumber;
            sj = sjnew;
         }
         if (sj != 0) parents[sj] = clusterNumber;
         merge(k,1) = (double)sj;
      }
   }
}


List HClustResult::toR(
      const HClustStats& hclustStats,
      const HClustOptions& hclustOptions,
      const DistanceStats& distStats)
{
   MESSAGE_2("[%010.3f] generating output matrix\n", clock()/(float)CLOCKS_PER_SEC);

   List result = List::create(
      _["merge"]  = merge,
      _["height"] = height,
      _["order"]  = order,
      _["labels"] = labels,
      _["call"]   = R_NilValue,
      _["method"] = R_NilValue,
      _["dist.method"] = dist_method,
      // _["links"]  = links,
      _["stats"] = List::create(
         _["method"] = hclustStats.toR(),
         _["distance"] = distStats.toR()
      ),
      _["control"] = List::create(
         _["method"] = hclustOptions.toR()
      )
   );

   result.attr("class") = "hclust";
   return result;
}
