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


#include "hclust2_single.h"



/* improvement idea:
 * add custom sort of _indices
 *
 * useful for Levenshtein distance
 * long strings should be put at the end
 * ?
 */


using namespace Rcpp;
using namespace std;
using namespace boost;
using namespace DataStructures;


// constructor (OK, we all know what this is, but I label it for faster in-code search)
HClustBiVpTreeSingle::HClustBiVpTreeSingle(Distance* dist, size_t maxNumberOfElementsInLeaves) :
   maxNumberOfElementsInLeaves(maxNumberOfElementsInLeaves),
   _root(NULL), _n(dist->getObjectCount()), _distance(dist),
   _indices(dist->getObjectCount()),
   neighborsCount(vector<size_t>(dist->getObjectCount(), 0)),
   minRadiuses(vector<double>(dist->getObjectCount(), -INFINITY)),
   // maxRadiuses(vector<double>(dist->getObjectCount(), INFINITY)),
   shouldFind(vector<bool>(dist->getObjectCount(), true)),
   nearestNeighbors(vector< deque<HeapNeighborItem> >(dist->getObjectCount())),
#ifdef GENERATE_STATS
   stats(HClustBiVpTreeStats()),
#endif
#ifdef USE_BOOST_DISJOINT_SETS
   ds(make_assoc_property_map(rank), make_assoc_property_map(parent))
#else
   ds(dist->getObjectCount())
#endif
{
#ifdef USE_BOOST_DISJOINT_SETS
   Rcpp::Rcout << "Warning: USE_BOOST_DISJOINT_SETS is defined in hclust2_single.h\n";
#endif
#if VERBOSE > 5
   Rprintf("[%010.3f] building vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   // maxNumberOfElementsInLeaves = 2; //(size_t)log2(_n);

   // starting _indices: random permutation of {0,1,...,_n-1}
   for(size_t i=0;i<_n;i++) _indices[i] = i;
   for(size_t i=_n-1; i>= 1; i--)
      swap(_indices[i], _indices[(size_t)(unif_rand()*(i+1))]);

   for(size_t i=0; i<_n; i++)
   {
#ifdef USE_BOOST_DISJOINT_SETS
     ds.make_set(i);
#endif
#ifdef MB_IMPROVEMENT
     clusters.emplace(i);
#endif // MB_IMPROVEMENT
   }

   _root = buildFromPoints(0, _n);
}


HClustBiVpTreeSingle::~HClustBiVpTreeSingle() {
// #if VERBOSE > 5
//       Rprintf("[%010.3f] destroying vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
// #endif
   if(_root) delete _root;
}


int HClustBiVpTreeSingle::chooseNewVantagePoint(size_t left, size_t right)
{
#if VANTAGE_POINT_SELECT_SCHEME == 1
   // idea by A. Fu et al., "Dynamic VP-tree indexing for n-nearest neighbor
   //    search given pair-wise distances"
   size_t numCandidates = VANTAGE_POINT_SELECT_SCHEME_1_NUMCANDIDATES;
   size_t numTest = VANTAGE_POINT_SELECT_SCHEME_1_NUMTEST;

   if (left + numCandidates + numTest > right )
      return left;

   // maximize variance
   size_t bestIndex = -1;
   double bestSigma = -INFINITY;
   for(size_t i=left; i<left+numCandidates; i++) {
      accumulators::accumulator_set< double,
         accumulators::features<accumulators::tag::variance> > acc;
      for (size_t j = left+numCandidates; j < left+numCandidates+numTest; ++j)
         acc( (*_distance)( _indices[i], _indices[j] ) );
      double curSigma = accumulators::variance(acc);
      if (curSigma > bestSigma) {
         bestSigma = curSigma;
         bestIndex = i;
      }
   }

   return bestIndex;
#elif VANTAGE_POINT_SELECT_SCHEME == 2
   // idea by T. Bozkaya and M. Ozsoyoglu, "Indexing large metric spaces
   //      for similarity search queries"

   // which one maximizes dist to _indices[left]?
   size_t bestIndex = left;
   double bestDist  = 0.0;
   for (size_t i=left+1; i<right; ++i) {
      double curDist = (*_distance)(_indices[left], _indices[i]);
      if (curDist > bestDist) {
         bestDist = curDist;
         bestIndex = i;
      }
   }
//       for (size_t i=left+2; i<right; ++i) {
//          double curDist = (*_distance)(_indices[left+1], _indices[i]);
//          if (curDist > bestDist) {
//             bestDist = curDist;
//             bestIndex = i;
//          }
//       }
//       for (size_t i=left+3; i<right; ++i) {
//          double curDist = (*_distance)(_indices[left+2], _indices[i]);
//          if (curDist > bestDist) {
//             bestDist = curDist;
//             bestIndex = i;
//          }
//       }
   return bestIndex;
#else
   // return random index == left one (sample is randomized already)
   return left;
#endif
}


HClustBiVpTreeNode* HClustBiVpTreeSingle::buildFromPoints(size_t left, size_t right)
{
#ifdef GENERATE_STATS
      ++stats.nodeCount;
#endif
   if(right - left <= maxNumberOfElementsInLeaves)
   {
      // for (size_t i=left; i<right; ++i) {
         // size_t j = _indices[(i+1 < right)?(i+1):left];
         // if (_indices[i] < j)
            // maxRadiuses[ _indices[i] ] = (*_distance)(_indices[i], j);
      // }

      return new HClustBiVpTreeNode(left, right);
   }

   size_t vpi_idx = chooseNewVantagePoint(left, right);
   std::swap(_indices[left], _indices[vpi_idx]);
   size_t vpi = _indices[left];

// #if VERBOSE > 7
//       ((EuclideanDistance*)_distance)->isVP[vpi] = true;
// #endif

   size_t median = ( right + left - 1 ) / 2;
   std::nth_element(_indices.begin() + left + 1, _indices.begin() + median,  _indices.begin() + right,
                    DistanceComparator(vpi, _distance));
   // std::sort(_indices.begin() + left+1, _indices.begin() + right,
                    // DistanceComparator(vpi, &_distance ));
   // printf("(%d,%d,%d)\n", left, median, right);
   // for (int i=left; i<right; ++i) printf("%d, ", _indices[i]+1);
   // printf("\n");
   HClustBiVpTreeNode* node = new HClustBiVpTreeNode(vpi, (*_distance)(vpi, _indices[median]));


   size_t middle1 = std::partition(_indices.begin() + left,  _indices.begin() + median + 1,  IndexComparator(vpi)) - _indices.begin();
   size_t middle2 = std::partition(_indices.begin() + median + 1,  _indices.begin() + right, IndexComparator(vpi)) - _indices.begin();
   // printf("(%d,%d,%d,%d,%d)\n", left, middle1, median, middle2, right);
   // for (int i=left; i<right; ++i) printf("%d, ", _indices[i]+1);
   // printf("\n");


   if (middle1 - left > 0)     node->ll = buildFromPoints(left, middle1);
   if (median+1 - middle1 > 0) node->lr = buildFromPoints(middle1, median + 1);
   if (middle2 - median-1 > 0) node->rl = buildFromPoints(median + 1, middle2);
   if (right-middle2 > 0)      node->rr = buildFromPoints(middle2, right);

   return node;
}


void HClustBiVpTreeSingle::getNearestNeighborsFromMinRadiusRecursive( HClustBiVpTreeNode* node, size_t index,
   size_t clusterIndex, double minR, double& maxR,
   std::priority_queue<HeapNeighborItem>& heap )
{
   // search within (minR, maxR]
   // if (node == NULL) return; // this should not happen
#ifdef GENERATE_STATS
   ++stats.nodeVisit;
#endif
   if (node->sameCluster) {
      if (node->vpindex == SIZE_MAX) {
         if (ds.find_set(_indices[node->left]) == clusterIndex) return;
      } else {
         if (ds.find_set(node->vpindex) == clusterIndex) return;
      }
   }

   if (node->vpindex == SIZE_MAX) // leaf
   {
      if (node->sameCluster) {
#ifdef MB_IMPROVEMENT
      size_t s = SIZE_MAX;
      if (mbimprovement)
         s = ds.find_set(_indices[node->left]);
      std::unordered_map<SortedPoint,double>::const_iterator distToClusterIterator;
      if (mbimprovement)
         distToClusterIterator = distClust.find(SortedPoint(s, clusterIndex));
      double distToCluster = INFINITY;
      if (mbimprovement && distToClusterIterator != distClust.end())
         distToCluster = distToClusterIterator->second;
#endif // MB_IMPROVEMENT
         for (size_t i=node->left; i<node->right; i++)
         {
            if(index >= _indices[i]) continue;
            double dist2 = (*_distance)(index, _indices[i]);
            if (dist2 > maxR || dist2 <= minR) continue;

#ifdef MB_IMPROVEMENT
            if (mbimprovement && dist2 > distToCluster) {
               //Rcout << "odrzucam!" << endl;
               continue;
            }
#endif // MB_IMPROVEMENT

            if (heap.size() >= maxNearestNeighborPrefetch) {
               if (dist2 < maxR) {
                  while (!heap.empty() && heap.top().dist == maxR) {
                     heap.pop();
                  }
               }
            }
            heap.push( HeapNeighborItem(_indices[i], dist2) );
            maxR = heap.top().dist;
#ifdef MB_IMPROVEMENT
            if (mbimprovement)
               distClust.emplace(SortedPoint(s, clusterIndex), dist2);
#endif // MB_IMPROVEMENT
         }
      }
      else {
         size_t commonCluster = ds.find_set(_indices[node->left]);
         for (size_t i=node->left; i<node->right; i++)
         {
            size_t currentCluster = ds.find_set(_indices[i]);
#ifdef MB_IMPROVEMENT
            std::unordered_map<SortedPoint,double>::const_iterator distToClusterIterator;
            if (mbimprovement)
               distToClusterIterator = distClust.find(SortedPoint(currentCluster, clusterIndex));
            double distToCluster = INFINITY;
            if (mbimprovement && distToClusterIterator != distClust.end())
               distToCluster = distToClusterIterator->second;
#endif // MB_IMPROVEMENT
            if (currentCluster != commonCluster) commonCluster = SIZE_MAX;
            if (currentCluster == clusterIndex) continue;

            if (index >= _indices[i]) continue;

            double dist2 = (*_distance)(index, _indices[i]);
            if (dist2 > maxR || dist2 <= minR) continue;
#ifdef MB_IMPROVEMENT
            if (mbimprovement && dist2 > distToCluster) continue;
#endif // MB_IMPROVEMENT

            if (heap.size() >= maxNearestNeighborPrefetch) {
               if (dist2 < maxR) {
                  while (!heap.empty() && heap.top().dist == maxR) {
                     heap.pop();
                  }
               }
            }
            heap.push( HeapNeighborItem(_indices[i], dist2) );
            maxR = heap.top().dist;
#ifdef MB_IMPROVEMENT
            if (mbimprovement)
               distClust.emplace(SortedPoint(currentCluster, clusterIndex), dist2);
#endif // MB_IMPROVEMENT
         }
         if (commonCluster != SIZE_MAX) node->sameCluster = true;
      }
      return;
   }
   // else // not a leaf
   double dist = (*_distance)(node->vpindex, index);

// this MB's improvement is not well-tested:
//       if (dist < maxR && dist > minR && index < node->vpindex && ds.find_set(node->vpindex) != clusterIndex) {
//
//          heap.push( HeapNeighborItem(node->vpindex, dist) );
//          maxR = heap.top().dist;
//       }


   if ( dist < node->radius ) {
      if ( dist - maxR <= node->radius && dist + node->radius > minR ) {

         if(node->ll && index <= node->vpindex)
            getNearestNeighborsFromMinRadiusRecursive( node->ll, index, clusterIndex, minR, maxR, heap );
         if(node->lr)
            getNearestNeighborsFromMinRadiusRecursive( node->lr, index, clusterIndex, minR, maxR, heap );
      }

      if ( dist + maxR >= node->radius ) {
         if(node->rl && index <= node->vpindex)
            getNearestNeighborsFromMinRadiusRecursive( node->rl, index, clusterIndex, minR, maxR, heap );
         if(node->rr)
            getNearestNeighborsFromMinRadiusRecursive( node->rr, index, clusterIndex, minR, maxR, heap );
      }

   } else /* ( dist >= node->radius ) */ {
      if ( dist + maxR >= node->radius ) {
         if(node->rl && index <= node->vpindex)
            getNearestNeighborsFromMinRadiusRecursive( node->rl, index, clusterIndex, minR, maxR, heap );
         if(node->rr)
            getNearestNeighborsFromMinRadiusRecursive( node->rr, index, clusterIndex, minR, maxR, heap );
      }

      if ( dist - maxR <= node->radius && dist + node->radius > minR ) {
         if(node->ll && index <= node->vpindex)
            getNearestNeighborsFromMinRadiusRecursive( node->ll, index, clusterIndex, minR, maxR, heap );
         if(node->lr)
            getNearestNeighborsFromMinRadiusRecursive( node->lr, index, clusterIndex, minR, maxR, heap );
      }
   }

   if (   !node->sameCluster
      && (!node->ll || node->ll->sameCluster)
      && (!node->lr || node->lr->sameCluster)
      && (!node->rl || node->rl->sameCluster)
      && (!node->rr || node->rr->sameCluster)  )
   {
      size_t commonCluster = SIZE_MAX;
      if (node->ll) {
         size_t currentCluster = ds.find_set((node->ll->vpindex == SIZE_MAX)?_indices[node->ll->left]:node->ll->vpindex);
         if (commonCluster == SIZE_MAX) commonCluster = currentCluster;
         else if (currentCluster != commonCluster) return;
      }
      if (node->lr) {
         size_t currentCluster = ds.find_set((node->lr->vpindex == SIZE_MAX)?_indices[node->lr->left]:node->lr->vpindex);
         if (commonCluster == SIZE_MAX) commonCluster = currentCluster;
         else if (currentCluster != commonCluster) return;
      }
      if (node->rl) {
         size_t currentCluster = ds.find_set((node->rl->vpindex == SIZE_MAX)?_indices[node->rl->left]:node->rl->vpindex);
         if (commonCluster == SIZE_MAX) commonCluster = currentCluster;
         else if (currentCluster != commonCluster) return;
      }
      if (node->rr) {
         size_t currentCluster = ds.find_set((node->rr->vpindex == SIZE_MAX)?_indices[node->rr->left]:node->rr->vpindex);
         if (commonCluster == SIZE_MAX) commonCluster = currentCluster;
         else if (currentCluster != commonCluster) return;
      }
      node->sameCluster = true;
   }
}


HeapNeighborItem HClustBiVpTreeSingle::getNearestNeighbor(size_t index)
{
#if VERBOSE > 5
   // Rprintf(".");
#endif
//   This is not faster:
//       while(!nearestNeighbors[index].empty())
//       {
//          size_t sx = ds.find_set(index);
//          auto res = nearestNeighbors[index].front();
//          size_t sy = ds.find_set(res.index);
//          nearestNeighbors[index].pop_front();
//          if (sx != sy) {
//             return res;
//          }
//          // else just go on removing items
//       }


   if(shouldFind[index] && nearestNeighbors[index].empty())
   {
      std::priority_queue<HeapNeighborItem> heap;
      size_t clusterIndex = ds.find_set(index);

      double _tau = INFINITY;//maxRadiuses[index];

//       THIS IS SLOWER:
//          double _tau = (*_distance)(index,
//             *(ds.getClusterMembers(ds.getClusterNext(clusterIndex)).cbegin())
//          );

//       THIS IS SLOWER TOO:
//          size_t test = (size_t)(index+unif_rand()*(_n-index));
//          if (ds.find_set(test) != clusterIndex)
//             _tau = (*_distance)(index, test);

#ifdef GENERATE_STATS
      ++stats.nnCals;
#endif
      getNearestNeighborsFromMinRadiusRecursive( _root, index, clusterIndex, minRadiuses[index], _tau, heap );
      while( !heap.empty() ) {
         nearestNeighbors[index].push_front(heap.top());
         heap.pop();
      }
      // maxRadiuses[index] = INFINITY;
      size_t newNeighborsCount = nearestNeighbors[index].size();

      neighborsCount[index] += newNeighborsCount;
      if(neighborsCount[index] > _n - index || newNeighborsCount == 0)
         shouldFind[index] = false;

      if(newNeighborsCount > 0)
         minRadiuses[index] = nearestNeighbors[index].back().dist;
   }

   if(!nearestNeighbors[index].empty())
   {
#ifdef GENERATE_STATS
      ++stats.nnCount;
#endif
      auto res = nearestNeighbors[index].front();
      nearestNeighbors[index].pop_front();
      return res;
   }
   else
   {
      return HeapNeighborItem(SIZE_MAX,-INFINITY);
      //stop("nie ma sasiadow!");
   }
}


NumericMatrix HClustBiVpTreeSingle::compute()
{
   NumericMatrix ret(_n-1, 2);
   priority_queue<HeapHierarchicalItem> pq;

   // INIT: Pre-fetch a few nearest neighbors for each point
#if VERBOSE > 5
   Rprintf("[%010.3f] prefetching NNs\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   for (size_t i=0;i<_n;i++)
   {
#if VERBOSE > 7
      if (i % 1024 == 0) Rprintf("\r             prefetch NN: %d/%d", i, _n-1);
#endif
      Rcpp::checkUserInterrupt(); // may throw an exception, fast op
      HeapNeighborItem hi=getNearestNeighbor(i);

      if (hi.index != SIZE_MAX)
      {
         pq.push(HeapHierarchicalItem(i, hi.index, hi.dist));
      }
   }
#if VERBOSE > 7
   Rprintf("\r             prefetch NN: %d/%d\n", _n-1, _n-1);
#endif
#if VERBOSE > 5
   Rprintf("[%010.3f] merging clusters\n", clock()/(float)CLOCKS_PER_SEC);
#endif

#ifdef MB_IMPROVEMENT
   int nsqrt = (int)sqrt((double)_n);
#endif  // MB_IMPROVEMENT

   size_t i = 0;
   while(i < _n - 1)
   {
      //Rcout << "iteracja " << i << endl;
      //Rcout << "pq size = " << pq.size()<< endl;
      HeapHierarchicalItem hhi = pq.top();
      pq.pop();

      size_t s1 = ds.find_set(hhi.index1);
      size_t s2 = ds.find_set(hhi.index2);
      if (s1 != s2)
      {
         Rcpp::checkUserInterrupt(); // may throw an exception, fast op

         ret(i,0)=(double)hhi.index1;
         ret(i,1)=(double)hhi.index2;
         ++i;
         ds.link(s1, s2);
#ifdef MB_IMPROVEMENT
         size_t s_new = ds.find_set(s1);
         size_t s_old;
         //Rcout << "przed usuwaniem z clusters" << endl;
         if(s1==s_new)
         {
         	clusters.erase(s2);
         	s_old = s2;
         }
         else
         {
         	clusters.erase(s1);
         	s_old = s1;
         }
         //Rcout << "clusters size = " << clusters.size() << endl;
         if(i >= _n - nsqrt)
         {
            //Rcout << "po sqrt" << endl;
            mbimprovement = true;
            //Rcout << "aktualizuje clusters dist" << endl;
            for ( auto it = clusters.begin(); it != clusters.end(); ++it )
            {
               if(*it == s_new || *it == s_old) continue;
               SortedPoint spold = SortedPoint(*it, s_old);
               SortedPoint spnew = SortedPoint(*it, s_new);
               auto dold = distClust.find(spold);
               auto dnew = distClust.find(spnew);
               if(dold != distClust.end() && dnew != distClust.end())
               {
                  //Rcout << "znaleziono obie" << endl;
                  dnew->second = min(dnew->second, dold->second);
                  //Rcout << "po aktu" << endl;
               }
               else if(dold != distClust.end())
               {
                  //Rcout << "znaleziono tylko stara" << endl;
                  distClust.emplace(spnew, dold->second);
                  //Rcout << "po aktu2" << endl;
               }

               if(dold != distClust.end())
               {
                  //Rcout << "znaleziono stara" << endl;
                  distClust.erase(spold);
                  //Rcout << "po aktu3" << endl;
               }
            }
            //Rcout << "po aktualizacji clusters dist" << endl;
         }
#endif // MB_IMPROVEMENT
      }
#if VERBOSE > 7
      if (i % 1024 == 0) Rprintf("\r             %d / %d", i+1, _n);
#endif

      // ASSERT: hhi.index1 < hhi.index2
      HeapNeighborItem hi=getNearestNeighbor(hhi.index1);
      if(hi.index != SIZE_MAX)
         pq.push(HeapHierarchicalItem(hhi.index1, hi.index, hi.dist));
   }
#if VERBOSE > 7
   Rprintf("\r             %d / %d\n", _n, _n);
#endif
#if VERBOSE > 5
   Rprintf("[%010.3f] generating output matrix\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   Rcpp::checkUserInterrupt();

   MergeMatrixGenerator mmg(ret.nrow());
   return mmg.generateMergeMatrix(ret);
}


void HClustBiVpTreeSingle::print(HClustBiVpTreeNode* n) {
   if (n->ll) {
      Rprintf("\"%llx\" -> \"%llx\" [label=\"LL\"];\n", (unsigned long long)n, (unsigned long long)(n->ll));
      print(n->ll);
   }
   if (n->lr) {
      Rprintf("\"%llx\" -> \"%llx\" [label=\"LR\"];\n", (unsigned long long)n, (unsigned long long)(n->lr));
      print(n->lr);
   }
   if (n->rl) {
      Rprintf("\"%llx\" -> \"%llx\" [label=\"RL\"];\n", (unsigned long long)n, (unsigned long long)(n->rl));
      print(n->rl);
   }
   if (n->rr) {
      Rprintf("\"%llx\" -> \"%llx\" [label=\"RR\"];\n", (unsigned long long)n, (unsigned long long)(n->rr));
      print(n->rr);
   }
   if (n->vpindex == SIZE_MAX) {
      for (size_t i=n->left; i<n->right; ++i)
         Rprintf("\"%llx\" -> \"%llu\" [arrowhead = diamond];\n", (unsigned long long)n, (unsigned long long)_indices[i]+1);
   }
   else {
      Rprintf("\"%llx\" [label=\"(%llu, %g)\"];\n", (unsigned long long)n, (unsigned long long)n->vpindex+1, n->radius);
   }
}


void HClustBiVpTreeSingle::print() {
   Rprintf("digraph vptree {\n");
   Rprintf("size=\"6,6\";\n");
   Rprintf("node [color=lightblue2, style=filled];");
   print(_root);
   Rprintf("}\n");
}


// [[Rcpp::export(".hclust2_single")]]
RObject hclust2_single(RObject distance, RObject objects, int maxNumberOfElementsInLeaves=2) {
#if VERBOSE > 5
   Rprintf("[%010.3f] starting timer\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   RObject result(R_NilValue);
   DataStructures::Distance* dist = DataStructures::Distance::createDistance(distance, objects);

   try {
      /* Rcpp::checkUserInterrupt(); may throw an exception */
      DataStructures::HClustBiVpTreeSingle hclust(dist, (int)maxNumberOfElementsInLeaves);
      RObject merge = hclust.compute();
      result = Rcpp::as<RObject>(List::create(
         _["merge"]  = merge,
         _["height"] = R_NilValue,
         _["order"]  = R_NilValue,
         _["labels"] = R_NilValue,
         _["call"]   = R_NilValue,
         _["method"] = "single",
         _["dist.method"] = R_NilValue
      ));
      result.attr("class") = "hclust";
      //hclust.print();
   }
   catch(...) {

   }

   if (dist) delete dist;
#if VERBOSE > 5
   Rprintf("[%010.3f] done\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   if (Rf_isNull(result)) stop("stopping on error or explicit user interrupt");
   return result;
}
