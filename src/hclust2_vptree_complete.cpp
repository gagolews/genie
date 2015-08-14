/* ************************************************************************* *
 *   This file is part of the `DataStructures` package.                      *
 *                                                                           *
 *   Copyright 2015 Maciej Bartoszuk, Anna Cena, Marek Gagolewski,           *
 *                                                                           *
 *   Parts of the code are taken from the 'CITAN' R package by M. Gagolewski *
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

#include "hclust2_vptree_complete.h"

using namespace Rcpp;
using namespace std;
using namespace boost;
using namespace DataStructures;

#ifdef VERBOSE
#undef VERBOSE
#define VERBOSE 0
#endif

#define UPDATEKK
#define CALCULATETIMEFORCLUSTERDISTANCE

int HClustVpTreeComplete::chooseNewVantagePoint(size_t left, size_t right)
{
#if VANTAGE_POINT_SELECT_SCHEME == 1
   // idea by A. Fu et al., "Dynamic VP-tree indexing for n-nearest neighbor
   //    search given pair-wise distances"
   size_t numCandidates = 5;
   size_t numTest = 12;

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


HClustVpTreeCompleteNode* HClustVpTreeComplete::buildFromPoints(size_t left, size_t right)
{
   if(right - left <= maxNumberOfElementsInLeaves)
   {
      for (size_t i=left; i<right; ++i) {
         size_t j = _indices[(i+1 < right)?(i+1):left];
         if (_indices[i] < j)
            maxRadiuses[ _indices[i] ] = (*_distance)(_indices[i], j);
      }

      return new HClustVpTreeCompleteNode(left, right);
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
   HClustVpTreeCompleteNode* node = new HClustVpTreeCompleteNode(vpi, (*_distance)(vpi, _indices[median]));


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



   /*

   size_t calculateHClustVpTreeCompleteNodeSize(HClustVpTreeCompleteNode* node)
   {
      return node->radiuses.size()*sizeof(double)
         + node->points.size()*sizeof(int)
         + node->children.size()*sizeof(HClustVpTreeCompleteNode*)
         + sizeof(HClustVpTreeCompleteNode);
   }

   size_t treeSize_rec(HClustVpTreeCompleteNode* node)
   {
      size_t size = calculateHClustVpTreeCompleteNodeSize(node);
      for(int i=0;i<node->childCount;i++)
      {
         size += treeSize_rec(node->children[i]);
      }
      return size;
   }

   int treeHeight_rec(HClustVpTreeCompleteNode* node)
   {
      int maxH = 0;
      for(int i=0;i<node->childCount;i++)
      {
         maxH = max(treeHeight_rec(node->children[i]), maxH);
      }
      return maxH+1;
   }
*/

void HClustVpTreeComplete::getNearestNeighborsFromMinRadiusRecursive( HClustVpTreeCompleteNode* node, size_t index,
   size_t clusterIndex, double minR, double& maxR,
   std::priority_queue<HeapNeighborItem>& heap )
{
   /*if(index == 3)
   {
      Rcout << "index 3 to klaster " << clusterIndex << endl;
      Rcout << "maxR = " << maxR << endl;
      Rcout << "minR = " << minR << endl;
   }*/
   // search within (minR, maxR]
   if (node == NULL) return;

   if (node->sameCluster) {
      if (node->vpindex == SIZE_MAX) {
         if (ds.find_set(_indices[node->left]) == clusterIndex) return;
      } else {
         if (ds.find_set(node->vpindex) == clusterIndex) return;
      }
   }

   if(node->vpindex == SIZE_MAX) // leaf
   {
      if (node->sameCluster) {
#ifdef MB_IMPROVEMENT
      size_t s = SIZE_MAX;
      if(mbimprovement)
         s = ds.find_set(_indices[node->left]);
      std::unordered_map<SortedPoint,double>::const_iterator distToClusterIterator;
      if(mbimprovement)
         distToClusterIterator = distClust.find(SortedPoint(s, clusterIndex));
      double distToCluster = INFINITY;
      if(mbimprovement && distToClusterIterator != distClust.end())
         distToCluster = distToClusterIterator->second;
#else
      size_t s = ds.find_set(_indices[node->left]);
#endif
      /*
      for(size_t i=node->left; i<node->right; i++)
      {
         if(index == 3)
         {
            Rcout << "bedac 4, moglbym napotkac "<<  _indices[i] + 1 << "z klastra "<< ds.find_set(_indices[i]) + 1 << endl;
         }
      }*/
      //if(KK.find(SortedPoint(clusterIndex, s)) != KK.end()) return; //very important line
      auto s1s2dist = KK.find(SortedPoint(clusterIndex, s));
      if(s1s2dist != KK.end())
      {
         /*if(index == 3)
         {
            Rcout << "index 3, a to klaster " << s<< " ktory odrzucam, bo znaleziona odleglosc to " << s1s2dist->second.dist << endl;
            Rcout << ds.find_set(8) << endl;
         }*/
         if(timestamp.find(clusterIndex)->second < s1s2dist->second.iter
               && timestamp.find(s)->second < s1s2dist->second.iter)
         return; //very important line
      }

         for(size_t i=node->left; i<node->right; i++)
         {
            /*if(index == 3)
            {
               Rcout << "bedac 4, napotkalem na "<<  _indices[i] + 1 << endl;
            }*/
            if(index >= _indices[i]) continue;
            double dist2 = (*_distance)(index, _indices[i]);
            if (dist2 > maxR || dist2 <= minR) continue;

#ifdef MB_IMPROVEMENT
            if(mbimprovement && dist2 > distToCluster) {
               //Rcout << "odrzucam!" << endl;
               continue;
            }
#endif
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
            if(mbimprovement)
               distClust.emplace(SortedPoint(s, clusterIndex), dist2);
#endif
         }
      }
      else {
         size_t commonCluster = ds.find_set(_indices[node->left]);
         for(size_t i=node->left; i<node->right; i++)
         {
            /*if(index == 3)
            {
               Rcout << "bedac 4, napotkalem na "<<  _indices[i] + 1 << endl;
            }*/
            size_t currentCluster = ds.find_set(_indices[i]);
#ifdef MB_IMPROVEMENT
            std::unordered_map<SortedPoint,double>::const_iterator distToClusterIterator;
            if(mbimprovement)
               distToClusterIterator = distClust.find(SortedPoint(currentCluster, clusterIndex));
            double distToCluster = INFINITY;
            if(mbimprovement && distToClusterIterator != distClust.end())
               distToCluster = distToClusterIterator->second;
#endif


            if (currentCluster != commonCluster) commonCluster = SIZE_MAX;
            if (currentCluster == clusterIndex) continue;

            auto s1s2dist = KK.find(SortedPoint(clusterIndex, currentCluster));

            if(s1s2dist != KK.end())
               {
               /*if(index == 3)
               {
                  Rcout << "index 3, a to klaster " << currentCluster<< "ktory odrzucam, bo znaleziona odleglosc to " << s1s2dist->second.dist << endl;
               }*/
               if(timestamp.find(clusterIndex)->second < s1s2dist->second.iter
                     && timestamp.find(currentCluster)->second < s1s2dist->second.iter)
               continue; //very important line
               }

            if (index >= _indices[i]) continue;

            double dist2 = (*_distance)(index, _indices[i]);
            if (dist2 > maxR || dist2 <= minR) continue;
#ifdef MB_IMPROVEMENT
            if(mbimprovement && dist2 > distToCluster) continue;
#endif

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
            if(mbimprovement)
               distClust.emplace(SortedPoint(currentCluster, clusterIndex), dist2);
#endif
         }
         if (commonCluster != SIZE_MAX)
            {
               /*for(size_t i=node->left; i<node->right; i++)
               {
                  Rcout << "ustawiam sameCluster na true, bo "<<  _indices[i] + 1 << "z klastra "<< ds.find_set(_indices[i]) + 1 << endl;
               }*/
               node->sameCluster = true;
            }
      }
      return;
   }
   // else // not a leaf
   double dist = (*_distance)(node->vpindex, index);
   /*if (dist < maxR && dist > minR && index < node->vpindex && ds.find_set(node->vpindex) != clusterIndex) {

      heap.push( HeapNeighborItem(node->vpindex, dist) );
      maxR = heap.top().dist;
   }*/
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

void HClustVpTreeComplete::print(HClustVpTreeCompleteNode* n) {
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

HeapNeighborItem HClustVpTreeComplete::getNearestNeighbor(size_t index)
{
#if VERBOSE > 5
   // Rprintf(".");
#endif
   if(shouldFind[index] && nearestNeighbors[index].empty())
   {
      std::priority_queue<HeapNeighborItem> heap;
      size_t clusterIndex = ds.find_set(index);

      double _tau = maxRadiuses[index];
//       THIS IS SLOWER:
//          size_t test = (size_t)(index+unif_rand()*(_n-index));
//          if (ds.find_set(test) != clusterIndex)
//             _tau = (*_distance)(index, test);

      getNearestNeighborsFromMinRadiusRecursive( _root, index, clusterIndex, minRadiuses[index], _tau, heap );
      while( !heap.empty() ) {
         nearestNeighbors[index].push_front(heap.top());
         heap.pop();
      }
      maxRadiuses[index] = INFINITY;
      size_t newNeighborsCount = nearestNeighbors[index].size();

      neighborsCount[index] += newNeighborsCount;
      if(neighborsCount[index] > _n - index || newNeighborsCount == 0)
         shouldFind[index] = false;

      if(newNeighborsCount > 0)
         minRadiuses[index] = nearestNeighbors[index].back().dist;
   }

   if(!nearestNeighbors[index].empty())
   {
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

#ifdef DEBUG
void printCounters()
{
   _distance.printCounters();
}
#endif

// constructor (OK, we all know what this is, but I label it for faster in-code search)
HClustVpTreeComplete::HClustVpTreeComplete(Distance* dist, size_t maxNumberOfElementsInLeaves) :
   maxNumberOfElementsInLeaves(maxNumberOfElementsInLeaves),
   _root(NULL), _n(dist->getObjectCount()), _distance(dist),
   _indices(dist->getObjectCount()),
   neighborsCount(vector<size_t>(dist->getObjectCount(), 0)),
   minRadiuses(vector<double>(dist->getObjectCount(), -INFINITY)),
   maxRadiuses(vector<double>(dist->getObjectCount(), INFINITY)),
   shouldFind(vector<bool>(dist->getObjectCount(), true)),
   nearestNeighbors(vector< deque<HeapNeighborItem> >(dist->getObjectCount())),
#ifdef USE_BOOST_DISJOINT_SETS
   ds(make_assoc_property_map(rank), make_assoc_property_map(parent))
#else
   ds(dist->getObjectCount())
#endif
{
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
     ;//ds.make_set(i);
#ifdef MB_IMPROVEMENT
     clusters.emplace(i);
#endif
   }

   _root = buildFromPoints(0, _n);
}

HClustVpTreeComplete::~HClustVpTreeComplete() {
#if VERBOSE > 5
   Rprintf("[%010.3f] destroying vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   if(_root) delete _root;
}


   /*size_t treeSize()
   {
      if(_root==NULL) return sizeof(VpTree);
      return sizeof(VpTree) + treeSize_rec(_root);
   }

   int treeHeight()
   {
      if(_root==NULL) return 0;
      return treeHeight_rec(_root);
   }*/

void HClustVpTreeComplete::print() {
   Rprintf("digraph vptree {\n");
   Rprintf("size=\"6,6\";\n");
   Rprintf("node [color=lightblue2, style=filled];");
   print(_root);
   Rprintf("}\n");
}

HClustVpTreeComplete::HeapHierarchicalItemMax HClustVpTreeComplete::calculateCluster2ClusterMaxDistance(size_t item1, size_t item2, size_t iter)
{
   size_t s1 = ds.find_set(item1);
   size_t s2 = ds.find_set(item2);
   HeapHierarchicalItemMax ret(s1,s2, -1, iter);

   for(size_t i=0;i<_n;++i)
   {
      size_t clusterRepresentant1 = ds.find_set(i);
      if(clusterRepresentant1 != s1)
         continue;
      for(size_t j=0;j<_n;++j)
      {
         size_t clusterRepresentant2 = ds.find_set(j);
         if(clusterRepresentant2 != s2)
            continue;
         double dist = (*_distance)(i, j);
         if(ret.dist < dist)
            ret.dist = dist;
      }
   }
   return ret;
}

NumericMatrix HClustVpTreeComplete::compute()
{
   Rcout << "wchodze do complete!" << endl;
   static double timeForClusterDistance = 0;
   static double timeForKKUpdate = 0;
   static double timeForMinSearch1 = 0;
   static double timeForMinSearch2 = 0;
   static double timeForMinSearchPrefetch = 0;
   static double timeForPriorityQueues = 0;
   NumericMatrix ret(_n-1, 2);
   priority_queue<HeapHierarchicalItem> pqMin;
   priority_queue<HeapHierarchicalItemMax> pqMax;

   //vector<bool> stillExists(_n, true);

   for(size_t i=0;i<_n;i++)
   {
      timestamp.emplace(i, 0);
   }

   // INIT: Pre-fetch a few nearest neighbors for each point
#if VERBOSE > 5
   Rprintf("[%010.3f] prefetching NNs\n", clock()/(float)CLOCKS_PER_SEC);
#endif
#if VERBOSE > 3
   int misses = 0;
#endif
   for(size_t i=0;i<_n;i++)
   {
      if (true) Rcpp::checkUserInterrupt(); // may throw an exception
#if VERBOSE > 7
      Rprintf("\r             prefetch NN: %d/%d", i, _n-1);
#endif
      double before = clock()/(double)CLOCKS_PER_SEC;
      HeapNeighborItem hi=getNearestNeighbor(i);
      double after = clock()/(double)CLOCKS_PER_SEC;
      timeForMinSearchPrefetch += after - before;
      if(hi.index != SIZE_MAX)
      {
         //Rcout <<"dla " << i << "najblizszym jest " << hi->index << endl;
         pqMin.push(HeapHierarchicalItem(i, hi.index, hi.dist));
      }
   }
#if VERBOSE > 7
   Rprintf("\n");
#endif
#if VERBOSE > 5
   Rprintf("[%010.3f] merging clusters\n", clock()/(float)CLOCKS_PER_SEC);
#endif

   size_t i = 0;
   size_t iter = 0;
   pqMax.push(HeapHierarchicalItemMax(SIZE_MAX, SIZE_MAX, INFINITY, SIZE_MAX));
   pqMin.push(HeapHierarchicalItem(SIZE_MAX, SIZE_MAX, INFINITY));
   bool awaria = false;
   while(i < _n - 1)
   {
      iter++;
      //if(iter > 3000) {awaria = true; break;}
#if VERBOSE > 10
      Rcout << "i " << i << endl;
      Rcout << "iteracja " << iter << endl;
      Rcout << "pqMin size = " << pqMin.size()<< endl;
      Rcout << "pqMax size = " << pqMax.size()<< endl;
#endif
      double beforeQueue = clock()/(double)CLOCKS_PER_SEC;
      HeapHierarchicalItem hhiMin = pqMin.top();
      HeapHierarchicalItemMax hhiMax = pqMax.top();
#if VERBOSE > 11
         Rcout << hhiMin.index1+1 << " " <<hhiMin.index2+1 <<" " << hhiMin.dist << " w pqMin" << endl;
#endif


#ifdef UPDATEKK
      if(hhiMax.iter != SIZE_MAX)
      {
         size_t s1 = ds.find_set(hhiMax.index1);
         size_t s2 = ds.find_set(hhiMax.index2);
         hhiMax.index1 = s1;
         hhiMax.index2 = s2;
         if(s1 == s2) {pqMax.pop(); continue;}
         auto pair1 = KK.find(SortedPoint(s1,s2));
         if(pair1 != KK.end())
         {

            if(pair1->second.dist > hhiMax.dist)
            {
#if VERBOSE > 11
               Rcout << "uaktualniam odleglosc w kolejce pqMax" << endl;
               Rcout << "BYLO:" <<hhiMax.dist << " " << hhiMax.iter << endl;
               Rcout << "JEST:" <<pair1->second.dist << " " << pair1->second.iter << endl;
#endif
               hhiMax.dist = pair1->second.dist;
               hhiMax.iter = pair1->second.iter;
               HeapHierarchicalItemMax tmp = HeapHierarchicalItemMax(hhiMax);
               pqMax.pop();
               pqMax.push(tmp);
#if VERBOSE > 11
               Rcout << "po uaktualnieniu odleglosci w kolejce pqMax" << endl;
#endif
               continue;
            }
         }
      }
#endif

      bool upToDate = false;
      if((hhiMax.iter == SIZE_MAX)
            || (hhiMax.iter> timestamp.find(hhiMax.index1)->second
                  && hhiMax.iter> timestamp.find(hhiMax.index2)->second))
      {
#if VERBOSE > 11
         Rcout << hhiMax.index1+1 << " " <<hhiMax.index2+1 <<" " << hhiMax.dist << " nadal aktualne" << endl;
#endif
         upToDate = true;
      }

      if(!upToDate)
      {
#ifndef UPDATEKK
         if(ds.find_set(hhiMax.index1) == ds.find_set(hhiMax.index2)) {pqMax.pop(); continue;}
#endif
#if VERBOSE > 11
         Rcout << hhiMax.index1+1 << " " <<hhiMax.index2+1 <<" " << hhiMax.dist << " nieaktualne" << endl;
#endif
#ifdef CALCULATETIMEFORCLUSTERDISTANCE
         double  before = clock()/(double)CLOCKS_PER_SEC;
#endif
         HeapHierarchicalItemMax hhim = calculateCluster2ClusterMaxDistance(hhiMax.index1, hhiMax.index2, iter);
#ifdef CALCULATETIMEFORCLUSTERDISTANCE
         double after = clock()/(double)CLOCKS_PER_SEC;
         timeForClusterDistance += after - before;
#endif
#if VERBOSE > 11
         Rcout << "najwieksza odleglosc klastrowa to" << hhim.dist << endl;
         Rcout << "wrzucam " << hhim.index1+1 << ", " << hhim.index2+1 << endl;
#endif
         KK[SortedPoint(hhiMax.index1, hhiMax.index2)] = KKItem(hhim.dist, hhim.iter);
         pqMax.pop();
         pqMax.push(hhim);
         continue;
      }
      double afterQueue = clock()/(double)CLOCKS_PER_SEC;
      timeForPriorityQueues = afterQueue - beforeQueue;
      if(hhiMin.dist < hhiMax.dist) //przetwarzamy element z PQmin
      {
         pqMin.pop();
         size_t s1 = ds.find_set(hhiMin.index1);
         size_t s2 = ds.find_set(hhiMin.index2);
#if VERBOSE > 11
         Rcout << "przetwarzam z PQmin, " << hhiMin.index1+1 << " " << hhiMin.index2+1 << endl;
         Rcout << "ich reprezentanci to " << s1+1 << " " << s2+1 << endl;
#endif
         if(s1!=s2) //czy aby na pewno konieczne?
         {
            if(ds.getClusterSize(ds.find_set(hhiMin.index1)) == 1 && ds.getClusterSize(ds.find_set(hhiMin.index2)) == 1)
            {
#if VERBOSE > 11
               Rcout << "1 elementowe zbiory to sa" << endl;
#endif
               ret(i,0)=(double)hhiMin.index1;
               ret(i,1)=(double)hhiMin.index2;
               ++i;
               ds.link(s1, s2);
               // size_t s = ds.find_set(hhiMin.index1);
               timestamp[hhiMin.index1] = iter;
               timestamp[hhiMin.index2] = iter;
               /*if(hhiMin.index1 != s)
                  stillExists[hhiMin.index1] = false;
               else
                  stillExists[hhiMin.index2] = false;*/
            }
            else
            {
#if VERBOSE > 11
               Rcout << "co najmniej 1 klaster" << endl;
#endif
               auto s1s2dist = KK.find(SortedPoint(s1,s2));
               if(s1s2dist != KK.end())
               {
                  if(timestamp.find(s1)->second < s1s2dist->second.iter
                        && timestamp.find(s2)->second < s1s2dist->second.iter)
                  {
#if VERBOSE > 11
                     Rcout << "odleglosc aktualna" << endl;
#endif
                     double before = clock()/(double)CLOCKS_PER_SEC;
                     HeapNeighborItem hi=getNearestNeighbor(hhiMin.index1);
                     double after = clock()/(double)CLOCKS_PER_SEC;
#if VERBOSE > 11
         Rcout << "Dla " << hhiMin.index1+1 << " znaleziono " <<hi.index+1 <<" jako najblizszego sasiada (dist = " << hi.dist << " )" << endl;
#endif

                     timeForMinSearch1 += after - before;
                     if(hi.index != SIZE_MAX)
                        pqMin.push(HeapHierarchicalItem(hhiMin.index1, hi.index, hi.dist));
                     continue;
                  }
#if VERBOSE > 11
                  Rcout << "odleglosc nieaktualna" << endl;
#endif
               }
#if VERBOSE > 11
               Rcout << "odleglosci nie ma albo nieaktualna" << endl;
#endif
#ifdef CALCULATETIMEFORCLUSTERDISTANCE
               double before = clock()/(double)CLOCKS_PER_SEC;
#endif
               HeapHierarchicalItemMax hhim = calculateCluster2ClusterMaxDistance(s1,s2,iter);
#ifdef CALCULATETIMEFORCLUSTERDISTANCE
               double after = clock()/(double)CLOCKS_PER_SEC;
               timeForClusterDistance += after - before;
#endif
#if VERBOSE > 11
               Rcout << "najwieksza odleglosc klastrowa to" << hhim.dist << endl;
               Rcout << "wrzucam " << hhim.index1+1 << ", " << hhim.index2+1 << endl;
#endif
               KK[SortedPoint(hhim.index1, hhim.index2)] = KKItem(hhim.dist, hhim.iter);
               pqMax.push(hhim);

            }
         }
         double before = clock()/(double)CLOCKS_PER_SEC;
         HeapNeighborItem hi=getNearestNeighbor(hhiMin.index1);
         double after = clock()/(double)CLOCKS_PER_SEC;
#if VERBOSE > 11
         Rcout << "Dla " << hhiMin.index1+1 << " znaleziono " <<hi.index+1 <<" jako najblizszego sasiada (dist = " << hi.dist << " )" << endl;
#endif
         timeForMinSearch2 += after - before;
         if(hi.index != SIZE_MAX)
            pqMin.push(HeapHierarchicalItem(hhiMin.index1, hi.index, hi.dist));
      }
      else // przetwarzamy element z PQmax
      {
#if VERBOSE > 11
         Rcout << "przetwarzam z PQMax" << endl;
#endif
         pqMax.pop();
         ds.union_set(hhiMax.index1, hhiMax.index2); //hhiMax.index1 is always a parent
         ret(i,0)=(double)hhiMax.index1;
         ret(i,1)=(double)hhiMax.index2;
         ++i;
         timestamp[hhiMax.index1] = iter;
         timestamp[hhiMax.index2] = iter;
#if VERBOSE > 11
         Rcout << "reprezentant teraz to " << hhiMax.index1+1 << endl;
#endif
#ifdef UPDATEKK
         //idz po wszystykich entry w KK, aktualizuj i dawaj timestamp nowy!

#if VERBOSE > 10
               Rcout << "KK size = " << KK.size() << endl;
#endif
         double before = clock()/(double)CLOCKS_PER_SEC;

         for ( auto it = KK.begin(); it != KK.end(); ++it )
         {
            size_t otherCluster;
            if(it->first.i == hhiMax.index1)
            {
               otherCluster = it->first.j;
            }
            else if(it->first.j == hhiMax.index1)
            {
               otherCluster = it->first.i;
            }
            else
            {
               continue;
            }
            auto pair2 = KK.find(SortedPoint(otherCluster, hhiMax.index2));
            if(pair2 != KK.end())
            {
               KK[SortedPoint(otherCluster, hhiMax.index1)] = KKItem(max(it->second.dist, pair2->second.dist), iter);
               KK.erase(SortedPoint(otherCluster, hhiMax.index2));
            }
         }
         double after = clock()/(double)CLOCKS_PER_SEC;
         timeForKKUpdate += after - before;
#endif
      }

         Rcpp::checkUserInterrupt(); // may throw an exception
      }
#if VERBOSE > 7
      if (i % 1024 == 0) Rprintf("\r             %d / %d / %d ", i+1, _n, misses);
#endif

#if VERBOSE > 3
   Rprintf("Total ignored NNs: %d\n", misses);
#endif
#if VERBOSE > 5
Rprintf("[%010.3f] generating output matrix\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   Rcpp::checkUserInterrupt();

#ifdef CALCULATETIMEFORCLUSTERDISTANCE
   Rcout << "time for calculating cluster distances = " << timeForClusterDistance << endl;
   Rcout << "time for updating KK = " << timeForKKUpdate << endl;
   Rcout << "timeForMinSearch1 = " << timeForMinSearch1 << endl;
   Rcout << "timeForMinSearch2 = " << timeForMinSearch2 << endl;
   Rcout << "timeForMinSearchPrefetch = " << timeForMinSearchPrefetch << endl;
   Rcout << "timeForPriorityQueues = " << timeForPriorityQueues << endl;
#endif

   MergeMatrixGenerator mmg(ret.nrow());
   if(!awaria)
   {
#if VERBOSE > 11
      Rcout << "wyszedlem poprawnie!" << endl;
#endif
      return mmg.generateMergeMatrix(ret);
   }
   else
   {
#if VERBOSE > 11
      Rcout << "byla awaria" << endl;
#endif
      return ret;
   }
}


// [[Rcpp::export(".hclust2_complete")]]
RObject hclust2_complete(RObject distance, RObject objects, int maxNumberOfElementsInLeaves=2) {
#if VERBOSE > 5
   Rprintf("[%010.3f] starting timer\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   RObject result(R_NilValue);
   DataStructures::Distance* dist = DataStructures::Distance::createDistance(distance, objects);

   try {
      /* Rcpp::checkUserInterrupt(); may throw an exception */
      DataStructures::HClustVpTreeComplete hclust(dist, (int)maxNumberOfElementsInLeaves);
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
