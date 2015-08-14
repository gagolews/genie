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

#include "hclust2_vptree_medoid.h"
//#define VERBOSE 17
/* Improvement ideas:
 *
 * 1. add custom sort of input objects
 *       useful for Levenshtein distance
 *       long strings should be put at the end
 */


using namespace Rcpp;
using namespace std;
using namespace boost;
using namespace DataStructures;


// constructor (OK, we all know what this is, but I label it for faster in-code search)
HClustBiVpTreeMedoid::HClustBiVpTreeMedoid(Distance* dist, RObject control) :
   opts(control), _root(NULL), _n(dist->getObjectCount()), _distance(dist),
   _indices(dist->getObjectCount()),
   //_indicesinv(dist->getObjectCount()),
   neighborsCount(vector<size_t>(dist->getObjectCount(), 0)),
   minRadiuses(vector<double>(dist->getObjectCount(), -INFINITY)),
   // maxRadiuses(vector<double>(dist->getObjectCount(), INFINITY)),
   shouldFind(vector<bool>(dist->getObjectCount(), true)),
   nearestNeighbors(vector< deque<HeapNeighborItem> >(dist->getObjectCount())),
   distances(vector<double>(_n)),
   //medoids(dist->getObjectCount()),
   medoidFound(dist->getObjectCount(), 0),
   nnback(dist->getObjectCount()),
#ifdef GENERATE_STATS
   stats(HClustTreeStats()),
#endif
   ds(dist->getObjectCount())
{
#if VERBOSE > 5
   Rprintf("[%010.3f] building vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
#endif

   // starting _indices: random permutation of {0,1,...,_n-1}
   for (size_t i=0;i<_n;i++)
      _indices[i] = i;
   for (size_t i=_n-1; i>= 1; i--)
      swap(_indices[i], _indices[(size_t)(unif_rand()*(i+1))]);

   _root = buildFromPoints(0, _n);
}


HClustBiVpTreeMedoid::~HClustBiVpTreeMedoid() {
// #if VERBOSE > 5
//       Rprintf("[%010.3f] destroying vp-tree\n", clock()/(float)CLOCKS_PER_SEC);
// #endif
   if(_root) delete _root;
}


size_t HClustBiVpTreeMedoid::chooseNewVantagePoint(size_t left, size_t right)
{
   if (opts.vpSelectScheme == 1) {
      // idea by Yianilos (original vp-tree paper)
      if (left + opts.vpSelectCand + opts.vpSelectTest > right)
         return left;

      // randomize:
      for (size_t i=left; i<left+opts.vpSelectCand+opts.vpSelectTest; ++i)
         std::swap(_indices[i], _indices[i+(size_t)(unif_rand()*(right-i))]);

      // maximize variance
      size_t bestIndex = -1;
      double bestSigma = -INFINITY;
      for (size_t i=left; i<left+opts.vpSelectCand; i++) {
         accumulators::accumulator_set< double,
            accumulators::features<accumulators::tag::variance> > acc;
         for (size_t j = left+opts.vpSelectCand; j < left+opts.vpSelectCand+opts.vpSelectTest; ++j)
            acc( (*_distance)( _indices[i], _indices[j] ) );
         double curSigma = accumulators::variance(acc);
         if (curSigma > bestSigma) {
            bestSigma = curSigma;
            bestIndex = i;
         }
      }

      return bestIndex;
   }
   else if (opts.vpSelectScheme == 2) {
      // idea by T. Bozkaya and M. Ozsoyoglu, "Indexing large metric spaces
      //      for similarity search queries"

      // randomize:
      std::swap(_indices[left], _indices[left+(size_t)(unif_rand()*(right-left))]);

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
   }
   else {
      // return random index
      // don'use left one (even if sample seems to be randomized already,
      // vp in subtrees is already on the left...)
      return left+(size_t)(unif_rand()*(right-left));
   }
}

//bool comparer_gt(int i,int j) { return (i>j); }


HClustBiVpTreeMedoidNode* HClustBiVpTreeMedoid::buildFromPoints(size_t left,
   size_t right)
{
#ifdef GENERATE_STATS
   ++stats.nodeCount;
#endif
   if (right - left <= opts.maxLeavesElems)
   {
   #ifdef GENERATE_STATS
      ++stats.leafCount;
   #endif
      HClustBiVpTreeMedoidNode* leaf = new HClustBiVpTreeMedoidNode(left, right);
      // std::sort(_indices.begin()+left, _indices.begin()+right, comparer_gt);
      // leaf->maxindex = _indices[left];
      // leaf->maxindex = _indices[left];
      // for (size_t i=left+1; i<right; ++i)
         // if (_indices[i] > leaf->maxindex)
            // leaf->maxindex = _indices[i];
 //     for (size_t i=left; i<right; ++i)
 //        _indicesinv[_indices[i]] = i;
      leaf->maxindex = right-1;
      return leaf;
   }

   size_t vpi_idx = chooseNewVantagePoint(left, right);
   std::swap(_indices[left], _indices[vpi_idx]);
   size_t vpi = _indices[left];
   //_indicesinv[vpi] = left;
   size_t median = (right + left) / 2;

   for (size_t i=left+1; i<right; ++i)
      distances[_indices[i]] = (*_distance)(vpi, _indices[i]);
   //std::sort(_indices.begin()+left+1, _indices.begin()+right, DistanceComparatorCached(&distances));
   std::nth_element(_indices.begin()+left+1, _indices.begin() + median, _indices.begin()+right, DistanceComparatorCached(&distances));

// slower -- computes some distances > 1 time
//    std::nth_element(_indices.begin() + left + 1, _indices.begin() + median,  _indices.begin() + right,
//                     DistanceComparator(vpi, _distance));
//    HClustBiVpTreeMedoidNode* node = new HClustBiVpTreeMedoidNode(vpi, left, left+1, (*_distance)(vpi, _indices[median]));

   /*if (left == 0) { // for the root node we get NNs for free
      for (size_t i=left+1; i<right; ++i) {
         nearestNeighbors[vpi].push_back(HeapNeighborItem(_indices[i], distances[_indices[i]]));
         ++neighborsCount[vpi];
      }
      if (neighborsCount[vpi] > _n - left) shouldFind[vpi] = false;
      if (neighborsCount[vpi] > 0) minRadiuses[vpi] = nearestNeighbors[vpi].back().dist;
   }*/
// doesn't work: some of NNs can be determined here
//    if (isLeftChild) {
//       bool distanceToParentVP = ((parentVP == SIZE_MAX)?0.0:(*_distance)(parentVP, vpi));
//       for (size_t i=left+1; i<right; ++i) {
//          if (distanceToParentVP+distances[_indices[i]] >= parentRadius) break;
//          nearestNeighbors[vpi].push_back(HeapNeighborItem(_indices[i], distances[_indices[i]]));
//          ++neighborsCount[vpi];
//       }
//       if (neighborsCount[vpi] > _n - left) shouldFind[vpi] = false;
//       if (neighborsCount[vpi] > 0) minRadiuses[vpi] = nearestNeighbors[vpi].back().dist;
//    }
//    else {
//       bool distanceToParentVP = ((parentVP == SIZE_MAX)?0.0:(*_distance)(parentVP, vpi));
//       for (size_t i=left+1; i<right; ++i) {
//          if (distances[_indices[i]] >= distanceToParentVP-parentRadius) break;
//          nearestNeighbors[vpi].push_back(HeapNeighborItem(_indices[i], distances[_indices[i]]));
//          ++neighborsCount[vpi];
//       }
//       if (neighborsCount[vpi] > _n - left) shouldFind[vpi] = false;
//       if (neighborsCount[vpi] > 0) minRadiuses[vpi] = nearestNeighbors[vpi].back().dist;
//    }

   HClustBiVpTreeMedoidNode* node = new HClustBiVpTreeMedoidNode(vpi, left, left+1, distances[_indices[median]]);

   node->maxindex = left;
   if (median - left > 0) { // don't include vpi
      node->childL = buildFromPoints(left+1, median+1);
      if (node->childL->maxindex > node->maxindex)
         node->maxindex = node->childL->maxindex;
   }
   if (right - median - 1 > 0) {
      node->childR = buildFromPoints(median+1, right);
      if (node->childR->maxindex > node->maxindex)
         node->maxindex = node->childR->maxindex;
   }

   return node;
}


void HClustBiVpTreeMedoid::getNearestNeighborsFromMinRadiusRecursive(
   HClustBiVpTreeMedoidNode* node, size_t index,
   size_t clusterIndex, double minR, double& maxR, NNHeap& nnheap)
{
   // search within (minR, maxR]
   // if (node == NULL) return; // this should not happen
#ifdef GENERATE_STATS
#ifdef _OPENMP
#pragma omp atomic
#endif
   ++stats.nodeVisit;
#endif

   if (!prefetch && node->sameCluster && clusterIndex == ds.find_set(node->left))
      return;

   if (node->vpindex == SIZE_MAX) { // leaf
      if (!prefetch && !node->sameCluster) {
         size_t commonCluster = ds.find_set(node->left);
         for (size_t i=node->left; i<node->right; ++i) {
            //RCOUT("Rozwazam punkt " << _indices[i], 3);
            size_t currentCluster = ds.find_set(i);
            if (currentCluster != commonCluster) commonCluster = SIZE_MAX;
            if (currentCluster == clusterIndex) continue;
            if (index == currentCluster) continue;
            if(currentCluster != i) continue;
            if(medoidFound[currentCluster]==timestamp) continue;
            //RCOUT("Punkt ten przeszedl przez warunki",3);
            double dist2 = (*_distance)(_indices[index], _indices[currentCluster]); // the slow part
            medoidFound[currentCluster] = timestamp;
            if (dist2 > maxR || dist2 <= minR) continue;
            //RCOUT("Punkt wrzucam do kolejki", 3);
            nnheap.insert(currentCluster, dist2, maxR);
         }
         if (commonCluster != SIZE_MAX)
            node->sameCluster = true; // set to true (btw, may be true already)
      }
      if(prefetch)
      {
         for (size_t i=node->left; i<node->right; ++i) 
         {
            //size_t currentCluster = ds.find_set(i);
            //size_t medoid = medoids[currentCluster];
            if (index == i) continue;
            //if(medoid != i) continue;
            //if(medoidFound[currentCluster]) continue;
            double dist2 = (*_distance)(_indices[index], _indices[i]); // the slow part
            if (dist2 > maxR || dist2 <= minR) continue;
            nnheap.insert(i, dist2, maxR);
         }
      }
      else /* node->sameCluster -- komentarze takie jak ten potrafia byc bardzo mylace! */ {
         //tutaj pojawia sie pytanie, czy moze tez sprawdzac, czy w tym node jest w ogole medoid. jak nie ma, to moze w ogole nie liczyc?
         size_t currentCluster = ds.find_set(node->left);
         if (index == currentCluster) return;
         if(medoidFound[currentCluster]==timestamp) return;
         double dist2 = (*_distance)(_indices[index], _indices[currentCluster]); // the slow part
         medoidFound[currentCluster] = timestamp;
         if (dist2 > maxR || dist2 <= minR) return;
         nnheap.insert(currentCluster, dist2, maxR);
      }
      return; // nothing more to do
   }
   // else // not a leaf

   // first visit the vantage point
   double dist;
   if(prefetch)
   {
      //size_t currentCluster = ds.find_set(node->left);  
      //size_t medoid = medoids[currentCluster];
      dist = (*_distance)(_indices[index], _indices[node->left]); // the slow part
      //RCOUT("Rozwazam punkt " << _indices[node->left] << ", jest on vp", 3);
      if (index < node->left && dist <= maxR && dist > minR) {
      //RCOUT("Wrzucam go", 3);
         nnheap.insert(node->left, dist, maxR);
      }
      //if(node->sameCluster) return;
   }
   else
   {
      size_t currentCluster = ds.find_set(node->left);  
      //size_t medoid = medoids[currentCluster];
      if(medoidFound[currentCluster] < timestamp && node->left == currentCluster)
      {
         double dist = (*_distance)(_indices[index], _indices[currentCluster]); // the slow part
         //RCOUT("Rozwazam punkt " << _indices[node->left] << ", jest on vp", 3);
         if (index != node->left && dist <= maxR && dist > minR &&
               currentCluster != clusterIndex) 
         {
            nnheap.insert(currentCluster, dist, maxR);
            medoidFound[currentCluster] = timestamp;
         }
      }
      if(node->sameCluster) return;
      dist = (*_distance)(_indices[index], _indices[node->left]); // the slow part
   }

   if (dist < node->radius) {
      if (dist - maxR <= node->radius && dist + node->radius > minR) {
         if (node->childL )//&& (!prefetch ||  index < node->childL->maxindex))
            getNearestNeighborsFromMinRadiusRecursive(node->childL, index, clusterIndex, minR, maxR, nnheap);
      }

      if (dist + maxR >= node->radius) {
         if (node->childR)// && (!prefetch ||  index < node->childR->maxindex))
            getNearestNeighborsFromMinRadiusRecursive(node->childR, index, clusterIndex, minR, maxR, nnheap);
      }
   }
   else /* ( dist >= node->radius ) */ {
      if (dist + maxR >= node->radius) {
         if (node->childR)//  && (!prefetch ||  index < node->childR->maxindex) )
            getNearestNeighborsFromMinRadiusRecursive(node->childR, index, clusterIndex, minR, maxR, nnheap);
      }

      if (dist - maxR <= node->radius && dist + node->radius > minR) {
         if (node->childL)// && (!prefetch ||  index < node->childL->maxindex))
            getNearestNeighborsFromMinRadiusRecursive(node->childL, index, clusterIndex, minR, maxR, nnheap);
      }
   }
   if (prefetch || node->sameCluster ||
      (node->childL && !node->childL->sameCluster) ||
      (node->childR && !node->childR->sameCluster)
   ) return;

   // otherwise check if node->sameCluster flag needs updating
   size_t commonCluster = ds.find_set(node->left);
   if (node->childL) {
      size_t currentCluster = ds.find_set(node->childL->left);
      if (currentCluster != commonCluster) return; // not ready yet
   }
   if (node->childR) {
      size_t currentCluster = ds.find_set(node->childR->left);
      if (currentCluster != commonCluster) return; // not ready yet
   }
   node->sameCluster = true;
}


HeapNeighborItem HClustBiVpTreeMedoid::getNearestNeighbor(size_t index)
{
   if (shouldFind[index] && nearestNeighbors[index].empty())
   {
      size_t clusterIndex = ds.find_set(index);
      double _tau = INFINITY;//maxRadiuses[index];

#ifdef GENERATE_STATS
#ifdef _OPENMP
#pragma omp atomic
#endif
      ++stats.nnCals;
#endif
      /*if(!prefetch)
      for (size_t i=0; i<_n; ++i)
         medoidFound[i]=false;*/
      NNHeap nnheap(opts.maxNNPrefetch);
      getNearestNeighborsFromMinRadiusRecursive(_root, index, clusterIndex, minRadiuses[index], _tau, nnheap);
      nnheap.fill(nearestNeighbors[index]);

      // maxRadiuses[index] = INFINITY;
      size_t newNeighborsCount = nearestNeighbors[index].size();

      neighborsCount[index] += newNeighborsCount;
      if (neighborsCount[index] > _n - index || newNeighborsCount == 0)
         shouldFind[index] = false;

      if (newNeighborsCount > 0)
         minRadiuses[index] = nearestNeighbors[index].back().dist;
   }

   if (!nearestNeighbors[index].empty())
   {
#ifdef GENERATE_STATS
#ifdef _OPENMP
#pragma omp atomic
#endif
      ++stats.nnCount;
#endif
      auto res = nearestNeighbors[index].front();
      nearestNeighbors[index].pop_front();
      return res;
   }
   else
   {
      return HeapNeighborItem(SIZE_MAX,-INFINITY);
   }
}

/*HClustBiVpTreeMedoid::HeapHierarchicalItemMedoid HClustBiVpTreeMedoid::calculateCluster2ClusterMedoidDistance(size_t item1, size_t item2, size_t iter)
{
   size_t s1 = ds.find_set(item1);
   size_t s2 = ds.find_set(item2);
   HeapHierarchicalItemMedoid ret(s1,s2, -1, iter);
   //auto med1 = medoids.find(s1);
   //auto med2 = medoids.find(s2);
   //ret.dist = (*_distance)(med1->second, med2->second);
   return 0;
//   return ret;
}*/

size_t HClustBiVpTreeMedoid::mergeTwoClusters(size_t s1, size_t s2)
{
   size_t medoid1 = s1;
   size_t medoid2 = s2;

   double sumMedoid1Cluster2 = 0;
   double sumMedoid2Cluster1 = 0;

   double R = (*_distance)(_indices[medoid1], _indices[medoid2]);

   for (auto element=ds.getClusterMembers(s2).begin(); element != ds.getClusterMembers(s2).end(); ++element)
   {
      sumMedoid1Cluster2 += (*_distance)(_indices[*element], _indices[medoid1]);    
   }

   for (auto element=ds.getClusterMembers(s1).begin(); element != ds.getClusterMembers(s1).end(); ++element)
   {
      sumMedoid2Cluster1 += (*_distance)(_indices[*element], _indices[medoid2]);    
   }
   RCOUT("sumMedoid1Cluster2 = " << sumMedoid1Cluster2,3);
   RCOUT("sumMedoid2Cluster1 = " << sumMedoid2Cluster1,3);
   //szukam kandydata na medoid z klasta 1
   size_t candidate1 = SIZE_MAX;
   double distCandidate1 = INFINITY;
   for (auto element=ds.getClusterMembers(s1).begin(); element != ds.getClusterMembers(s1).end(); ++element)
   {
      if((*_distance)(_indices[*element], _indices[medoid1]) > R || (*_distance)(_indices[*element], _indices[medoid2]) > R) continue;
      double sumDist = 0.0;
      bool badOne = false;
      for (auto element2=ds.getClusterMembers(s2).begin(); element2 != ds.getClusterMembers(s2).end(); ++element2)
      {
         sumDist += (*_distance)(_indices[*element], _indices[*element2]);
         if(sumDist > sumMedoid1Cluster2) 
         {
            badOne = true;
            break;
         }
      }
      if(badOne) continue;
      for (auto element2=ds.getClusterMembers(s1).begin(); element2 != ds.getClusterMembers(s1).end(); ++element2)
      {
         sumDist += (*_distance)(_indices[*element], _indices[*element2]);
      }
      RCOUT("sumDist = " << sumDist, 3);
      if(sumDist < distCandidate1)
      {
         candidate1 = *element;
         distCandidate1 = sumDist;
      }
   }
   //szukam kandydata na medoid z klasta 2
   size_t candidate2 = SIZE_MAX;
   double distCandidate2 = INFINITY;
   for (auto element=ds.getClusterMembers(s2).begin(); element != ds.getClusterMembers(s2).end(); ++element)
   {
      if((*_distance)(_indices[*element], _indices[medoid1]) > R || (*_distance)(_indices[*element], _indices[medoid2]) > R) continue;
      double sumDist = 0.0;
      bool badOne = false;
      for (auto element2=ds.getClusterMembers(s1).begin(); element2 != ds.getClusterMembers(s1).end(); ++element2)
      {
         sumDist += (*_distance)(_indices[*element], _indices[*element2]);
         if(sumDist > sumMedoid2Cluster1) 
         {
            badOne = true;
            break;
         }
      }
      if(badOne) continue;
      for (auto element2=ds.getClusterMembers(s2).begin(); element2 != ds.getClusterMembers(s2).end(); ++element2)
      {
         sumDist += (*_distance)(_indices[*element], _indices[*element2]);
      }
      RCOUT("sumDist = " << sumDist, 3);
      if(sumDist < distCandidate2)
      {
         candidate2 = *element;
         distCandidate2 = sumDist;
      }
   }
   if(distCandidate1 < distCandidate2)
   {
      return candidate1;
   }
   else if(distCandidate1 == distCandidate2)
   {
      if(_indices[candidate1] < _indices[candidate2])
         return candidate1;
      else
         return candidate2;
   }
   else
   {
      return candidate2;
   }
}

size_t HClustBiVpTreeMedoid::medoidForCluster(size_t s)
{
   size_t indexMedoid = -1;
   /*double medoidSumDists = INFINITY;
   for (auto element = ds.getClusterMembers(s).begin(); element != ds.getClusterMembers(s).end(); ++element)
   {
      double sumdists = 0.0;
      Rcout << ".";
      for (auto element2 = ds.getClusterMembers(s).begin(); element2 != ds.getClusterMembers(s).end(); ++element2)
      {
         sumdists += (*_distance)(_indices[(*element)], _indices[(*element2)]);
      }     
      if (sumdists < medoidSumDists)
      {
         indexMedoid = (*element);
         medoidSumDists = sumdists;
      }
      if(sumdists == medoidSumDists)
      {
         if(_indices[*element] < _indices[indexMedoid])
            indexMedoid = *element;
      }
   }
   Rcout << endl;
   return indexMedoid;*/

   indexMedoid = *(ds.getClusterMembers(s).begin());
   for (auto element = ds.getClusterMembers(s).begin(); element != ds.getClusterMembers(s).end(); ++element)
   {
      if(_indices[indexMedoid] > _indices[*element])
         indexMedoid = *element;
   }
   return indexMedoid;
}

NumericMatrix HClustBiVpTreeMedoid::compute()
{
   //Rcout << "wchodze do medoid!" << endl;
   NumericMatrix ret(_n-1, 2);
   priority_queue<HeapHierarchicalItem> pq;

   //for (size_t i=0; i<_n; ++i)
   //   medoids[i] = i;

   // INIT: Pre-fetch a few nearest neighbors for each point
#if VERBOSE > 1
   Rprintf("[%010.3f] prefetching NNs\n", clock()/(float)CLOCKS_PER_SEC);
#endif

   prefetch = true;
#ifdef _OPENMP
   omp_set_dynamic(0); /* the runtime will not dynamically adjust the number of threads */
   omp_lock_t writelock;
   omp_init_lock(&writelock);
   #pragma omp parallel for schedule(dynamic)
#endif
   for (size_t i=0; i<_n; i++)
   {
#ifndef _OPENMP
   Rcpp::checkUserInterrupt(); // may throw an exception, fast op, not thread safe
#endif
      HeapNeighborItem hi=getNearestNeighbor(i);

      if (hi.index != SIZE_MAX)
      {
#if VERBOSE > 7 && !defined(_OPENMP)
         Rprintf("\r             prefetch NN: %d/%d", i, _n-1);
#endif
#ifdef _OPENMP
         omp_set_lock(&writelock);
#endif
         RCOUT("dla " << _indices[i]+1 << " najblizszym sasiadem jest " << _indices[hi.index]+1, 5);
         pq.push(HeapHierarchicalItem(i, hi.index, hi.dist));
         nnback[hi.index].push_front(i);
#ifdef _OPENMP
         omp_unset_lock(&writelock);
#endif
      }
   }
#ifdef _OPENMP
   omp_destroy_lock(&writelock);
#endif
#if VERBOSE > 7
   Rprintf("\r             prefetch NN: %d/%d\n", _n-1, _n-1);
#endif
#if VERBOSE > 1
   Rprintf("[%010.3f] merging clusters\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   bool awaria = false;
   prefetch = false;
   size_t i = 0;
   timestamp = 1;
   while(true)
   {
      Rcout << "iteracja " << i << endl;
      Rcout << "pq size = " << pq.size()<< endl;
      if(pq.size() > 200*_n) {awaria = true; break;}
      HeapHierarchicalItem hhi = pq.top();
      pq.pop();

      RCOUT("robie ds find dla " << _indices[hhi.index1]+1 << " " << _indices[hhi.index2]+1, 5);
      size_t s1 = ds.find_set(hhi.index1);
      size_t s2 = ds.find_set(hhi.index2);
      RCOUT("ds find udany", 5);
      if (s1 != hhi.index1 || s2 != hhi.index2) continue;
      
      Rcpp::checkUserInterrupt(); // may throw an exception, fast op

      ret(i,0)=(double)_indices[hhi.index1];
      ret(i,1)=(double)_indices[hhi.index2];
      //if(opts.medoidUpdateMethod == 1)
      //{
      //   medoids[s2] = medoids[s1] = mergeTwoClusters(s1, s2);
      //   //RCOUT("metoda sprytna zwrocila " << _indices[medoids[s1]]+1, 3);
      //}
      size_t s3;
      if(true || opts.medoidUpdateMethod == 0)
      {
         RCOUT("Rozpoczynam metode naiwna", 8);
         s3 = medoidForCluster(s1);
         //medoids[s2] = medoids[s1] = medoidNaive;
         RCOUT("Koncze metode naiwna", 8);
      }
      RCOUT("robie ds link dla " << _indices[s1]+1 << " " << _indices[s2]+1 << " " << _indices[s3]+1 ,5);
      ds.link(s1, s2, s3);
      RCOUT("ds link udany",5);

      ++i;
      timestamp++;
      RCOUT("po polaczeniu " << _indices[s1]+1<< " i " << _indices[s2]+1  << " nowym medoidem jest " << _indices[s3]+1, 5);
      if (i == _n-1) break; /* avoid computing unnecessary nn */
      
      HeapNeighborItem hi=getNearestNeighbor(s3);
      if (hi.index != SIZE_MAX)
      {
         nnback[hi.index].push_front(s3);
         RCOUT("1: dla " << _indices[s3]+1  << " najblizszym sasiadem jest " << _indices[hi.index]+1, 5);
         pq.push(HeapHierarchicalItem(s3, hi.index, hi.dist));
         //RCOUT("Znalazlem najblizszego sasiada",8);
      }
      if (s1 != s3)
      {
         for (auto it = nnback[s1].begin(); it != nnback[s1].end(); ++it)
         {
            if(ds.find_set(*it) != *it) continue;
            if(*it == s3) continue;
            double d_old = (*_distance)(*it, s1); // czy to nie bedzie liczone zbyt wiele razy to samo?
            double d_new = (*_distance)(*it, s3);
            size_t p;
            double d_p;
            if (d_new <= d_old)
            {
               p = s3;      
               d_p = d_new;
            }
            else
            {
               HeapNeighborItem hi2=getNearestNeighbor(*it);
               p = hi2.index;
               d_p = hi2.dist;
            }
            if(p != SIZE_MAX)
            {
               RCOUT("2: dla " << _indices[*it]+1  << " najblizszym sasiadem jest " << _indices[hi.index]+1, 5);
               nnback[p].push_front(*it);
               pq.push(HeapHierarchicalItem(*it, p, d_p));
            }
         }
         nnback[s1].clear();
      }
      if (s2 != s3)
      {
         for (auto it = nnback[s2].begin(); it != nnback[s2].end(); ++it)
         {
            if(ds.find_set(*it) != *it) continue;
            if(*it == s3) continue;
            double d_old = (*_distance)(*it, s2); // czy to nie bedzie liczone zbyt wiele razy to samo?
            double d_new = (*_distance)(*it, s3);
            size_t p;
            double d_p;
            if (d_new <= d_old)
            {
               p = s3;      
               d_p = d_new;
            }
            else
            {
               HeapNeighborItem hi2=getNearestNeighbor(*it);
               p = hi2.index;
               d_p = hi2.dist;
            }
            if(p != SIZE_MAX)
            {
               RCOUT("3: dla " << _indices[*it]+1  << " najblizszym sasiadem jest " << _indices[hi.index]+1, 5);
               nnback[p].push_front(*it);
               pq.push(HeapHierarchicalItem(*it, p, d_p));
            }
         }
         nnback[s2].clear();
      }
   }
#if VERBOSE > 7
      if (i % 1024 == 0) Rprintf("\r             %d / %d", i+1, _n);
#endif
   
#if VERBOSE > 7
   Rprintf("\r             %d / %d\n", _n, _n);
#endif
   Rcpp::checkUserInterrupt();

#if VERBOSE > 5
   Rprintf("[%010.3f] generating output matrix\n", clock()/(float)CLOCKS_PER_SEC);
#endif
   if(!awaria)
   {
      MergeMatrixGenerator mmg(ret.nrow());
      return mmg.generateMergeMatrix(ret);
   }
   else
      return ret;
}


void HClustBiVpTreeMedoid::print(HClustBiVpTreeMedoidNode* n) {
   if (n->childL) {
      Rprintf("\"%llx\" -> \"%llx\" [label=\"L\"];\n",
         (unsigned long long)n, (unsigned long long)(n->childL));
      print(n->childL);
   }
   if (n->childR) {
      Rprintf("\"%llx\" -> \"%llx\" [label=\"R\"];\n",
         (unsigned long long)n, (unsigned long long)(n->childR));
      print(n->childR);
   }

   if (n->vpindex == SIZE_MAX) {
      for (size_t i=n->left; i<n->right; ++i)
         Rprintf("\"%llx\" -> \"%llu\" [arrowhead = diamond];\n", (unsigned long long)n, (unsigned long long)_indices[i]+1);
   }
   else {
      Rprintf("\"%llx\" [label=\"(%llu, %g)\"];\n", (unsigned long long)n, (unsigned long long)n->vpindex+1, n->radius);
   }
}


void HClustBiVpTreeMedoid::print() {
   Rprintf("digraph vptree {\n");
   Rprintf("size=\"6,6\";\n");
   Rprintf("node [color=lightblue2, style=filled];");
   print(_root);
   Rprintf("}\n");
}


// [[Rcpp::export(".hclust2_medoid")]]
RObject hclust2_medoid(RObject distance, RObject objects, RObject control=R_NilValue) {
#if VERBOSE > 5
   Rprintf("[%010.3f] starting timer\n", clock()/(double)CLOCKS_PER_SEC);
#endif
   RObject result(R_NilValue);
   DataStructures::Distance* dist = DataStructures::Distance::createDistance(distance, objects);

   try {
      /* Rcpp::checkUserInterrupt(); may throw an exception */
      DataStructures::HClustBiVpTreeMedoid hclust(dist, control);
      RObject merge = hclust.compute();
      result = Rcpp::as<RObject>(List::create(
         _["merge"]  = merge,
         _["height"] = R_NilValue,
         _["order"]  = R_NilValue,
         _["labels"] = R_NilValue,
         _["call"]   = R_NilValue,
         _["method"] = "medoid",
         _["dist.method"] = R_NilValue,
         _["stats"] = List::create(
            _["vptree"] = hclust.getStats().toR(),
            _["distance"] = dist->getStats().toR()
         ),
         _["control"] = List::create(
            _["vptree"] = hclust.getOptions().toR()
         )
      ));
      result.attr("class") = "hclust";
      //hclust.print();
   }
   catch(...) {

   }

   if (dist) delete dist;
#if VERBOSE > 5
   Rprintf("[%010.3f] done\n", clock()/(double)CLOCKS_PER_SEC);
#endif
   if (Rf_isNull(result)) stop("stopping on error or explicit user interrupt");
   return result;
}

