// MIT License Terms (http://en.wikipedia.org/wiki/MIT_License)
//
// Copyright (C) 2011 by Oleg Alexandrov
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
#ifndef DTREE_H
#define DTREE_H

#include <vector>
#include <algorithm>
#include <cfloat> // defines DBL_MAX
#include <geomUtils.h>

// Trees for storing double precision (as opposed to integer) geometry.

// boxTree
// * Put Manhattan boxes in a tree.
// * Fast access to all boxes intersecting a given rectangular region.
// * Implemented as a template. The user defines the box class. Can be
//   for example a plain box, or a box with id meant to store an
//   edge (as in edgeTree below).
//
// Usage:
// -----
//   void formTreeOfBoxes(std::vector<Box> & Boxes // This vector is reordered inside
//                        );
//
//   void getBoxesInRegion(double xl, double yl, double xh, double yh, // Inputs
//                         vector<Box> & outBoxes                      // Outputs
//                         );
//

// edgeTree
// * Uses boxTree.
// * Put the edges of a given set of polygons in a tree.
// * Fast access to all edges intersecting a given rectangular region.
// * Fast access to the closest edge to a given point.
//
// Usage:
// -----
// void putPolyEdgesInTree(const polygonset & poly);
// void findPolyEdgesInBox(// inputs
//                         double xl, double yl,
//                         double xh, double yh,
//                         // outputs
//                         vector<seg> & edgesInBox
//                         );
// void findClosestEdgeToPoint(// inputs
//                             double x0, double y0,
//                             // outputs
//                             seg    & closestEdge,
//                             double & closestDist,
//                             // The location on the closest
//                             // edge where closestDist is achieved
//                             double & closestX, double & closestY
//                             );

namespace utils {
  class dPoly;
}

template <typename Box>
inline bool leftLessThan (Box P, Box Q){
  return (P.xl + P.xh) < (Q.xl + Q.xh); // x midpoint comparison
}

template <typename Box>
inline bool botLessThan (Box P, Box Q){
  return (P.yl + P.yh) < (Q.yl + Q.yh); // y midpoint comparison
}

template <typename Box>
inline bool lexLessThan(Box P, Box Q){
  if (P.xl < Q.xl) return true; if (P.xl > Q.xl) return false;
  if (P.yl < Q.yl) return true; if (P.yl > Q.yl) return false;
  if (P.xh < Q.xh) return true; if (P.xh > Q.xh) return false;
  if (P.yh < Q.yh) return true; if (P.yh > Q.yh) return false;
  return false;
}

// The MS compiler does not like these
// template <typename Box>
// inline bool operator==(Box P, Box Q){
//   return (P.xl == Q.xl && P.xh == Q.xh && P.yl == Q.yl && P.yh == Q.yh);
// }
// template <typename Box>
// inline bool operator!=(Box P, Box Q){
//   return (!(P == Q));
// }

template <typename Box>
struct boxNode{
  boxNode<Box> * left;
  boxNode<Box> * right;
  Box       Rect;
  bool      isLeftRightSplit;
  double    maxInLeftChild, minInRightChild;
  boxNode<Box>(): left(NULL), right(NULL), isLeftRightSplit(false),
                  maxInLeftChild(-DBL_MAX), minInRightChild(DBL_MAX){}
};

template <typename Box>
class boxTree{

public:
  boxTree();
  void formTreeOfBoxes(// Boxes will be reordered but otherwise
                       // unchanged inside this function
                       std::vector<Box> & Boxes);

  void getBoxesInRegion(double xl, double yl, double xh, double yh, // Inputs
                        std::vector<Box> & outBoxes                 // Outputs
                        );

  boxNode<Box> * getTreeRoot(){
    return m_root;
  }

private:

  void formTreeOfBoxesInternal(Box * Boxes, int numBoxes,
                               bool isLeftRightSplit,
                               boxNode<Box> *&  root
                               );


  void getBoxesInRegionInternal(//Inputs
                                double xl, double yl, double xh, double yh,
                                boxNode<Box> * root,
                                // Outputs
                                std::vector<Box> & outBoxes);

  void reset();
  boxNode<Box> * getNewboxNode();

  // Will get nodes from this pool (for performance reasons)
  std::vector< boxNode<Box> > m_nodePool;

  int m_freeNodeIndex;
  boxNode<Box> * m_root;

};

template <typename Box>
boxTree<Box>::boxTree(){
  reset();
  return;
}

template <typename Box>
void boxTree<Box>::reset(){
  m_freeNodeIndex = 0;
  m_root          = NULL;
  m_nodePool.clear();
  return;
}

template <typename Box>
boxNode<Box> * boxTree<Box>::getNewboxNode(){
  // Get a node from the pool
  boxNode<Box> * ptr = vecPtr(m_nodePool) + m_freeNodeIndex;
  assert( m_freeNodeIndex < (int)m_nodePool.size() );
  m_freeNodeIndex++;
  return ptr;
}

template <typename Box>
void boxTree<Box>::formTreeOfBoxes(// Boxes will be reordered but otherwise
                                   // unchanged inside this function
                                   std::vector<Box> & Boxes){

  reset();
  int numBoxes = Boxes.size();
  m_nodePool.resize(numBoxes);

  bool isLeftRightSplit = true;
  formTreeOfBoxesInternal(vecPtr(Boxes), numBoxes, isLeftRightSplit, m_root);
  return;
}

template <typename Box>
void boxTree<Box>::formTreeOfBoxesInternal(Box * Boxes, int numBoxes,
                                           bool isLeftRightSplit,
                                           boxNode<Box> *&  root
                                           ){


  // To do: Implement this without recursion.

  // To do: No need to store the box B in the tree. Store just a
  // pointer to B as sorting the left and right halves does not change
  // the midpoint box B.

  // To do: No need even for the tree, it can be stored in-place in
  // Boxes, in the same way as an in-place heap is stored.

  assert(numBoxes >= 0);

  if (numBoxes == 0){
    root = NULL;
    return;
  }

  root = getNewboxNode();
  root->isLeftRightSplit = isLeftRightSplit;

  if (isLeftRightSplit){ // Left-most boxes go in the left child
    std::sort(Boxes, Boxes + numBoxes, leftLessThan<Box>);
  }else{               // Bottom-most boxes go in the left child
    std::sort(Boxes, Boxes + numBoxes, botLessThan<Box>);
  }

  int mid = numBoxes/2;
  assert( 0 <= mid && mid < numBoxes);

  root->Rect = Boxes[mid]; // Must happen after sorting

  // A box is not a point, it has non-zero width. A box whose midpoint
  // is to the left of the root box midpoint, may still overlap with
  // the region to the right of the root box midpoint. As such, when
  // searching for boxes in the tree, we will need to consider how far
  // to the right the left subtree goes, and how far to the left the
  // right subtree goes.

  double maxInLeftChild = -DBL_MAX, minInRightChild = DBL_MAX;
  if (isLeftRightSplit){

    for (int s = 0; s < mid; s++)
      if (Boxes[s].xh >= maxInLeftChild) maxInLeftChild = Boxes[s].xh;

    for (int s = mid + 1; s < numBoxes; s++)
      if (Boxes[s].xl <= minInRightChild) minInRightChild = Boxes[s].xl;

  }else{

    for (int s = 0; s < mid; s++)
      if (Boxes[s].yh >= maxInLeftChild) maxInLeftChild = Boxes[s].yh;

    for (int s = mid + 1; s < numBoxes; s++)
      if (Boxes[s].yl <= minInRightChild) minInRightChild = Boxes[s].yl;

  }

  root->maxInLeftChild  = maxInLeftChild;
  root->minInRightChild = minInRightChild;

  // At the next split we will split perpendicularly to the direction
  // of the current split
  formTreeOfBoxesInternal( Boxes,            mid,
                           !isLeftRightSplit, root->left  );
  formTreeOfBoxesInternal( Boxes + mid + 1,  numBoxes - mid - 1,
                           !isLeftRightSplit, root->right );

  return;
}

template <typename Box>
void boxTree<Box>::getBoxesInRegion(// Inputs
                                    double xl, double yl, double xh, double yh,
                                    // Output
                                    std::vector<Box> & outBoxes
                                    ){

  outBoxes.clear();
  if (xl > xh || yl > yh || m_root == NULL) return;
  getBoxesInRegionInternal(xl, yl, xh, yh, m_root, // Inputs
                           outBoxes                // Outputs
                           );
  return;
}

template <typename Box>
void boxTree<Box>::getBoxesInRegionInternal(// Inputs
                                            double xl, double yl,
                                            double xh, double yh,
                                            boxNode<Box> * root,
                                            // Outputs
                                            std::vector<Box> & outBoxes
                                            ){

  assert (root != NULL);

  const Box & B = root->Rect; // alias

  if (utils::boxesIntersect(B.xl, B.yl, B.xh, B.yh, xl, yl, xh, yh))
    outBoxes.push_back(B);

  if (root->isLeftRightSplit){

    if (root->left != NULL  && xl <= root->maxInLeftChild)
      getBoxesInRegionInternal(xl, yl, xh, yh, root->left,  outBoxes);
    if (root->right != NULL && xh >= root->minInRightChild)
      getBoxesInRegionInternal(xl, yl, xh, yh, root->right, outBoxes);

  }else{

    if (root->left != NULL  && yl <= root->maxInLeftChild)
      getBoxesInRegionInternal(xl, yl, xh, yh, root->left,  outBoxes);
    if (root->right != NULL && yh >= root->minInRightChild)
      getBoxesInRegionInternal(xl, yl, xh, yh, root->right, outBoxes);

  }

  return;
}

class edgeTree{

public:

  void putPolyEdgesInTree(const utils::dPoly & poly);

  void findPolyEdgesInBox(// inputs
                          double xl, double yl,
                          double xh, double yh,
                          // outputs
                          std::vector<utils::seg> & edgesInBox);

  void findClosestEdgeToPoint(// inputs
                              double x0, double y0,
                              // outputs
                              utils::seg & closestEdge,
                              double     & closestDist,
                              // The location on the closest
                              // edge where closestDist is achieved
                              double & closestX, double & closestY
                              );

private:

  inline void edgeToBox(//inputs
                        double bx, double by, double ex, double ey,
                        // outputs
                        utils::dRectWithId & R
                        ){
    // Given an edge, compute its bounding box and the info necessary info
    // to reverse this transformation, see boxToEdge.
    int id = 0;
    if (bx > ex){ std::swap(bx, ex); id |= 1; } // id = id | 01
    if (by > ey){ std::swap(by, ey); id |= 2; } // id = id | 10
    R = utils::dRectWithId(bx, by, ex, ey, id);
    return;
  }

  inline void boxToEdge(// inputs
                        const utils::dRectWithId & R,
                        // outputs
                        double & bx, double & by, double & ex, double & ey
                        ){
    // Recover the edge based on the bounding box and id. The inverse of edgeToBox.
    bx = R.xl; by = R.yl; ex = R.xh; ey = R.yh;
    int id = R.id;
    if (id & 1) std::swap (bx, ex);
    if (id & 2) std::swap (by, ey);
  }

  void findClosestEdgeToPointInternal(// inputs
                                      double x0, double y0,
                                      boxNode<utils::dRectWithId> * root,
                                      // outputs
                                      utils::seg & closestEdge,
                                      double     & closestDist);
  // Internal data structures
  boxTree<utils::dRectWithId>      m_boxTree;
  std::vector<utils::dRectWithId>  m_allEdges;
  std::vector<utils::dRectWithId>  m_boxesInRegion;

};


#endif
