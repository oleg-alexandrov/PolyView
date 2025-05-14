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
#ifndef DPOLY_H
#define DPOLY_H

#include <vector>
#include <map>
#include <baseUtils.h>
#include <geomUtils.h>
#include "dTree.h"
#include "kdTree.h"

namespace utils {
  
enum AnnoType {
  fileAnno = 0, vertAnno, polyAnno, layerAnno, angleAnno, lastAnno
};

// A class holding a set of polygons in double precision
class dPoly {

public:

  dPoly(){
    reset();
  }

  void reset();

  bool read_pol_or_cnt_format(std::string filename,
                              std::string type,
                              bool isPointCloud = false);

  bool readPoly(std::string filename,
                bool isPointCloud = false
                );

  void writePoly(std::string filename, std::string defaultColor = "yellow");
  void bdBoxCenter(double & mx, double & my) const;

  void appendPolygon(int numVerts,
                     const double * xv,
                     const double * yv,
                     bool isPolyClosed,
                     const std::string & color,
                     const std::string & layer
                     );

  void appendPolygons(const dPoly & poly);

  void appendRectangle(double xl, double yl, double xh, double yh,
                       bool isPolyClosed,
                       const std::string & color, const std::string & layer
                       );

  void getEdges(std::vector<seg> &allEdges) const;

  void setRectangle(double xl, double yl, double xh, double yh,
                    bool isPolyClosed,
                    const std::string & color, const std::string & layer
                    );

  bool isXYRect();

  void clipPointCloud(const dRect &clip_box,
                      dPoly & clippedPoly, // output
                      const std::vector<int> *selected);

  void clipPolygons(const dRect &clip_box,
                    dPoly & clippedPoly, // output
                    const std::vector<int> *selected);

  void clipAll(double clip_xll, double clip_yll,
                double clip_xur, double clip_yur,
                dPoly & clippedPoly, // output
                const std::vector<int> *selected = nullptr // optional input, if provided only clip selected polygons
  );

  void clipAnno(const dRect &clip_box,
                dPoly & clippedPoly);

  void copyAnno(dPoly & clippedPoly);

  std::vector<dPoint> getAcuteAngleLocs(double angle);
  std::vector<dPoint> getNonManhLocs() const;
  std::vector<dPoint> getNon45Locs() const;
  std::vector<dPoint> getDuplicates() const;

  void shift(double shift_x, double shift_y);
  void rotate(double angle);
  void scale(double scale);
  void transformMarkedPolys(std::vector<int> const& mark, const linTrans & T);
  void transformMarkedAnnos(std::vector<int> const& mark, const linTrans & T);

  void transformMarkedPolysAroundPt(std::vector<int> const& mark,
                                    const matrix2 & M, dPoint P);
  void transformMarkedAnnosAroundPt(std::vector<int> const& mark,
                                    const matrix2 & M, dPoint P);

  void applyTransform(double a11, double a12, double a21, double a22,
                      double sx, double sy,
                      linTrans & T); // save the transform here
  
  void applyTransformAroundBdBoxCenter(double a11, double a12,
                                       double a21, double a22,
                                       linTrans & T
                                       );

  const int    * get_numVerts         () const { return vecPtr(m_numVerts); }
  const double * get_xv               () const { return vecPtr(m_xv);       }
  const double * get_yv               () const { return vecPtr(m_yv);       }

  // Non-const versions of the above
  int    * get_numVerts         () { return vecPtr(m_numVerts); }
  double * get_xv               () { return vecPtr(m_xv);       }
  double * get_yv               () { return vecPtr(m_yv);       }
  
  int get_numPolys                    () const { return m_numPolys;                }
  int get_totalNumVerts               () const { return m_totalNumVerts;           }
  std::vector<char> get_isPolyClosed  () const { return m_isPolyClosed;            }
  std::vector<std::string> get_colors () const { return m_colors;                  }
  std::vector<std::string> get_layers () const { return m_layers;                  }

  void set_color(std::string color);

  void set_isPolyClosed(bool isPolyClosed);

  void eraseMarkedPolys(std::vector<int> const& mark);
  void eraseMarkedAnnos(std::vector<int> const& mark);

  void erasePolysIntersectingBox(double xll, double yll, double xur, double yur);
  void eraseAnnosIntersectingBox(double xll, double yll, double xur, double yur);
  void appendAndShiftMarkedPolys(// Inputs
                                 std::vector<int> & mark,
                                 double shift_x, double shift_y
                                 );
  void set_isPointCloud(bool isPointCloud){ m_isPointCloud = isPointCloud; }
  bool isPointCloud() const { return m_isPointCloud;}

  void set_pointCloud(const std::vector<dPoint> & P, std::string color,
                      std::string layer);
  void buildGrid(double xl, double yl, double xh, double yh,
                 double gridSize, std::string gridColor);

  void markPointsInBox(// Inputs
      double xll, double yll,
      double xur, double yur,
      // Outputs
      std::vector<int> & mark) const;

  void markPolysIntersectingBox(// Inputs
                                double xll, double yll,
                                double xur, double yur,
                                // Outputs
                                std::vector<int> & mark) const;

  void markAnnosIntersectingBox(// Inputs
                                double xll, double yll,
                                double xur, double yur,
                                // Outputs
                                std::vector<int> & mark) const;
  
  //void replaceOnePoly(int polyIndex, int numV, const double* x, const double* y);
  // Annotations
  std::vector<anno>&  get_annotations()  { return m_annotations;}
  std::vector<anno>&  get_vertIndexAnno(){ return m_vertIndexAnno;}
  std::vector<anno>&  get_polyIndexAnno(){ return m_polyIndexAnno;}
  std::vector<anno>&  get_layerAnno()    { return m_layerAnno;}
  std::vector<anno>&  get_angleAnno()    { return m_angleAnno;}

  const std::vector<anno>&  get_annotations()  const { return m_annotations;}
  const std::vector<anno>&  get_vertIndexAnno()const { return m_vertIndexAnno;}
  const std::vector<anno>&  get_polyIndexAnno()const { return m_polyIndexAnno;}
  const std::vector<anno>&  get_layerAnno()    const {return m_layerAnno;}

  void set_annotations(const std::vector<anno> & A);
  void set_layerAnno(const std::vector<anno> & annotations);
  void set_vertIndexAnno(const std::vector<anno> & annotations);
  void set_polyIndexAnno(const std::vector<anno> & annotations);

  void compVertFullIndexAnno();
  void compVertIndexAnno();
  void compPolyIndexAnno();
  void compLayerAnno();
  void updateBoundingBox() const;

  const dRect& bdBox() const;
  void bdBox(double & xll, double & yll, double & xur, double & yur) const;
  void annoBdBox(double & xll, double & yll, double & xur, double & yur) const;
  dRect annoBdBox() const;

  void bdBoxes(std::vector<double> & xll, std::vector<double> & yll,
               std::vector<double> & xur, std::vector<double> & yur) const;

  void setPolygon(int numVerts,
                  const double * xv,
                  const double * yv,
                  bool isPolyClosed,
                  const std::string & color,
                  const std::string & layer
                  );

  void eraseAnno(int annoIndex);

  void findClosestAnnotation(// inputs
                             double x0, double y0,
                             // outputs
                             int & annoIndex,
                             double & min_dist
                             ) const;

  void findClosestPolyVertex(// inputs
                             double x0, double y0,
                             // outputs
                             int & polyIndex,
                             int & vertIndex,
                             double & min_x, double & min_y,
                             double & min_dist
                             ) const;

  std::pair<std::complex<double>, std::complex<double>>
  getClosestPolyEdge( double x0, double y0, double &minDist) const;

  void findClosestPolyEdge(//inputs
                           double x0, double y0,
                           // outputs
                           int & polyIndex, int & vertIndex,
                           double & minX, double & minY, double & minDist
                           ) const;

  void eraseOnePoly(int polyIndex);
  void insertVertex(int polyIndex, int vertIndex,
                    double x, double y);
  void eraseVertex(int polyIndex, int vertIndex);
  void changeVertexValue(int polyIndex, int vertIndex, double x, double y);
  void shiftEdge(int polyIndex, int vertIndex, double shift_x, double shift_y);
  void shiftOnePoly(int polyIndex, double shift_x, double shift_y);
  void shiftOneAnno(int index, double shift_x, double shift_y);
  void shiftMarkedPolys(std::vector<int> const & mark, double shift_x, double shift_y);
  void shiftMarkedAnnos(std::vector<int> const & amark, double shift_x, double shift_y);

  void reverseMarkedPolys(std::vector<int> const & mark);

  void extractOnePoly(int polyIndex,       // input
                      dPoly & poly,  // output
					  int start_index = -1) const; // if start index provided do not re-compute
                      
  void extractMarkedPolys(std::vector<int> const& mark, // input
                          dPoly & polys) const;           // output
  
  // Reverse orientation of all polygons
  void reverse();
  
  void reverseOnePoly(int polyIndex);
  void sortFromLargestToSmallest(bool counter_cc);

  void sortBySizeAndMaybeAddBigContainingRect(// inputs
                                              double bigXll, double bigYll,
                                              double bigXur, double bigYur,
                                              bool counter_cc);

  void enforce45();

  std::vector<int> getPolyIdsInBox(const dRect &box) const;

  // Pointer to a structure storing an image and its positioning
  // information in world coordinates. The image can be shown in the
  // background of the given polygon. The user is responsible
  // for managing this correctly. This approach makes it possible
  // for this polygon class not to know a lot of GUI-related constructs,
  // and for not having to allocate memory for an image, which can be
  // an issue if polygons are copied around, such as in an undo buffer.
  void * img;

  // This will create the bounding box tree and return a pointer to it
  // If three is already created it will just return the pointer
  const boxTree< dRectWithId>  * getBoundingBoxTree() const;
  const kdTree * getPointTree() const;

private:

  // Clear pre-computed data if geometry changes
  void clearExtraData();

  bool getColorInCntFile(const std::string & line, std::string & color);
  std::vector<anno> &  get_annoByType(AnnoType annoType);
  void set_annoByType(const std::vector<anno> & annotations, AnnoType annoType);
  const std::vector<int> & getStartingIndices() const;
  // If isPointCloud is true, treat each point as a set of unconnected points
  bool                     m_isPointCloud;
  bool                     m_has_color_in_file;
  std::vector<double>      m_xv;
  std::vector<double>      m_yv;
  std::vector<int>         m_numVerts;
  int                      m_numPolys;
  int                      m_totalNumVerts;
  std::vector<char>        m_isPolyClosed;
  std::vector<std::string> m_colors;
  std::vector<std::string> m_layers;
  std::vector<anno>        m_annotations;
  std::vector<anno>        m_vertIndexAnno; // Anno showing vertex index
  std::vector<anno>        m_polyIndexAnno; // Anno showing poly index in a set of polys
  std::vector<anno>        m_layerAnno;     // Anno showing layer number
  std::vector<anno>        m_angleAnno;
  // The following items are used for performance.
  // They should be cleared if geometry changes by calling clearExtraData()
  mutable boxTree< dRectWithId> m_boundingBoxTree;
  mutable kdTree m_pointTree;

  mutable std::vector<int>  m_startingIndices;
  mutable dRect m_BoundingBox;
};

} // end namespace utils
#endif
