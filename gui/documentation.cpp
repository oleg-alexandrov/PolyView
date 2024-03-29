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
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cassert>
#include <cstring>
#include <gui/utils.h>
using namespace std;

std::string utils::getDocText(){

  // Begin auto-generated text
char docText[] ="<div style=\"color:rgb(0, 0, 0)\">\n"
"<div dir=\"ltr\">\n"
"<div style=\"background-color:rgb(255, 255, 255)\" align=\"center\"><a href=:pvLogo.png imageanchor=\"1\"><img alt=\"PolyView logo\" src=:pvLogo.png border=\"0\"></a></div>\n"
"<div style=\"background-color:rgb(255, 255, 255)\" align=\"left\">\n"
"<br>\n"
"<br>\n"
"\n"
"<p>PolyView is a free and open source cross-platform program designed to\n"
"quickly load, view, and edit multiple files containing polygons. It\n"
"can zoom and pan, show polygons as edges, points, and filled, display\n"
"text labels, print the coordinates of vertices and measure distances,\n"
"change the order in which polygons are displayed, choose which\n"
"polygons to show, etc.</p>\n"
"\n"
"<p>PolyView can add, delete, move, rotate, and scale polygons, add,\n"
"delete, and move vertices and edges, add and delete text labels, and\n"
"cut polygons to a box.</p>\n"
"\n"
"<h1><a id=\"user-content-download\" class=\"anchor\" href=\"#download\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Download</h1>\n"
"\n"
"<ul>\n"
"<li><a href=\"https://github.com/oleg-alexandrov/PolyView/releases\">https://github.com/oleg-alexandrov/PolyView/releases</a></li>\n"
"</ul>\n"
"\n"
"<p>On Linux and OSX the PolyView program is installed in the <em>bin</em>\n"
"directory.  On Windows it is in the base directory. The Windows\n"
"version is large because of the size of the Qt libraries.</p>\n"
"\n"
"<h1><a id=\"user-content-license\" class=\"anchor\" href=\"#license\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>License</h1>\n"
"\n"
"<p>PolyView is available under the MIT license and can be used for any\n"
"purpose, academic or commercial.</p>\n"
"\n"
"<h1><a id=\"user-content-compiling\" class=\"anchor\" href=\"#compiling\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Compiling</h1>\n"
"\n"
"<p>PolyView is written in C++. It was successfully compiled on Linux with\n"
"g++, on OSX with clang, and on Windows with MinGW. It was written with\n"
"portability in mind and it should be possible to compile it on any\n"
"platform and compiler.</p>\n"
"\n"
"<p>Its only dependency is the Qt 4 library (e.g., Qt 4.8). Instructions\n"
"on how to compile Qt are given at the end of this document.</p>\n"
"\n"
"<p>To compile PolyView on Linux or OSX, run:</p>\n"
"\n"
"<pre><code>qmake polyview.pro\n"
"make\n"
"</code></pre>\n"
"\n"
"<h1><a id=\"user-content-documentation\" class=\"anchor\" href=\"#documentation\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Comparison with XGRAPH</h1>\n"
"\n"
"<p>PolyView was inspired by xgraph (<a\n"
"href=\"http://www.xgraph.org/\">http://www.xgraph.org/</a>).  The latter\n"
"is a general purpose data plotter and has a lot of features PolyView\n"
"lacks. PolyView has extra features for viewing and editing polygons.\n"
"\n"
"<p>PolyView is (as of 2007) more responsive for polygons with very\n"
"many vertices, and it can handle zooming to small regions of polygons\n"
"with large floating point vertex coordinates without overflowing and\n"
"showing incorrect results. Credit for responsiveness goes to Qt, and\n"
"issues with overlow required careful handling.\n"
"\n"
"<p>Lastly, PolyView is open-source and under a liberal license, and can\n"
"be improved in a collaborative manner.\n"
"\n"
"<h1><a id=\"user-content-documentation\" class=\"anchor\" href=\"#documentation\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Documentation</h1>\n"
"\n"
"<h2><a id=\"user-content-file-format\" class=\"anchor\" href=\"#file-format\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>File format</h2>\n"
"\n"
"<p>PolyView uses the XGRAPH file format, in which the x and y\n"
"coordinates of each polygon are stored as two columns in a plain text\n"
"file followed by the \"NEXT\" statement. Here is a sample polygon file\n"
"storing a rectangle and a triangle.</p>\n"
"\n"
"<pre><code>color = red\n"
"6.5 5.0\n"
"7.2 5.0\n"
"7.2 7.0\n"
"6.5 7.0\n"
"6.5 5.0\n"
"NEXT\n"
"color = green\n"
"8.0 3.0\n"
"9.0 3.0\n"
"8.0 4.0\n"
"8.0 3.0\n"
"NEXT\n"
"anno 7.0 5.5 text label\n"
"</code></pre>\n"
"\n"
"<p>Polygon vertices are stored in double precision.</p>\n"
"\n"
"<p>Notice that the first vertex of each polygon is repeated again before\n"
"the \"NEXT\" statement, to signal to the viewer that the polygon is\n"
"closed. PolyView can also display polygonal lines, in addition to\n"
"polygons, when the last vertex need not equal the first one. The last line\n"
"above shows how to place a text label, with its coordinates and text.</p>\n"
"\n"
"<h2><a id=\"user-content-functionality\" class=\"anchor\" href=\"#functionality\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Functionality</h2>\n"
"\n"
"<h3><a id=\"user-content-mouse-buttons\" class=\"anchor\" href=\"#mouse-buttons\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Mouse buttons</h3>\n"
"\n"
"<p>The left mouse button snaps to the closest polygon vertex and prints\n"
"its coordinates in the terminal used to start PolyView. A subsequent\n"
"click also prints the distance from the previous vertex to the current\n"
"one.  The middle mouse button prints the coordinates of where the\n"
"mouse clicked, without snapping to the closest vertex.  Dragging the\n"
"mouse from lower-left to upper-right zooms in, and doing it in reverse\n"
"zooms out.  Dragging the mouse while keeping the Control key pressed\n"
"creates a highlight which can be used to cut the polygons to the\n"
"highlight or to paste/move/delete them.</p>\n"
"\n"
"<h3><a id=\"user-content-keyboard-shortcuts\" class=\"anchor\" href=\"#keyboard-shortcuts\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Keyboard shortcuts</h3>\n"
"\n"
"<p>Panning is accomplished by using the arrow keys, zooming is done with\n"
"the plus and minus keys. Many other actions are bound to keyboard\n"
"keys, as can be seen when invoking them from the menu.</p>\n"
"\n"
"<h3><a id=\"user-content-command-box\" class=\"anchor\" href=\"#command-box\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Command box</h3>\n"
"\n"
"<p>Many GUI operations (such as zoom) echo their action as a command in\n"
"the terminal. That command can be pasted in the command box at the\n"
"bottom of the PolyView GUI to reproduce the original operation. This\n"
"provides a basic level of scripting and reproducibility.</p>\n"
"\n"
"<h3><a id=\"user-content-menu-functions\" class=\"anchor\" href=\"#menu-functions\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Menu functions</h3>\n"
"\n"
"<h4><a id=\"user-content-file-menu\" class=\"anchor\" href=\"#file-menu\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>File menu</h4>\n"
"\n"
"<ul>\n"
"<li>Load a polygon file in addition to existing files</li>\n"
"<li>Save the polygons as one file</li>\n"
"<li>Save the polygons as individual files</li>\n"
"<li>Overwrite the existing polygons</li>\n"
"</ul>\n"
"\n"
"<h4><a id=\"user-content-view-menu\" class=\"anchor\" href=\"#view-menu\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>View menu</h4>\n"
"\n"
"<ul>\n"
"<li>Choose which files to hide/show</li>\n"
"<li>Zoom and pan</li>\n"
"<li>Reset the view to contain all polygons with a small padding</li>\n"
"<li>Change the order in which the polygons are displayed</li>\n"
"<li>Show/hide annotations (text labels)</li>\n"
"<li>Show the polygons filled with color</li>\n"
"<li>Show the vertices only (no edges)</li>\n"
"<li>Show the index of each vertex</li>\n"
"<li>Show the layer ids (if present)</li>\n"
"</ul>\n"
"\n"
"<h4><a id=\"user-content-edit-menu\" class=\"anchor\" href=\"#edit-menu\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Edit menu</h4>\n"
"\n"
"<ul>\n"
"<li>Undo and redo</li>\n"
"<li>Create a polygon with integer vertices and edge angles multiple of 45 degrees</li>\n"
"<li>Enforce integer vertices and edge angles multiple of 45 degrees</li>\n"
"<li>Create a polygon with arbitrary angles</li>\n"
"</ul>\n"
"\n"
"<h4><a id=\"user-content-transform-menu\" class=\"anchor\" href=\"#transform-menu\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Transform menu</h4>\n"
"\n"
"<ul>\n"
"<li>Translate/rotate/scale the polygons</li>\n"
"</ul>\n"
"\n"
"<h4><a id=\"user-content-selection-menu\" class=\"anchor\" href=\"#selection-menu\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Selection menu</h4>\n"
"\n"
"<ul>\n"
"<li>Create a highlight (the polygons in the highlight are automatically selected and copied to a buffer)</li>\n"
"<li>Cut the polygons to the current highlight</li>\n"
"<li>Delete the selected polygons</li>\n"
"<li>Paste the selected polygons</li>\n"
"<li>Move the selected polygons (use Shift-Mouse)</li>\n"
"<li>Deselect all polygons and delete all highlights</li>\n"
"<li>Translate/rotate/scale/transform selected polygons</li>\n"
"</ul>\n"
"\n"
"<h4><a id=\"user-content-grid-menu\" class=\"anchor\" href=\"#grid-menu\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Grid menu</h4>\n"
"\n"
"<ul>\n"
"<li>Show/hide the grid</li>\n"
"<li>Enforce that all polygon edges have angles multiple of 45 degrees and snap the vertices to the grid</li>\n"
"<li>Set the grid size</li>\n"
"<li>Set the grid linewidth</li>\n"
"<li>Set the grid color</li>\n"
"</ul>\n"
"\n"
"<h4><a id=\"user-content-diff-menu\" class=\"anchor\" href=\"#diff-menu\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Diff menu</h4>\n"
"\n"
"<ul>\n"
"<li>Change the colors of polygons so that polygons from different files have different colors</li>\n"
"<li>Enter diff mode (a mode in which two similar polygon files can be compared)</li>\n"
"<li>Show the next/previous difference between given two given polygon files (starting with the largest difference)</li>\n"
"</ul>\n"
"\n"
"<h4><a id=\"user-content-options-menu\" class=\"anchor\" href=\"#options-menu\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Options menu</h4>\n"
"\n"
"<ul>\n"
"<li>Set the linewidth of polygon edges</li>\n"
"<li>Set the background color</li>\n"
"</ul>\n"
"\n"
"<h4><a id=\"user-content-right-click-menu\" class=\"anchor\" href=\"#right-click-menu\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Right-click menu</h4>\n"
"\n"
"<ul>\n"
"<li>Show and save a mark at the current point</li>\n"
"<li>Create a polygon with integer vertices and edge angles multiple of 45 degrees</li>\n"
"<li>Enforce integer vertices and edge angles multiple of 45 degrees</li>\n"
"<li>Create a polygon with arbitrary angles</li>\n"
"<li>Delete the polygon at mouse cursor</li>\n"
"<li>Enter move polygons mode (use Shift-Mouse to move a polygon; if some polygons are selected using a highlight, then all selected polygons will be moved)</li>\n"
"<li>Enter move vertices mode (use Shift-Mouse to move a vertex)</li>\n"
"<li>Enter move edges mode (use Shift-Mouse to move an edge)</li>\n"
"<li>Insert a vertex on the edge closest to the mouse cursor</li>\n"
"<li>Delete the vertex closest to the mouse cursor</li>\n"
"<li>Copy the polygon closest to the mouse cursor</li>\n"
"<li>Paste the polygon closest to the mouse cursor</li>\n"
"<li>Reverse the orientation of the polygon closest to the mouse cursor</li>\n"
"<li>Insert text label</li>\n"
"<li>Delete text label</li>\n"
"<li>Enter align mode (a mode in which, given two polygon files, the second polygon file is kept fixed, while the first one can be interactively translated using Shift-Mouse and rotated/flipped from the right-click menu until it aligns to the second one)</li>\n"
"</ul>\n"
"\n"
"<h1><a id=\"user-content-command-line-options\" class=\"anchor\" href=\"#command-line-options\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Command line options</h1>\n"
"\n"
"<p>PolyView will open simultaneously all polygon files supplied as inputs\n"
"on the command line. Various command line options can modify how the\n"
"polygons are displayed.</p>\n"
"\n"
"<ul>\n"
"<li>-h | -help    Show the available command line options</li>\n"
"<li>-geo[metry] <code>width</code>x<code>height</code>  The window size in pixels (for example, 800x800)</li>\n"
"<li>-bg | -backgroundColor <code>color</code>    Color of the background (the default is black)</li>\n"
"<li>-c | -color <code>color</code>   All polygons after this option will show up in the given color (the default is to use the colors specified in the polygon files) </li>\n"
"<li>-nc | -noColorOverride  All polygons after this option will show up in the color specified in the file.</li>\n"
"<li>-fs | -fontSize <code>integer</code> The text font size in pixels</li>\n"
"<li>-lw | -lineWidth <code>integer</code>    All polygons after this option will show up with given linewidth</li>\n"
"<li>-p | -points  All polygons after this option will show up as vertices rather than edges (a subsequent -p option undoes this behavior)</li>\n"
"<li>-cp | -closedPoly All polygons after this option will show up as closed (the last vertex is connected to the first one)</li>\n"
"<li>-ncp | -nonClosedPoly  Interpret the polygons after this option as polygonal lines (the last vertex is not connected to the first one)</li>\n"
"<li>-f | -filledPoly  All polygons after this option will show up as filled</li>\n"
"<li>-cw | -clockwisePoly  Polygons oriented clockwise are assumed to have positive area</li>\n"
"<li>-nf | -nonFilledPoly  All polygons after this option will show up as not filled</li>\n"
"<li>-grid on | off    Turn on/off the grid display</li>\n"
"<li>-gridSize <code>integer</code>   Grid size</li>\n"
"<li>-gridWidth <code>integer</code>  Grid width in pixels</li>\n"
"<li>-gridColor <code>color</code>    Grid color</li>\n"
"</ul>\n"
"\n"
"<h1><a id=\"user-content-compiling-qt\" class=\"anchor\" href=\"#compiling-qt\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Compiling Qt</h1>\n"
"\n"
"<p>PolyView was tested to compile on Linux, OSX, and Windows with Qt 4.8 but other 4.x versions should work as well.</p>\n"
"\n"
"<p>Qt can be installed on Ubuntu with the command:</p>\n"
"\n"
"<pre><code>apt-get install qmake\n"
"</code></pre>\n"
"\n"
"<p>To compile Qt 4.8 from source, get the source code from:</p>\n"
"\n"
"<ul>\n"
"<li> <a href=\"https://download.qt.io/archive/qt/4.8/\">https://download.qt.io/archive/qt/4.8/</a></li>\n"
"</ul>\n"
"\n"
"<p>After unzipping the archive and ensuring that the compiler is in the\n"
"path, run:</p>\n"
"\n"
"<ul>\n"
"<li>configure -opensource -fast -confirm-license -nomake demos -nomake examples -nomake docs -nomake translations -no-webkit -no-script -no-scripttools -no-openssl -no-libjpeg -no-libmng  -no-libtiff -no-cups -no-nis -no-opengl -no-openvg -no-phonon -no-phonon-backend -no-sql-psql -no-dbus</li>\n"
"</ul>\n"
"\n"
"<h1><a id=\"user-content-author\" class=\"anchor\" href=\"#author\" aria-hidden=\"true\"><span class=\"octicon octicon-link\"></span></a>Author</h1>\n"
"\n"
"<p>Oleg Alexandrov (<a href=\"mailto:oleg.alexandrov@gmail.com\">oleg.alexandrov@gmail.com</a>)</p>\n";
  // End auto-generated text
 return docText;
}
