![ScreenShot](gui/pvLogo.png)

PolyView is a free and open source cross-platform program designed to
quickly load and visualize multiple sets of polygon files from the
command line. It can zoom and pan, show polygons as edges, points, and
filled, display text labels, print the coordinates of vertices and
measure distances, change the order in which polygons are displayed,
choose which polygons to show, etc.

PolyView can also do basic polygon editing, such as
creating/deleting/moving/rotating/scaling polygons,
adding/removing/modifying vertices, as well as cutting polygons to a
box.

PolyView accepts a lengthy list of command line arguments, and the
view can be also manipulated interactively from the GUI once the
polygons are loaded. The most frequent actions are bound to keyboard
shortcuts. A complete list of features and functions is shown below.

# License

PolyView is available under the MIT license and can be used for any
purpose, academic or commercial.

# Compiling

PolyView is written in C++. It was successfully compiled on Linux and
OSX using g++ and on Windows with MinGW. It was written with
portability in mind and it should be possible to compile it on any
platform and compiler.

Its only dependency is the Qt 4 development library, which can be
installed on Ubuntu with the command:

```
apt-get install qmake
```

To compile on Linux or OSX:

```
cd gui
qmake polyview.pro
make
```

The 'qmake' command will generate a makefile which will be used by
'make' to build the executable. The 'gui' directory has a sample
makefile as generated by 'qmake'.


# Documentation

## File format

PolyView was inspired by xgraph, and it uses the xgraph file format,
in which the x and y coordinates of each polygon are stored as two
columns in a plain text file followed by the "NEXT" statement. Here is
a sample polygon file storing a rectangle and a triangle.

```
color = red
60 50
70 50
70 70
60 70
60 50
NEXT
color = green
80 30
90 30
80 40
80 30
NEXT
```

Notice that the first vertex of each polygon is repeated again before
the "NEXT" statement, to signal to the viewer that the polygon is
closed.  

## GUI functionality

### Mouse buttons

The left mouse button snaps to the closest polygon vertex and prints
its coordinates. A subsequent click also prints the distance from the
previous vertex to the current one.  The middle mouse button prints
the coordinates of where the mouse clicked, without snapping to the
closest vertex.  Dragging the mouse from lower-left to upper-right
zooms in, and doing it in reverse zooms out.  Dragging the mouse while
keeping the Control key pressed creates a highlight which can be used
to cut the polygons to the highlight or to paste/move/delete them.

### Command box

Many GUI operations (such as zoom) echo their action as a command in
the terminal used to start PolyView. That command can be pasted in the
command box at the bottom of the PolyView GUI to reproduce the
original operation. This provides a basic level of scripting and
reproducibility.

### Menu functions 

#### File menu

* Load a polygon file in addition to existing files
* Save the polygons as one file
* Save the polygons as individual files
* Overwrite the existing polygons

#### View menu

* Choose which files to hide/show
* Zoom and pan
* Reset the view to contain all polygons with a small padding
* Change the order in which the polygons are displayed
* Show/hide annotations (text labels)
* Show the polygons filled with color
* Show the vertices only (no edges)
* Show the index of each vertex
* Show the layer ids (if present)

#### Edit menu

* Undo and redo
* Create a polygon with integer vertices and edge angles multiple of 45 degrees
* Enforce integer vertices and edge angles multiple of 45 degrees
* Create a polygon with arbitrary angles

#### Transform menu

* Translate/rotate/scale the polygons

#### Selection menu

* Create a highlight (the polygons in the highlight are automatically selected and copied to a buffer)
* Cut the polygons to the current highlight
* Delete the selected polygons
* Paste the selected polygons
* Move the selected polygons (use Shift-Mouse)
* Deselect all polygons and delete all highlights
* Translate/rotate/scale/transform selected polygons

#### Grid menu

* Show/hide the grid
* Enforce that all polygon edges have angles multiple of 45 degrees and snap the vertices to the grid
* Set the grid size
* Set the grid linewidth
* Set the grid color

#### Diff menu

* Change the colors of polygons so that polygons from different files have different colors
* Enter diff mode (a mode in which two similar polygon files can be compared)
* Show the next/previous difference between given two given polygon files (starting with the largest difference)

#### Options menu

* Set the linewidth of polygon edges
* Set the background color

#### Right-click menu

* Show and save a mark at the current point
* Create a polygon with integer vertices and edge angles multiple of 45 degrees
* Enforce integer vertices and edge angles multiple of 45 degrees
* Create a polygon with arbitrary angles
* Delete the polygon at mouse cursor
* Enter align mode (a mode in which, given two polygon files, the second polygon file is kept fixed, while the first one can be interactively translated using Shift-Mouse and rotated/flipped from the right-click menu until it aligns to the second one)
* Enter move polygons mode (use Shift-Mouse to move a polygon; if some polygons are selected using a highlight, then all selected polygons will be moved)
* Enter move vertices mode (use Shift-Mouse to move a vertex)
* Enter move edges mode (use Shift-Mouse to move an edge)
* Insert a vertex on the edge closest to the mouse cursor
* Delete the vertex closest to the mouse cursor
* Copy the polygon closest to the mouse cursor
* Paste the polygon closest to the mouse cursor
* Reverse the orientation of the polygon closest to the mouse cursor

# Command line options

PolyView will open simultaneously all polygon files supplied as inputs
on the command line. Various command line options can modify how the
polygons are displayed.

* -h | -help	Show the available command line options
* -geo[metry] <width>x<height>	The window size in pixels (for example, 800x800)
* -bg | -backgroundColor<color>	Color of the background (the default is black)
* -c | -color<color>	All polygons after this option will show up in the given color (the default is to use the colors specified in the polygon files) 
* -fs | -fontSize<integer>	The text font size in pixels
* -lw | -lineWidth<integer>	All polygons after this option will show up with given linewidth
* -p | -points	All polygons after this option will show up as vertices rather than edges (a subsequent -p option undoes this behavior)
* -cp | -closedPoly	All polygons after this option will show up as closed (the last vertex is connected to the first one)
* -nc | -nonClosedPoly	Interpret the polygons after this option as polygonal lines (the last vertex is not connected to the first one)
* -f | -filledPoly	All polygons after this option will show up as filled
* -nf | -nonFilledPoly	All polygons after this option will show up as not filled
* -grid on | off	Turn on/off the grid display
* -gridSize<integer>	Grid size
* -gridWidth<integer>	Grid width in pixels
* -gridColor<color>	Grid color

# Author

Oleg Alexandrov (oleg.alexandrov@gmail.com)
