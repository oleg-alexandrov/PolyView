TEMPLATE        = app
TARGET          = polyview
QT             += qt3support
RESOURCES     = polyview.qrc
SOURCES = mainProg.cpp polyView.cpp appWindow.cpp chooseFilesDlg.cpp utils.cpp documentation.cpp ../geom/dPoly.cpp ../geom/cutPoly.cpp ../geom/geomUtils.cpp ../geom/polyUtils.cpp ../geom/edgeUtils.cpp ../geom/dTree.cpp ../geom/kdTree.cpp
HEADERS = polyView.h appWindow.h chooseFilesDlg.h utils.h ../geom/dPoly.h ../geom/cutPoly.h  ../geom/geomUtils.h ../geom/polyUtils.h ../geom/edgeUtils.h ../geom/dTree.h ../geom/kdTree.h ../geom/baseUtils.h


