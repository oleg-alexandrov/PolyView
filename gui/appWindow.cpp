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
#include <qapplication.h>
// #include <QMenu>
#include <qlabel.h>
#include <QtGui>
#include <qmenubar.h>
#include <qmessagebox.h>
#include <qstatusbar.h>
#include <qlayout.h>
#include <QMenu>
#include <QEvent>
#include <QTextEdit>
#include <QUrl>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <appWindow.h>
#include <chooseFilesDlg.h>
#include <polyView.h>
#include <utils.h>

using namespace std;

cmdLine::cmdLine(QWidget* parent): QLineEdit(parent){}
cmdLine::~cmdLine(){}

appWindow::appWindow(QWidget* parent, std::string progName,
                     const cmdLineOptions & options,
                     int windowWidX, int windowWidY
                     ): QMainWindow(parent, Qt::Window), m_poly(NULL){

  installEventFilter(this);

  m_progName = progName;
  resize(windowWidX, windowWidY);

  createMenusAndMainWidget(options);

  // Command line
  m_cmdLine = new cmdLine(this);
  m_cmdLine->setAlignment(Qt::AlignLeft);
  m_cmdLine->setFocusPolicy(Qt::StrongFocus);
  connect( m_cmdLine, SIGNAL( returnPressed() ),
          this, SLOT( procCmdLine() ) );
  QStatusBar * status = statusBar();
  QRect Rp = status->rect();
  m_cmdLine->setGeometry(Rp);
  status->addWidget(m_cmdLine, 1);
  m_cmdHist.clear();
  m_histPos = 0;

  return;
}

void appWindow::closeEvent(QCloseEvent *){
  forceQuit();
}

void appWindow::forceQuit(){
  exit(0); // A fix for an older buggy version of Qt
}

bool appWindow::eventFilter(QObject *obj, QEvent *event){

  // If the alt or control key is hit, and the focus is on the command
  // line widget, move the focus to the m_poly widget.
  if ( event->type() == QEvent::KeyPress ||
       event->type() == QEvent::ShortcutOverride){
    QKeyEvent * keyEvent = (QKeyEvent*)event;
    Qt::KeyboardModifiers modifiers = keyEvent->modifiers();
    if ( ( modifiers & Qt::ControlModifier ) ||
         ( modifiers & Qt::AltModifier ) ){
      if (m_cmdLine->hasFocus() && m_poly != NULL) m_poly->setFocus();
    }
  }

  if (obj == m_poly) {
    // Avoid repainting on these events
    if (event->type() == QEvent::FocusIn          ||
        event->type() == QEvent::FocusOut         ||
        event->type() == QEvent::WindowDeactivate ||
        event->type() == QEvent::Leave
        ){
      return true;
    }
  }

  //cout << "Other event: " << (int)event->type() << endl;
  return QWidget::eventFilter(obj, event);
}

appWindow::~appWindow(){

  if (m_poly != NULL){ delete m_poly; m_poly = NULL; }

}

void appWindow::procCmdLine(){
  string cmd = m_cmdLine->text().toStdString();
  m_cmdHist.push_back(cmd);
  m_poly->runCmd(cmd);
  m_cmdLine->setText("");
  m_histPos = m_cmdHist.size();
  //m_poly->setFocus();
}

void appWindow::insertCmdFromHist(){

  int numHistItems = m_cmdHist.size();
  m_histPos = min(m_histPos, numHistItems);
  m_histPos = max(0, m_histPos);

  if (m_histPos < numHistItems){
    m_cmdLine->setText(m_cmdHist[m_histPos].c_str());
  }else{
    m_cmdLine->setText("");
  }

  return;
}

void appWindow::shiftUp (){

  if (m_poly->hasFocus()){
    m_poly->shiftUp ();
  }else if(m_cmdLine->hasFocus()){
    m_histPos--;
    insertCmdFromHist();
  }

}

void appWindow::shiftDown (){

  if (m_poly->hasFocus()){
    m_poly->shiftDown ();
  }else if(m_cmdLine->hasFocus()){
    m_histPos++;
    insertCmdFromHist();
  }

}

void appWindow::createMenusAndMainWidget(const cmdLineOptions & opt){

  // There is some twisted logic here. First initialize the menus,
  // then create the main widget, then finish creating the menus.
  // This is a workaround for a bug in a certain versions of Qt. If
  // the main widget is created before the menus then it gets
  // incorrect geometry.

  QMenuBar* menu = menuBar();
  QMenu* file = new QMenu(QString("File"));
  menu->addMenu(file);

  // Central widget
  m_poly = new polyView (this, opt);
  m_poly->setFocusPolicy(Qt::StrongFocus);
  m_poly->setFocus();
  setCentralWidget(m_poly);

  file->addAction("Open", m_poly, SLOT(openPoly()), Qt::CTRL+Qt::Key_O);
  file->addAction("Save as a combined file", m_poly, SLOT(saveOnePoly()),
                   Qt::CTRL+Qt::Key_S);
  file->addAction("Save as individual files", m_poly,
                   SLOT(saveAsMultiplePolys()), Qt::ALT+Qt::Key_S);
  file->addAction("Overwrite current files", m_poly,
                   SLOT(overwriteMultiplePolys()), Qt::CTRL+Qt::Key_W);
  file->addAction("Reload polygons from disk", m_poly,
                   SLOT(reloadPolys()), Qt::ALT+Qt::Key_R);
  file->addAction("Exit", this, SLOT(forceQuit()), Qt::Key_Q);

  QMenu* view = new QMenu("View", menu );
  menu->addMenu( view);
  view->addAction(chooseFilesDlg::selectFilesTag(), m_poly, SLOT(chooseFilesToShow()));
  view->addAction("Zoom out",             m_poly, SLOT(zoomOut()),      Qt::Key_Minus);
  view->addAction("Zoom in",              m_poly, SLOT(zoomIn()),       Qt::Key_Equal);
  view->addAction("Move left",            m_poly, SLOT(shiftLeft()),    Qt::Key_Left);
  view->addAction("Move right",           m_poly, SLOT(shiftRight()),   Qt::Key_Right);
  view->addAction("Move up",              this,   SLOT(shiftUp()),      Qt::Key_Up);
  view->addAction("Move down",            this,   SLOT(shiftDown()),    Qt::Key_Down);
  view->addAction("Reset view",           m_poly, SLOT(resetView()),    Qt::Key_R);
  view->addAction("Change display order", m_poly, SLOT(changeOrder()),  Qt::Key_O);
  view->addAction("Toggle annotations",   m_poly, SLOT(toggleAnno()),   Qt::Key_A);
  view->addAction("Toggle filled",        m_poly, SLOT(toggleFilled()), Qt::Key_F);
  view->addAction("Toggle points display",
                   m_poly, SLOT(togglePE()),          Qt::Key_P);
  view->addAction("Toggle show vertex indices",
                   m_poly, SLOT(toggleVertIndexAnno()), Qt::Key_V);
  view->addAction("Toggle show layers", m_poly, SLOT(toggleLayerAnno()), Qt::Key_L);

  QMenu* edit = new QMenu("Edit", menu );
  menu->addMenu( edit);
  edit->addAction("Undo",  m_poly, SLOT(undo()), Qt::Key_Z);
  edit->addAction("Redo",  m_poly, SLOT(redo()), Qt::ALT + Qt::Key_Z);

  edit->addAction("Create poly with int vertices and 45x deg angles",
                   m_poly, SLOT(create45DegreeIntPoly()), Qt::Key_N);
  edit->addAction("Enforce int vertices and 45x deg angles", m_poly, SLOT(enforce45()),
                   Qt::CTRL+Qt::Key_4);
  edit->addAction("Create arbitrary polygon",
                   m_poly, SLOT(createArbitraryPoly()), Qt::CTRL+Qt::Key_N);
  edit->addAction("Merge polygons (buggy)",
                   m_poly, SLOT(mergePolys()), Qt::CTRL+Qt::Key_M);

  QMenu* selection = new QMenu("Selection", menu );
  menu->addMenu( selection);
  selection->addAction("Create highlight", m_poly, SLOT(createHlt()));
  selection->addAction("Cut polygons to highlight",
                        m_poly, SLOT(cutToHlt()), Qt::Key_C);
  selection->addAction("Delete selected polygons",
                        m_poly, SLOT(deleteSelectedPolys()), Qt::CTRL+Qt::Key_D);
  selection->addAction("Paste selected polygons",
                        m_poly, SLOT(pasteSelectedPolys()), Qt::CTRL+Qt::Key_V);
  selection->addAction("Move selected polygons",
                        m_poly, SLOT(moveSelectedPolys()));
  selection->addAction("Deselect polygons/delete highlights",  m_poly,
                        SLOT(deselectPolysDeleteHlts()));

  QMenu* transform = new QMenu("Transform", menu );
  menu->addMenu( transform);
  transform->addAction("Translate selected polygons", m_poly,
                        SLOT(translateSelectedPolys()), Qt::CTRL+Qt::Key_T );
  transform->addAction("Rotate selected polygons",    m_poly,
                        SLOT(rotateSelectedPolys()), Qt::CTRL+Qt::Key_R );
  transform->addAction("Scale selected polygons",
                        m_poly, SLOT(scaleSelectedPolys()) );
  transform->addAction("Transform selected polygons", m_poly,
                        SLOT(transformSelectedPolys()) );
#if 0
  selection->addSeparator();
  transform->insertItem("Translate polygons", m_poly, SLOT(translatePolys()),
                        Qt::CTRL+Qt::Key_T);
  transform->insertItem("Rotate polygons around origin", m_poly, SLOT(rotatePolys()),
                        Qt::CTRL+Qt::Key_R);
  transform->insertItem("Scale polygons around origin", m_poly, SLOT(scalePolys()),
                        Qt::CTRL+Qt::Key_X);
#endif

  QMenu* grid = new QMenu("Grid", menu );
  menu->addMenu( grid);
  grid->addAction("Toggle poly grid", m_poly, SLOT(toggleShowGrid()), Qt::Key_G);
  grid->addAction("Enforce 45x deg angles and snap to grid", m_poly, SLOT(enforce45AndSnapToGrid()));
  grid->addAction("Set grid size", m_poly, SLOT(setGridSize()));
  grid->addAction("Set grid linewidth", m_poly, SLOT(setGridWidth()));
  grid->addAction("Set grid color", m_poly, SLOT(setGridColor()));

  QMenu* diff = new QMenu("Diff", menu );
  menu->addMenu( diff);
  diff->addAction("Toggle different colors", m_poly, SLOT(toggleDifferentColors()), Qt::SHIFT + Qt::Key_D);
  diff->addAction("Toggle show poly diff", m_poly, SLOT(toggleShowPolyDiff()), Qt::Key_D);
  diff->addAction("Show next diff", m_poly, SLOT(plotNextDiff()), Qt::Key_K);
  diff->addAction("Show prev diff", m_poly, SLOT(plotPrevDiff()), Qt::Key_J);

  QMenu* options = new QMenu("Options", menu );
  menu->addMenu( options);
  options->addAction("Set line width", m_poly, SLOT(setLineWidth()));
  options->addAction("Set background color", m_poly, SLOT(setBgColor()));

  QMenu* help = new QMenu("Help", menu );
  menu->addMenu( help);
  
  // Hide the doc as it is hard to regenerate the html each
  // time. GitHub has the latest documentation.
  //help->addAction("Show documentation", this, SLOT(showDoc()));

  help->addAction("About", this, SLOT(about()));

  return;
}

void appWindow::showDoc(){
  m_docWindow.setText(utils::getDocText());
  m_docWindow.show();
  return;
}

void appWindow::about(){

  QMessageBox msgBox;
  msgBox.setText("Author: Oleg Alexandrov");
  msgBox.setInformativeText("See the documentation on GitHub.");
  msgBox.setStandardButtons(QMessageBox::Ok);
  msgBox.setDefaultButton(QMessageBox::Ok);
  msgBox.exec();

  return;
}

docWindow::docWindow(QWidget *){
  resize(900, 800);
  m_textArea = new QTextEdit (this);
  setCentralWidget (m_textArea);
}

docWindow::~docWindow(){
  delete m_textArea;
  m_textArea = NULL;
}

void docWindow::setText(const std::string & docText){

  QImage img(":/pvLogo.png");
  m_textArea->document()->addResource(QTextDocument::ImageResource,
                                      QUrl("pvLogo.png" ), img);
  m_textArea->clear();
  m_textArea->setReadOnly(true);
  m_textArea->setCurrentFont(QFont("Monospace", 10));
  m_textArea->insertHtml(docText.c_str());
  m_textArea->moveCursor(QTextCursor::Start);

  return;
}
