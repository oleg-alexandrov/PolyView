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
#include <Q3PopupMenu>
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
                     ): QMainWindow(parent, progName.c_str()), m_poly(NULL){

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
  Q3PopupMenu* file = new Q3PopupMenu( menu );
  menu->insertItem("File", file);

  // Central widget
  m_poly = new polyView (this, opt);
  m_poly->setFocusPolicy(Qt::StrongFocus);
  m_poly->setFocus();
  setCentralWidget(m_poly);

  file->insertItem("Open", m_poly, SLOT(openPoly()), Qt::CTRL+Qt::Key_O);
  file->insertItem("Save as a combined file", m_poly, SLOT(saveOnePoly()),
                   Qt::CTRL+Qt::Key_S);
  file->insertItem("Save as individual files", m_poly,
                   SLOT(saveAsMultiplePolys()), Qt::ALT+Qt::Key_S);
  file->insertItem("Overwrite current files", m_poly,
                   SLOT(overwriteMultiplePolys()), Qt::CTRL+Qt::Key_W);
  file->insertItem("Reload polygons from disk", m_poly,
                   SLOT(reloadPolys()), Qt::ALT+Qt::Key_R);
  file->insertItem("Exit", this, SLOT(forceQuit()), Qt::Key_Q);

  Q3PopupMenu* view = new Q3PopupMenu( menu );
  menu->insertItem("View", view);
  view->insertItem(chooseFilesDlg::selectFilesTag(), m_poly, SLOT(chooseFilesToShow()));
  view->insertItem("Zoom out",             m_poly, SLOT(zoomOut()),      Qt::Key_Minus);
  view->insertItem("Zoom in",              m_poly, SLOT(zoomIn()),       Qt::Key_Equal);
  view->insertItem("Move left",            m_poly, SLOT(shiftLeft()),    Qt::Key_Left);
  view->insertItem("Move right",           m_poly, SLOT(shiftRight()),   Qt::Key_Right);
  view->insertItem("Move up",              this,   SLOT(shiftUp()),      Qt::Key_Up);
  view->insertItem("Move down",            this,   SLOT(shiftDown()),    Qt::Key_Down);
  view->insertItem("Reset view",           m_poly, SLOT(resetView()),    Qt::Key_R);
  view->insertItem("Change display order", m_poly, SLOT(changeOrder()),  Qt::Key_O);
  view->insertItem("Toggle annotations",   m_poly, SLOT(toggleAnno()),   Qt::Key_A);
  view->insertItem("Toggle filled",        m_poly, SLOT(toggleFilled()), Qt::Key_F);
  view->insertItem("Toggle points display",
                   m_poly, SLOT(togglePE()),          Qt::Key_P);
  view->insertItem("Toggle show vertex indices",
                   m_poly, SLOT(toggleVertIndexAnno()), Qt::Key_V);
  view->insertItem("Toggle show layers", m_poly, SLOT(toggleLayerAnno()), Qt::Key_L);

  Q3PopupMenu* edit = new Q3PopupMenu( menu );
  menu->insertItem("Edit", edit);
  edit->insertItem("Undo",  m_poly, SLOT(undo()), Qt::Key_Z);
  edit->insertItem("Redo",  m_poly, SLOT(redo()), Qt::ALT + Qt::Key_Z);

  edit->insertItem("Create poly with int vertices and 45x deg angles",
                   m_poly, SLOT(create45DegreeIntPoly()), Qt::Key_N);
  edit->insertItem("Enforce int vertices and 45x deg angles", m_poly, SLOT(enforce45()),
                   Qt::CTRL+Qt::Key_4);
  edit->insertItem("Create arbitrary polygon",
                   m_poly, SLOT(createArbitraryPoly()), Qt::CTRL+Qt::Key_N);
  edit->insertItem("Merge polygons (buggy)",
                   m_poly, SLOT(mergePolys()), Qt::CTRL+Qt::Key_M);

  Q3PopupMenu* selection = new Q3PopupMenu( menu );
  menu->insertItem("Selection", selection);
  selection->insertItem("Create highlight", m_poly, SLOT(createHlt()));
  selection->insertItem("Cut polygons to highlight",
                        m_poly, SLOT(cutToHlt()), Qt::Key_C);
  selection->insertItem("Delete selected polygons",
                        m_poly, SLOT(deleteSelectedPolys()), Qt::CTRL+Qt::Key_D);
  selection->insertItem("Paste selected polygons",
                        m_poly, SLOT(pasteSelectedPolys()), Qt::CTRL+Qt::Key_V);
  selection->insertItem("Move selected polygons",
                        m_poly, SLOT(moveSelectedPolys()));
  selection->insertItem("Deselect polygons/delete highlights",  m_poly,
                        SLOT(deselectPolysDeleteHlts()));

  Q3PopupMenu* transform = new Q3PopupMenu( menu );
  menu->insertItem("Transform", transform);
  transform->insertItem("Translate selected polygons", m_poly,
                        SLOT(translateSelectedPolys()), Qt::CTRL+Qt::Key_T );
  transform->insertItem("Rotate selected polygons",    m_poly,
                        SLOT(rotateSelectedPolys()), Qt::CTRL+Qt::Key_R );
  transform->insertItem("Scale selected polygons",
                        m_poly, SLOT(scaleSelectedPolys()) );
  transform->insertItem("Transform selected polygons", m_poly,
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

  Q3PopupMenu* grid = new Q3PopupMenu( menu );
  menu->insertItem("Grid", grid);
  grid->insertItem("Toggle poly grid", m_poly, SLOT(toggleShowGrid()), Qt::Key_G);
  grid->insertItem("Enforce 45x deg angles and snap to grid", m_poly, SLOT(enforce45AndSnapToGrid()));
  grid->insertItem("Set grid size", m_poly, SLOT(setGridSize()));
  grid->insertItem("Set grid linewidth", m_poly, SLOT(setGridWidth()));
  grid->insertItem("Set grid color", m_poly, SLOT(setGridColor()));

  Q3PopupMenu* diff = new Q3PopupMenu( menu );
  menu->insertItem("Diff", diff);
  diff->insertItem("Toggle different colors", m_poly, SLOT(toggleDifferentColors()), Qt::SHIFT + Qt::Key_D);
  diff->insertItem("Toggle show poly diff", m_poly, SLOT(toggleShowPolyDiff()), Qt::Key_D);
  diff->insertItem("Show next diff", m_poly, SLOT(plotNextDiff()), Qt::Key_K);
  diff->insertItem("Show prev diff", m_poly, SLOT(plotPrevDiff()), Qt::Key_J);

  Q3PopupMenu* options = new Q3PopupMenu( menu );
  menu->insertItem("Options", options);
  options->insertItem("Set line width", m_poly, SLOT(setLineWidth()));
  options->insertItem("Set background color", m_poly, SLOT(setBgColor()));

  Q3PopupMenu* help = new Q3PopupMenu( menu );
  menu->insertItem("Help", help);
  help->insertItem("Show documentation", this, SLOT(showDoc()));
  help->insertItem("About", this, SLOT(about()));

  return;
}

void appWindow::showDoc(){
  m_docWindow.setText(utils::getDocText());
  m_docWindow.setCaption(this->caption()); // Borrow the caption from the parent
  m_docWindow.show();
  return;
}

void appWindow::about(){

  string aboutStr = string("About ") + m_progName;
  static QMessageBox* about
    = new QMessageBox( aboutStr.c_str(),
                       "© Oleg Alexandrov        ", // extra space to make window bigger
                       QMessageBox::NoIcon, 1, 0, 0, this, 0, FALSE );
  about->setButtonText( 1, "OK" );
  about->show();

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
