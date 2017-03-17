#include "MainWindow.hh"

#include <QApplication>

using namespace image_undistort_gui ;

int main( int argc , char ** argv )
{
  QApplication app( argc , argv ) ; 

  MainWindow win ;
  win.show() ; 
  
  return app.exec() ;
}