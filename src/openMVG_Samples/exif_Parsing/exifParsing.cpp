#include "openMVG/exif/exif_IO_EasyExif.hpp"
using namespace openMVG::exif;

#include "third_party/cmdLine/cmdLine.h"
#include <memory>

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sInputImage;

  cmd.add( make_option('i', sInputImage, "imafile") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << ' '
      << "[-i|--imafile path] "
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--imafile " << sInputImage << std::endl;
  
  std::unique_ptr<Exif_IO> exif_io( new Exif_IO_EasyExif( sInputImage ) );

  std::cout << "width : " << exif_io->getWidth() << std::endl;
  std::cout << "height : " << exif_io->getHeight() << std::endl;
  std::cout << "focal : " << exif_io->getFocal() << std::endl;
  std::cout << "brand : " << exif_io->getBrand() << std::endl;
  std::cout << "model : " << exif_io->getModel() << std::endl;
  return EXIT_SUCCESS;
}



