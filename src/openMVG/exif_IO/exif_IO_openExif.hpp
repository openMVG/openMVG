#ifndef EXIF_IO_OPENEXIF_HPP
#define EXIF_IO_OPENEXIF_HPP

#include "exif_IO.hpp"

#include <dependencies/OpenExif/src/ExifImageFile.h>
#include <iostream>

class Exif_IO_OpenExif: public Exif_IO
{
  public:
    Exif_IO_OpenExif() {}

    Exif_IO_OpenExif( const std::string & sFileName )
    {
      open( sFileName );
    }

    size_t getWidth() const
    {
      return info_.width;
    }

    size_t getHeight() const
    {
      return info_.height;
    }

    float getFocal() const
    {
      float *f;
      return (f = getTag<float>(EXIFTAG_FOCALLENGTH)) == NULL ? -1 : *f;
    }

    std::string getBrand() const
    {
      std::string *s; return (s = getTag<std::string>(EXIFTAG_MAKE)) == NULL ? "Not found" : *s;
    }

    std::string getModel() const
    {
      std::string *s; return (s = getTag<std::string>(EXIFTAG_MODEL)) == NULL ? "Not found" : *s;
    }

    std::string getLensModel() const
    {
      return "";
    }

    /** Open the file for checking and parsing */
    bool open( const std::string & sFileName )
    {
      sfilename = sFileName;
      /////////////////////////////////////////////////////
      //Get data exif
      /////////////////////////////////////////////////////
      app1PathsTags_.clear();
      app3PathsTags_.clear();
      appSegs_.clear();

      bool isOpenWell = false;

      // Open the file in read-only mode and verify that it succeeds
      if (inImageFile_.open( sFileName.c_str(), "r" ) == EXIF_OK)
      {
        //-- Read EXIF data
        {
          // Get all the AppSeg 1 - "Exif" tags and output them
          inImageFile_.getAllTags( 0xFFE1, "Exif", app1PathsTags_ );

          // Get all the AppSeg 3 - "Meta" tags and output them
          inImageFile_.getAllTags( 0xFFE3, "Meta", app3PathsTags_ );

          //Now, recognition of any other app segments:
          // Get a vector with all the application segments in the file
          appSegs_ = inImageFile_.getAllAppSegs();

          // Now, lets output any COM marker data
          comList_ = (ExifComMarkerList*)inImageFile_.getComData();

          // And finally, let's output the SOF info
          inImageFile_.getImageInfo(info_);
        }
        isOpenWell = true;
      }
      return isOpenWell;
    }

    /**Verify if the file has metadata*/
    bool doesHaveExifInfo() const
    {
      bool bRet = false;
      // Create instance of ExifImageFile
      ExifImageFile inImageFile;

      // Open the file in read-only mode and verify that it succeeds
      if (inImageFile.open( sfilename.c_str(), "r" ) == EXIF_OK)
      {
        if( inImageFile.close() != EXIF_OK )
        {
          std::cerr << "Error: Could not close" << sfilename << std::endl;
        }
        else
        {
          bRet = true;
        }
      }
      else
        std::cerr << "File can't be opened" << std::endl;
      return bRet;
    }

    /** Print all data*/
    std::string allExifData() const
    {
       /////////////////////////////////////////////////////
      //Print data exif
      /////////////////////////////////////////////////////

      // Display "Exif tags"
      if (app1PathsTags_.begin() != app1PathsTags_.end())
        std::cout << "\nApp1 - \"Exif\" entries:" << std::endl;
      for (ExifPathsTags::const_iterator crntPathsTags = app1PathsTags_.begin();
           crntPathsTags != app1PathsTags_.end();
           crntPathsTags++ )
      {
        ExifIFDPath::const_iterator crntPath = (*crntPathsTags).first.begin();
        ExifIFDPath::const_iterator endPath = (*crntPathsTags).first.end();
        while( crntPath != endPath )
        {
          std::cout << "IFD: " << (*crntPath).first
               << "  Idx: " << (*crntPath).second << std::endl;
          crntPath++;
        }

        ExifTags::const_iterator crnt = (*crntPathsTags).second.begin();
        ExifTags::const_iterator end = (*crntPathsTags).second.end();

        std::cout << "Tag#\tType\tCount\tValue" << std::endl;
        while( crnt != end )
        {
          ExifTagEntry* tag = *(crnt);
          tag->print();
          std::cout << std::endl;
          crnt++;
        }
      }

      // Display Meta tags
      if (app3PathsTags_.begin() != app3PathsTags_.end())
        std::cout << "\nApp3 - \"Meta\" entries:" << std::endl;
      for (ExifPathsTags::const_iterator crntPathsTags = app3PathsTags_.begin();
           crntPathsTags != app3PathsTags_.end();
           crntPathsTags++ )
      {
        ExifIFDPath::const_iterator crntPath = (*crntPathsTags).first.begin();
        ExifIFDPath::const_iterator endPath = (*crntPathsTags).first.end();
        while( crntPath != endPath )
        {
          std::cout << "IFD: " << (*crntPath).first
               << "  Idx: " << (*crntPath).second << std::endl;
          crntPath++;
        }

        ExifTags::const_iterator crnt = (*crntPathsTags).second.begin();
        ExifTags::const_iterator end = (*crntPathsTags).second.end();
        std::cout << "Tag#\tType\tCount\tValue" << std::endl;
        while( crnt != end )
        {
          ExifTagEntry* tag = *(crnt);
          tag->print();
          std::cout << std::endl;
          crnt++;
        }
      }

      // Application Segments
      size_t numOfAppSegs = appSegs_.size();

      std::cout << "\n\nNumber of Application Segments "
           << ": " << numOfAppSegs << std::endl
           << "Marker\tLength\tIdent" << std::endl;

      // Loop through the application segments outputting their marker,
      // length and identifier.
      for ( size_t i = 0; i < numOfAppSegs; ++i )
      {
          std::cout << appSegs_[i]->getAppSegmentMarker() << "\t"
               << appSegs_[i]->getLength() << "\t"
               << appSegs_[i]->getAppIdent() << std::endl;
      }

      // And finally, let's output the SOF info
      std::cout << "Image Information:\n";
      std::cout << "\twidth:\t\t" << info_.width << std::endl
           << "\theight:\t\t" << info_.height << std::endl
           << "\tchannels:\t" << info_.numChannels << std::endl
           << "\tbit depth:\t" << info_.precision << std::endl;

      std::cout << std::endl
           << "Summation" << std::endl
           << " Exif data occurences   : " << app1PathsTags_.size() << std::endl
           << " Meta data occurences   : " << app3PathsTags_.size() << std::endl
           << " App segment occurences : " << appSegs_.size() << std::endl
           << " Com marker  occurences : " << comList_->size() << std::endl;
	  return std::string("");
    }


  private:
    /**Get the value which numTag designate*/
    template <typename T>
    T* getTag(const int numTag)const
    {
      for (ExifPathsTags::const_iterator crntPathsTags = app1PathsTags_.begin();
          crntPathsTags != app1PathsTags_.end();
          crntPathsTags++ )
      {
        ExifTags::const_iterator crnt = (*crntPathsTags).second.begin();
        ExifTags::const_iterator end = (*crntPathsTags).second.end();
        while( crnt != end )
        {
          ExifTagEntry* tag = *(crnt);

          if (tag->getTagNum() == numTag)
          {
            return &(static_cast<ExifTagEntryT<T>*>(tag))->getValue();
          }
          crnt++;
        }
      }
      return NULL;
    }

    //Variables
    std::string sfilename;
    ExifImageFile outImageFile_;
    ExifImageFile inImageFile_;
    ExifPathsTags app1PathsTags_;     // App Seg 1 "Exif tags"
    ExifPathsTags app3PathsTags_;     // App Seg 3 "Meta tags"
    std::vector<ExifAppSegment*> appSegs_; // App segments.
    ExifComMarkerList* comList_;     // COM marker data
    ExifImageInfo info_;
};
#endif //EXIF_IO_OPENEXIF_HPP

