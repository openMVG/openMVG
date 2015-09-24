#ifndef DATASHEET_HPP
#define DATASHEET_HPP

#include "openMVG/stl/split.hpp"
#include <iterator>
#include <algorithm>

// Database structure
struct Datasheet
{
  Datasheet()
  {}

  Datasheet( const std::string& brand,
             const std::string& model,
             const double& sensorSize ):
                  _brand(brand),
                  _model(model),
                  _sensorSize(sensorSize)
  {}

  bool operator==(const Datasheet& ds) const
  {
    bool isEqual = false;
    std::vector<std::string> vec_brand;
    stl::split(ds._brand, " ", vec_brand);
    std::string brandlower = _brand;

    std::transform(brandlower.begin(), brandlower.end(),
      brandlower.begin(), ::tolower);

    for ( std::vector<std::string>::const_iterator iter_brand = vec_brand.begin();
            iter_brand != vec_brand.end();
            iter_brand++ )
    {
      std::string brandlower2 = *iter_brand;
      std::transform(brandlower2.begin(), brandlower2.end(),
        brandlower2.begin(), ::tolower);
      //std::cout << brandlower << "\t" << brandlower2 << std::endl;
      if ( brandlower.compare( brandlower2 ) == 0 )
      {
        std::vector<std::string> vec_model1;
        stl::split(ds._model, " ", vec_model1);
        std::vector<std::string> vec_model2;
        stl::split(_model, " ", vec_model2);
        bool isAllFind = true;
        for ( std::vector<std::string>::const_iterator iter_model1 = vec_model1.begin();
            iter_model1 != vec_model1.end();
            iter_model1++ )
        {
          bool hasDigit = false;
          for(std::string::const_iterator c = (*iter_model1).begin(); c != (*iter_model1).end(); ++c )
          {
            if(isdigit(*c))
            {
              hasDigit = true;
              break;
            }
          }
          if ( hasDigit )
          {
            std::string modellower1 = *iter_model1;
            for ( int index = 0; index < modellower1.length(); index++ )
            {
              modellower1[index] = tolower(modellower1[index]);
            }
            bool isFind = false;
            for ( std::vector<std::string>::const_iterator iter_model2 = vec_model2.begin();
                  iter_model2 != vec_model2.end();
                  iter_model2++ )
            {
              std::string modellower2 = *iter_model2;
              for ( int index = 0; index < modellower2.length(); index++ )
              {
                modellower2[index] = tolower(modellower2[index]);
              }
              if ( modellower2.compare( modellower1 ) == 0 )
              {
                isFind = true;
              }
            }
            if ( !isFind )
            {
              isAllFind = false;
              break;
            }
          }
        }
        if ( isAllFind )
          isEqual = true;
      }
    }
    return isEqual;
  }

  std::string _brand;
  std::string _model;
  double _sensorSize;
};


#endif // DATASHEET_HPP

