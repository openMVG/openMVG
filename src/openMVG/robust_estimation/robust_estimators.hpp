#pragma once

namespace openMVG {
namespace robust {

enum EROBUST_ESTIMATOR
{
  ROBUST_ESTIMATOR_START = 0,
  ROBUST_ESTIMATOR_ACRANSAC = 1,        //< A-Contrario Ransac.
  ROBUST_ESTIMATOR_RANSAC = 2,          //< Classic Ransac.
  ROBUST_ESTIMATOR_LSMEDS = 3,          //< Variant of RANSAC using Least Median of Squares.
  ROBUST_ESTIMATOR_LORANSAC = 4,        //< LO-Ransac.
  ROBUST_ESTIMATOR_MAXCONSENSUS = 5,    //< Naive implementation of RANSAC without noise and iteration reduction options.
  ROBUST_ESTIMATOR_END
};

inline std::string EROBUST_ESTIMATOR_enumToString(EROBUST_ESTIMATOR estimator)
{
  switch(estimator)
  {
    case ROBUST_ESTIMATOR_ACRANSAC:
      return "acransac";
    case ROBUST_ESTIMATOR_RANSAC:
      return "ransac";
    case ROBUST_ESTIMATOR_LSMEDS:
      return "lsmeds";
    case ROBUST_ESTIMATOR_LORANSAC:
      return "loransac";
    case ROBUST_ESTIMATOR_MAXCONSENSUS:
      return "maxconsensus";
    case ROBUST_ESTIMATOR_START:
    case ROBUST_ESTIMATOR_END:
      break;
  }
  throw std::out_of_range("Invalid Ransac type Enum");
}

inline EROBUST_ESTIMATOR EROBUST_ESTIMATOR_stringToEnum(const std::string& estimator)
{
  if(estimator == "acransac")
    return ROBUST_ESTIMATOR_ACRANSAC;
  if(estimator == "ransac")
    return ROBUST_ESTIMATOR_RANSAC;
  if(estimator == "lsmeds")
    return ROBUST_ESTIMATOR_LSMEDS;
  if(estimator == "loransac")
    return ROBUST_ESTIMATOR_LORANSAC;
  if(estimator == "maxconsensus")
    return ROBUST_ESTIMATOR_MAXCONSENSUS;
  throw std::out_of_range("Invalid Ransac type string " + estimator);
}

} //namespace robust
} //namespace openMVG
