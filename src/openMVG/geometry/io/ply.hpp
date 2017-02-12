#ifndef OPENMVG_GEOMETRY_IO_PLY_HPP
#define OPENMVG_GEOMETRY_IO_PLY_HPP

#include "openMVG/numeric/numeric.h"

#include <string>
#include <vector>

namespace openMVG
{
namespace geometry
{
  namespace io
  {
    /**
    * @brief Export a set of points to a PLY file 
    * @param pts A list of points 
    * @param path Path of the file to save 
    * @param binary Indicate if the file should be written to binary format 
    * @retval true if load is successful 
    * @retval false if load fails 
    */
    bool PLYWrite( const std::vector<Vec3> &pts,
                   const std::string &path,
                   const bool binary = true );

    /**
    * @brief Export a set of points to a PLY file 
    * @param pts A list of points 
    * @param nor A list of normal (one for each points) 
    * @param col A list of color (one for each points) 
    * @param path Path of the file to save 
    * @param binary Indicate if the file should be written to binary format 
    * @retval true if load is successful 
    * @retval false if load fails 
    */
    bool PLYWrite( const std::vector<Vec3> &pts,
                   const std::vector<Vec3> *nor,
                   const std::vector<Vec3uc> *col,
                   const std::string &path,
                   const bool binary = true );

    /**
    * @brief Load a ply file and gets is content points 
    * @param path Path of the file to load 
    * @param[out] pts The list of points in the given file 
    * @retval true if load is successful 
    * @retval false if load fails 
    */
    bool PLYRead( const std::string &path,
                  std::vector<Vec3> &pts );

    /**
    * @brief Load a ply file and gets is content points 
    * @param path Path of the file to load 
    * @param[out] pts The list of points in the given file 
    * @param[out] nor The list of normals in the given file 
    * @param[out] col The list of colors in the given file 
    * @note If nor or col are equal to nullptr, nothing is retrived 
    * @note If file does not contains nor or col, corresponding vectors will be empty  
    * @retval true if load is successful 
    * @retval false if load fails 
    */
    bool PLYRead( const std::string &path,
                  std::vector<Vec3> &pts,
                  std::vector<Vec3> *nor,
                  std::vector<Vec3uc> *col );

  } // namespace io
} // namespace geometry
} // namespace openMVG

#endif