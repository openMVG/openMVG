// Copyright (c) 2013 libmv authors.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//
// Author: mierle@gmail.com (Keir Mierle)
//         sergey.vfx@gmail.com (Sergey Sharybin)
//
// This is an example application which contains bundle adjustment code used
// in the Libmv library and Blender. It reads problems from files passed via
// the command line and runs the bundle adjuster on the problem.
//
// File with problem a binary file, for which it is crucial to know in which
// order bytes of float values are stored in. This information is provided
// by a single character in the beginning of the file. There're two possible
// values of this byte:
//  - V, which means values in the file are stored with big endian type
//  - v, which means values in the file are stored with little endian type
//
// The rest of the file contains data in the following order:
//   - Space in which markers' coordinates are stored in
//   - Camera intrinsics
//   - Number of cameras
//   - Cameras
//   - Number of 3D points
//   - 3D points
//   - Number of markers
//   - Markers
//
// Markers' space could either be normalized or image (pixels). This is defined
// by the single character in the file. P means markers in the file is in image
// space, and N means markers are in normalized space.
//
// Camera intrinsics are 8 described by 8 float 8.
// This values goes in the following order:
//
//   - Focal length, principal point X, principal point Y, k1, k2, k3, p1, p2
//
// Every camera is described by:
//
//   - Image for which camera belongs to (single 4 bytes integer value).
//   - Column-major camera rotation matrix, 9 float values.
//   - Camera translation, 3-component vector of float values.
//
// Image number shall be greater or equal to zero. Order of cameras does not
// matter and gaps are possible.
//
// Every 3D point is decribed by:
//
//  - Track number point belongs to (single 4 bytes integer value).
//  - 3D position vector, 3-component vector of float values.
//
// Track number shall be greater or equal to zero. Order of tracks does not
// matter and gaps are possible.
//
// Finally every marker is described by:
//
//  - Image marker belongs to single 4 bytes integer value).
//  - Track marker belongs to single 4 bytes integer value).
//  - 2D marker position vector, (two float values).
//
// Marker's space is used a default value for refine_intrinsics command line
// flag. This means if there's no refine_intrinsics flag passed via command
// line, camera intrinsics will be refined if markers in the problem are
// stored in image space and camera intrinsics will not be refined if markers
// are in normalized space.
//
// Passing refine_intrinsics command line flag defines explicitly whether
// refinement of intrinsics will happen. Currently, only none and all
// intrinsics refinement is supported.
//
// There're existing problem files dumped from blender stored in folder
// ../data/libmv-ba-problems.

#include <cstdio>
#include <fcntl.h>
#include <sstream>
#include <string>
#include <vector>

#ifdef _MSC_VER
#  include <io.h>
#  define open _open
#  define close _close
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#include <unistd.h>

// O_BINARY is not defined on unix like platforms, as there is no
// difference between binary and text files.
#define O_BINARY 0

#endif

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "gflags/gflags.h"
#include "glog/logging.h"

typedef Eigen::Matrix<double, 3, 3> Mat3;
typedef Eigen::Matrix<double, 6, 1> Vec6;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector4d Vec4;

using std::vector;

DEFINE_string(input, "", "Input File name");
DEFINE_string(refine_intrinsics, "", "Camera intrinsics to be refined. "
              "Options are: none, radial.");

namespace {

// A EuclideanCamera is the location and rotation of the camera
// viewing an image.
//
// image identifies which image this camera represents.
// R is a 3x3 matrix representing the rotation of the camera.
// t is a translation vector representing its positions.
struct EuclideanCamera {
  EuclideanCamera() : image(-1) {}
  EuclideanCamera(const EuclideanCamera &c) : image(c.image), R(c.R), t(c.t) {}

  int image;
  Mat3 R;
  Vec3 t;
};

// A Point is the 3D location of a track.
//
// track identifies which track this point corresponds to.
// X represents the 3D position of the track.
struct EuclideanPoint {
  EuclideanPoint() : track(-1) {}
  EuclideanPoint(const EuclideanPoint &p) : track(p.track), X(p.X) {}
  int track;
  Vec3 X;
};

// A Marker is the 2D location of a tracked point in an image.
//
// x and y is the position of the marker in pixels from the top left corner
// in the image identified by an image. All markers for to the same target
// form a track identified by a common track number.
struct Marker {
  int image;
  int track;
  double x, y;
};

// Cameras intrinsics to be bundled.
//
// BUNDLE_RADIAL actually implies bundling of k1 and k2 coefficients only,
// no bundling of k3 is possible at this moment.
enum BundleIntrinsics {
  BUNDLE_NO_INTRINSICS = 0,
  BUNDLE_FOCAL_LENGTH = 1,
  BUNDLE_PRINCIPAL_POINT = 2,
  BUNDLE_RADIAL_K1 = 4,
  BUNDLE_RADIAL_K2 = 8,
  BUNDLE_RADIAL = 12,
  BUNDLE_TANGENTIAL_P1 = 16,
  BUNDLE_TANGENTIAL_P2 = 32,
  BUNDLE_TANGENTIAL = 48,
};

// Denotes which blocks to keep constant during bundling.
// For example it is useful to keep camera translations constant
// when bundling tripod motions.
enum BundleConstraints {
  BUNDLE_NO_CONSTRAINTS = 0,
  BUNDLE_NO_TRANSLATION = 1,
};

// The intrinsics need to get combined into a single parameter block; use these
// enums to index instead of numeric constants.
enum {
  OFFSET_FOCAL_LENGTH,
  OFFSET_PRINCIPAL_POINT_X,
  OFFSET_PRINCIPAL_POINT_Y,
  OFFSET_K1,
  OFFSET_K2,
  OFFSET_K3,
  OFFSET_P1,
  OFFSET_P2,
};

// Returns a pointer to the camera corresponding to a image.
EuclideanCamera *CameraForImage(vector<EuclideanCamera> *all_cameras,
                                const int image) {
  if (image < 0 || image >= all_cameras->size()) {
    return NULL;
  }
  EuclideanCamera *camera = &(*all_cameras)[image];
  if (camera->image == -1) {
    return NULL;
  }
  return camera;
}

const EuclideanCamera *CameraForImage(
    const vector<EuclideanCamera> &all_cameras,
    const int image) {
  if (image < 0 || image >= all_cameras.size()) {
    return NULL;
  }
  const EuclideanCamera *camera = &all_cameras[image];
  if (camera->image == -1) {
    return NULL;
  }
  return camera;
}

// Returns maximal image number at which marker exists.
int MaxImage(const vector<Marker> &all_markers) {
  if (all_markers.size() == 0) {
    return -1;
  }

  int max_image = all_markers[0].image;
  for (int i = 1; i < all_markers.size(); i++) {
    max_image = std::max(max_image, all_markers[i].image);
  }
  return max_image;
}

// Returns a pointer to the point corresponding to a track.
EuclideanPoint *PointForTrack(vector<EuclideanPoint> *all_points,
                              const int track) {
  if (track < 0 || track >= all_points->size()) {
    return NULL;
  }
  EuclideanPoint *point = &(*all_points)[track];
  if (point->track == -1) {
    return NULL;
  }
  return point;
}

// Reader of binary file which makes sure possibly needed endian
// conversion happens when loading values like floats and integers.
//
// File's endian type is reading from a first character of file, which
// could either be V for big endian or v for little endian.  This
// means you need to design file format assuming first character
// denotes file endianness in this way.
class EndianAwareFileReader {
 public:
  EndianAwareFileReader(void) : file_descriptor_(-1) {
    // Get an endian type of the host machine.
    union {
      unsigned char bytes[4];
      uint32_t value;
    } endian_test = { { 0, 1, 2, 3 } };
    host_endian_type_ = endian_test.value;
    file_endian_type_ = host_endian_type_;
  }

  ~EndianAwareFileReader(void) {
    if (file_descriptor_ > 0) {
      close(file_descriptor_);
    }
  }

  bool OpenFile(const std::string &file_name) {
    file_descriptor_ = open(file_name.c_str(), O_RDONLY | O_BINARY);
    if (file_descriptor_ < 0) {
      return false;
    }
    // Get an endian tpye of data in the file.
    unsigned char file_endian_type_flag = Read<unsigned char>();
    if (file_endian_type_flag == 'V') {
      file_endian_type_ = kBigEndian;
    } else if (file_endian_type_flag == 'v') {
      file_endian_type_ = kLittleEndian;
    } else {
      LOG(FATAL) << "Problem file is stored in unknown endian type.";
    }
    return true;
  }

  // Read value from the file, will switch endian if needed.
  template <typename T>
  T Read(void) const {
    T value;
    CHECK_GT(read(file_descriptor_, &value, sizeof(value)), 0);
    // Switch endian type if file contains data in different type
    // that current machine.
    if (file_endian_type_ != host_endian_type_) {
      value = SwitchEndian<T>(value);
    }
    return value;
  }
 private:
  static const long int kLittleEndian = 0x03020100ul;
  static const long int kBigEndian = 0x00010203ul;

  // Switch endian type between big to little.
  template <typename T>
  T SwitchEndian(const T value) const {
    if (sizeof(T) == 4) {
      unsigned int temp_value = static_cast<unsigned int>(value);
      return ((temp_value >> 24)) |
             ((temp_value << 8) & 0x00ff0000) |
             ((temp_value >> 8) & 0x0000ff00) |
             ((temp_value << 24));
    } else if (sizeof(T) == 1) {
      return value;
    } else {
      LOG(FATAL) << "Entered non-implemented part of endian switching function.";
    }
  }

  int host_endian_type_;
  int file_endian_type_;
  int file_descriptor_;
};

// Read 3x3 column-major matrix from the file
void ReadMatrix3x3(const EndianAwareFileReader &file_reader,
                   Mat3 *matrix) {
  for (int i = 0; i < 9; i++) {
    (*matrix)(i % 3, i / 3) = file_reader.Read<float>();
  }
}

// Read 3-vector from file
void ReadVector3(const EndianAwareFileReader &file_reader,
                 Vec3 *vector) {
  for (int i = 0; i < 3; i++) {
    (*vector)(i) = file_reader.Read<float>();
  }
}

// Reads a bundle adjustment problem from the file.
//
// file_name denotes from which file to read the problem.
// camera_intrinsics will contain initial camera intrinsics values.
//
// all_cameras is a vector of all reconstructed cameras to be optimized,
// vector element with number i will contain camera for image i.
//
// all_points is a vector of all reconstructed 3D points to be optimized,
// vector element with number i will contain point for track i.
//
// all_markers is a vector of all tracked markers existing in
// the problem. Only used for reprojection error calculation, stay
// unchanged during optimization.
//
// Returns false if any kind of error happened during
// reading.
bool ReadProblemFromFile(const std::string &file_name,
                         double camera_intrinsics[8],
                         vector<EuclideanCamera> *all_cameras,
                         vector<EuclideanPoint> *all_points,
                         bool *is_image_space,
                         vector<Marker> *all_markers) {
  EndianAwareFileReader file_reader;
  if (!file_reader.OpenFile(file_name)) {
    return false;
  }

  // Read markers' space flag.
  unsigned char is_image_space_flag = file_reader.Read<unsigned char>();
  if (is_image_space_flag == 'P') {
    *is_image_space = true;
  } else if (is_image_space_flag == 'N') {
    *is_image_space = false;
  } else {
    LOG(FATAL) << "Problem file contains markers stored in unknown space.";
  }

  // Read camera intrinsics.
  for (int i = 0; i < 8; i++) {
    camera_intrinsics[i] = file_reader.Read<float>();
  }

  // Read all cameras.
  int number_of_cameras = file_reader.Read<int>();
  for (int i = 0; i < number_of_cameras; i++) {
    EuclideanCamera camera;

    camera.image = file_reader.Read<int>();
    ReadMatrix3x3(file_reader, &camera.R);
    ReadVector3(file_reader, &camera.t);

    if (camera.image >= all_cameras->size()) {
      all_cameras->resize(camera.image + 1);
    }

    (*all_cameras)[camera.image].image = camera.image;
    (*all_cameras)[camera.image].R = camera.R;
    (*all_cameras)[camera.image].t = camera.t;
  }

  LOG(INFO) << "Read " << number_of_cameras << " cameras.";

  // Read all reconstructed 3D points.
  int number_of_points = file_reader.Read<int>();
  for (int i = 0; i < number_of_points; i++) {
    EuclideanPoint point;

    point.track = file_reader.Read<int>();
    ReadVector3(file_reader, &point.X);

    if (point.track >= all_points->size()) {
      all_points->resize(point.track + 1);
    }

    (*all_points)[point.track].track = point.track;
    (*all_points)[point.track].X = point.X;
  }

  LOG(INFO) << "Read " << number_of_points << " points.";

  // And finally read all markers.
  int number_of_markers = file_reader.Read<int>();
  for (int i = 0; i < number_of_markers; i++) {
    Marker marker;

    marker.image = file_reader.Read<int>();
    marker.track = file_reader.Read<int>();
    marker.x = file_reader.Read<float>();
    marker.y = file_reader.Read<float>();

    all_markers->push_back(marker);
  }

  LOG(INFO) << "Read " << number_of_markers << " markers.";

  return true;
}

// Apply camera intrinsics to the normalized point to get image coordinates.
// This applies the radial lens distortion to a point which is in normalized
// camera coordinates (i.e. the principal point is at (0, 0)) to get image
// coordinates in pixels. Templated for use with autodifferentiation.
template <typename T>
inline void ApplyRadialDistortionCameraIntrinsics(const T &focal_length_x,
                                                  const T &focal_length_y,
                                                  const T &principal_point_x,
                                                  const T &principal_point_y,
                                                  const T &k1,
                                                  const T &k2,
                                                  const T &k3,
                                                  const T &p1,
                                                  const T &p2,
                                                  const T &normalized_x,
                                                  const T &normalized_y,
                                                  T *image_x,
                                                  T *image_y) {
  T x = normalized_x;
  T y = normalized_y;

  // Apply distortion to the normalized points to get (xd, yd).
  T r2 = x*x + y*y;
  T r4 = r2 * r2;
  T r6 = r4 * r2;
  T r_coeff = 1.0 + k1 * r2 + k2 * r4 + k3 * r6;
  T xd = x * r_coeff + 2.0 * p1 * x * y + p2 * (r2 + 2.0 * x * x);
  T yd = y * r_coeff + 2.0 * p2 * x * y + p1 * (r2 + 2.0 * y * y);

  // Apply focal length and principal point to get the final image coordinates.
  *image_x = focal_length_x * xd + principal_point_x;
  *image_y = focal_length_y * yd + principal_point_y;
}

// Cost functor which computes reprojection error of 3D point X
// on camera defined by angle-axis rotation and it's translation
// (which are in the same block due to optimization reasons).
//
// This functor uses a radial distortion model.
struct OpenCVReprojectionError {
  OpenCVReprojectionError(const double observed_x, const double observed_y)
      : observed_x(observed_x), observed_y(observed_y) {}

  template <typename T>
  bool operator()(const T* const intrinsics,
                  const T* const R_t,  // Rotation denoted by angle axis
                                       // followed with translation
                  const T* const X,    // Point coordinates 3x1.
                  T* residuals) const {
    // Unpack the intrinsics.
    const T& focal_length      = intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = intrinsics[OFFSET_PRINCIPAL_POINT_Y];
    const T& k1                = intrinsics[OFFSET_K1];
    const T& k2                = intrinsics[OFFSET_K2];
    const T& k3                = intrinsics[OFFSET_K3];
    const T& p1                = intrinsics[OFFSET_P1];
    const T& p2                = intrinsics[OFFSET_P2];

    // Compute projective coordinates: x = RX + t.
    T x[3];

    ceres::AngleAxisRotatePoint(R_t, X, x);
    x[0] += R_t[3];
    x[1] += R_t[4];
    x[2] += R_t[5];

    // Compute normalized coordinates: x /= x[2].
    T xn = x[0] / x[2];
    T yn = x[1] / x[2];

    T predicted_x, predicted_y;

    // Apply distortion to the normalized points to get (xd, yd).
    // TODO(keir): Do early bailouts for zero distortion; these are expensive
    // jet operations.
    ApplyRadialDistortionCameraIntrinsics(focal_length,
                                          focal_length,
                                          principal_point_x,
                                          principal_point_y,
                                          k1, k2, k3,
                                          p1, p2,
                                          xn, yn,
                                          &predicted_x,
                                          &predicted_y);

    // The error is the difference between the predicted and observed position.
    residuals[0] = predicted_x - observed_x;
    residuals[1] = predicted_y - observed_y;

    return true;
  }

  const double observed_x;
  const double observed_y;
};

// Print a message to the log which camera intrinsics are gonna to be optimized.
void BundleIntrinsicsLogMessage(const int bundle_intrinsics) {
  if (bundle_intrinsics == BUNDLE_NO_INTRINSICS) {
    LOG(INFO) << "Bundling only camera positions.";
  } else {
    std::string bundling_message = "";

#define APPEND_BUNDLING_INTRINSICS(name, flag) \
    if (bundle_intrinsics & flag) { \
      if (!bundling_message.empty()) { \
        bundling_message += ", "; \
      } \
      bundling_message += name; \
    } (void)0

    APPEND_BUNDLING_INTRINSICS("f",      BUNDLE_FOCAL_LENGTH);
    APPEND_BUNDLING_INTRINSICS("px, py", BUNDLE_PRINCIPAL_POINT);
    APPEND_BUNDLING_INTRINSICS("k1",     BUNDLE_RADIAL_K1);
    APPEND_BUNDLING_INTRINSICS("k2",     BUNDLE_RADIAL_K2);
    APPEND_BUNDLING_INTRINSICS("p1",     BUNDLE_TANGENTIAL_P1);
    APPEND_BUNDLING_INTRINSICS("p2",     BUNDLE_TANGENTIAL_P2);

    LOG(INFO) << "Bundling " << bundling_message << ".";
  }
}

// Print a message to the log containing all the camera intriniscs values.
void PrintCameraIntrinsics(const char *text, const double *camera_intrinsics) {
  std::ostringstream intrinsics_output;

  intrinsics_output << "f=" << camera_intrinsics[OFFSET_FOCAL_LENGTH];

  intrinsics_output <<
    " cx=" << camera_intrinsics[OFFSET_PRINCIPAL_POINT_X] <<
    " cy=" << camera_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

#define APPEND_DISTORTION_COEFFICIENT(name, offset) \
  { \
    if (camera_intrinsics[offset] != 0.0) { \
      intrinsics_output << " " name "=" << camera_intrinsics[offset];  \
    } \
  } (void)0

  APPEND_DISTORTION_COEFFICIENT("k1", OFFSET_K1);
  APPEND_DISTORTION_COEFFICIENT("k2", OFFSET_K2);
  APPEND_DISTORTION_COEFFICIENT("k3", OFFSET_K3);
  APPEND_DISTORTION_COEFFICIENT("p1", OFFSET_P1);
  APPEND_DISTORTION_COEFFICIENT("p2", OFFSET_P2);

#undef APPEND_DISTORTION_COEFFICIENT

  LOG(INFO) << text << intrinsics_output.str();
}

// Get a vector of camera's rotations denoted by angle axis
// conjuncted with translations into single block
//
// Element with index i matches to a rotation+translation for
// camera at image i.
vector<Vec6> PackCamerasRotationAndTranslation(
    const vector<Marker> &all_markers,
    const vector<EuclideanCamera> &all_cameras) {
  vector<Vec6> all_cameras_R_t;
  int max_image = MaxImage(all_markers);

  all_cameras_R_t.resize(max_image + 1);

  for (int i = 0; i <= max_image; i++) {
    const EuclideanCamera *camera = CameraForImage(all_cameras, i);

    if (!camera) {
      continue;
    }

    ceres::RotationMatrixToAngleAxis(&camera->R(0, 0),
                                     &all_cameras_R_t[i](0));
    all_cameras_R_t[i].tail<3>() = camera->t;
  }

  return all_cameras_R_t;
}

// Convert cameras rotations fro mangle axis back to rotation matrix.
void UnpackCamerasRotationAndTranslation(
    const vector<Marker> &all_markers,
    const vector<Vec6> &all_cameras_R_t,
    vector<EuclideanCamera> *all_cameras) {
  int max_image = MaxImage(all_markers);

  for (int i = 0; i <= max_image; i++) {
    EuclideanCamera *camera = CameraForImage(all_cameras, i);

    if (!camera) {
      continue;
    }

    ceres::AngleAxisToRotationMatrix(&all_cameras_R_t[i](0),
                                     &camera->R(0, 0));
    camera->t = all_cameras_R_t[i].tail<3>();
  }
}

void EuclideanBundleCommonIntrinsics(const vector<Marker> &all_markers,
                                     const int bundle_intrinsics,
                                     const int bundle_constraints,
                                     double *camera_intrinsics,
                                     vector<EuclideanCamera> *all_cameras,
                                     vector<EuclideanPoint> *all_points) {
  PrintCameraIntrinsics("Original intrinsics: ", camera_intrinsics);

  ceres::Problem::Options problem_options;
  ceres::Problem problem(problem_options);

  // Convert cameras rotations to angle axis and merge with translation
  // into single parameter block for maximal minimization speed
  //
  // Block for minimization has got the following structure:
  //   <3 elements for angle-axis> <3 elements for translation>
  vector<Vec6> all_cameras_R_t =
    PackCamerasRotationAndTranslation(all_markers, *all_cameras);

  // Parameterization used to restrict camera motion for modal solvers.
  ceres::SubsetParameterization *constant_transform_parameterization = NULL;
  if (bundle_constraints & BUNDLE_NO_TRANSLATION) {
      std::vector<int> constant_translation;

      // First three elements are rotation, last three are translation.
      constant_translation.push_back(3);
      constant_translation.push_back(4);
      constant_translation.push_back(5);

      constant_transform_parameterization =
        new ceres::SubsetParameterization(6, constant_translation);
  }

  int num_residuals = 0;
  bool have_locked_camera = false;
  for (int i = 0; i < all_markers.size(); ++i) {
    const Marker &marker = all_markers[i];
    EuclideanCamera *camera = CameraForImage(all_cameras, marker.image);
    EuclideanPoint *point = PointForTrack(all_points, marker.track);
    if (camera == NULL || point == NULL) {
      continue;
    }

    // Rotation of camera denoted in angle axis followed with
    // camera translaiton.
    double *current_camera_R_t = &all_cameras_R_t[camera->image](0);

    problem.AddResidualBlock(new ceres::AutoDiffCostFunction<
        OpenCVReprojectionError, 2, 8, 6, 3>(
            new OpenCVReprojectionError(
                marker.x,
                marker.y)),
        NULL,
        camera_intrinsics,
        current_camera_R_t,
        &point->X(0));

    // We lock the first camera to better deal with scene orientation ambiguity.
    if (!have_locked_camera) {
      problem.SetParameterBlockConstant(current_camera_R_t);
      have_locked_camera = true;
    }

    if (bundle_constraints & BUNDLE_NO_TRANSLATION) {
      problem.SetParameterization(current_camera_R_t,
                                  constant_transform_parameterization);
    }

    num_residuals++;
  }
  LOG(INFO) << "Number of residuals: " << num_residuals;

  if (!num_residuals) {
    LOG(INFO) << "Skipping running minimizer with zero residuals";
    return;
  }

  BundleIntrinsicsLogMessage(bundle_intrinsics);

  if (bundle_intrinsics == BUNDLE_NO_INTRINSICS) {
    // No camera intrinsics are being refined,
    // set the whole parameter block as constant for best performance.
    problem.SetParameterBlockConstant(camera_intrinsics);
  } else {
    // Set the camera intrinsics that are not to be bundled as
    // constant using some macro trickery.

    std::vector<int> constant_intrinsics;
#define MAYBE_SET_CONSTANT(bundle_enum, offset) \
    if (!(bundle_intrinsics & bundle_enum)) { \
      constant_intrinsics.push_back(offset); \
    }
    MAYBE_SET_CONSTANT(BUNDLE_FOCAL_LENGTH,    OFFSET_FOCAL_LENGTH);
    MAYBE_SET_CONSTANT(BUNDLE_PRINCIPAL_POINT, OFFSET_PRINCIPAL_POINT_X);
    MAYBE_SET_CONSTANT(BUNDLE_PRINCIPAL_POINT, OFFSET_PRINCIPAL_POINT_Y);
    MAYBE_SET_CONSTANT(BUNDLE_RADIAL_K1,       OFFSET_K1);
    MAYBE_SET_CONSTANT(BUNDLE_RADIAL_K2,       OFFSET_K2);
    MAYBE_SET_CONSTANT(BUNDLE_TANGENTIAL_P1,   OFFSET_P1);
    MAYBE_SET_CONSTANT(BUNDLE_TANGENTIAL_P2,   OFFSET_P2);
#undef MAYBE_SET_CONSTANT

    // Always set K3 constant, it's not used at the moment.
    constant_intrinsics.push_back(OFFSET_K3);

    ceres::SubsetParameterization *subset_parameterization =
      new ceres::SubsetParameterization(8, constant_intrinsics);

    problem.SetParameterization(camera_intrinsics, subset_parameterization);
  }

  // Configure the solver.
  ceres::Solver::Options options;
  options.use_nonmonotonic_steps = true;
  options.preconditioner_type = ceres::SCHUR_JACOBI;
  options.linear_solver_type = ceres::ITERATIVE_SCHUR;
  options.use_inner_iterations = true;
  options.max_num_iterations = 100;
  options.minimizer_progress_to_stdout = true;

  // Solve!
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  std::cout << "Final report:\n" << summary.FullReport();

  // Copy rotations and translations back.
  UnpackCamerasRotationAndTranslation(all_markers,
                                      all_cameras_R_t,
                                      all_cameras);

  PrintCameraIntrinsics("Final intrinsics: ", camera_intrinsics);
}
}  // namespace

int main(int argc, char **argv) {
  CERES_GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  if (FLAGS_input.empty()) {
    LOG(ERROR) << "Usage: libmv_bundle_adjuster --input=blender_problem";
    return EXIT_FAILURE;
  }

  double camera_intrinsics[8];
  vector<EuclideanCamera> all_cameras;
  vector<EuclideanPoint> all_points;
  bool is_image_space;
  vector<Marker> all_markers;

  if (!ReadProblemFromFile(FLAGS_input,
                           camera_intrinsics,
                           &all_cameras,
                           &all_points,
                           &is_image_space,
                           &all_markers)) {
    LOG(ERROR) << "Error reading problem file";
    return EXIT_FAILURE;
  }

  // If there's no refine_intrinsics passed via command line
  // (in this case FLAGS_refine_intrinsics will be an empty string)
  // we use problem's settings to detect whether intrinsics
  // shall be refined or not.
  //
  // Namely, if problem has got markers stored in image (pixel)
  // space, we do full intrinsics refinement. If markers are
  // stored in normalized space, and refine_intrinsics is not
  // set, no refining will happen.
  //
  // Using command line argument refine_intrinsics will explicitly
  // declare which intrinsics need to be refined and in this case
  // refining flags does not depend on problem at all.
  int bundle_intrinsics = BUNDLE_NO_INTRINSICS;
  if (FLAGS_refine_intrinsics.empty()) {
    if (is_image_space) {
      bundle_intrinsics = BUNDLE_FOCAL_LENGTH | BUNDLE_RADIAL;
    }
  } else {
    if (FLAGS_refine_intrinsics == "radial") {
      bundle_intrinsics = BUNDLE_FOCAL_LENGTH | BUNDLE_RADIAL;
    } else if (FLAGS_refine_intrinsics != "none") {
      LOG(ERROR) << "Unsupported value for refine-intrinsics";
      return EXIT_FAILURE;
    }
  }

  // Run the bundler.
  EuclideanBundleCommonIntrinsics(all_markers,
                                  bundle_intrinsics,
                                  BUNDLE_NO_CONSTRAINTS,
                                  camera_intrinsics,
                                  &all_cameras,
                                  &all_points);

  return EXIT_SUCCESS;
}
