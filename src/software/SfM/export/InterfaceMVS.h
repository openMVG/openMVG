#ifndef _INTERFACE_MVS_H_
#define _INTERFACE_MVS_H_


// I N C L U D E S /////////////////////////////////////////////////

#include <fstream>
#include <string>
#include <cctype>
#include <limits>


// D E F I N E S ///////////////////////////////////////////////////

#define MVSI_PROJECT_ID "MVSI" // identifies the project stream
#define MVSI_PROJECT_VER ((uint32_t)6) // identifies the version of a project stream

// set a default namespace name if none given
#ifndef _INTERFACE_NAMESPACE
#define _INTERFACE_NAMESPACE MVS
#endif

// uncomment to enable custom OpenCV data types
// (should be uncommented if OpenCV is not available)
#if !defined(_USE_OPENCV) && !defined(_USE_CUSTOM_CV)
#define _USE_CUSTOM_CV
#endif

// set to disable custom NO_ID declaration
#ifndef _DISABLE_NO_ID
#define _INTERFACE_NO_ID
#endif


// S T R U C T S ///////////////////////////////////////////////////

#ifdef _USE_CUSTOM_CV

namespace cv {

// simple cv::Point3_
template<typename Type>
class Point3_
{
public:
	typedef Type value_type;

	inline Point3_() {}
	inline Point3_(Type _x, Type _y, Type _z) : x(_x), y(_y), z(_z) {}
	#ifdef _USE_EIGEN
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF_VECTORIZABLE_FIXED_SIZE(Type,3)
	typedef Eigen::Matrix<Type,3,1> EVec;
	typedef Eigen::Map<EVec> EVecMap;
	template<typename Derived>
	inline Point3_(const Eigen::EigenBase<Derived>& rhs) { operator EVecMap () = rhs; }
	template<typename Derived>
	inline Point3_& operator = (const Eigen::EigenBase<Derived>& rhs) { operator EVecMap () = rhs; return *this; }
	inline operator const EVecMap () const { return EVecMap((Type*)this); }
	inline operator EVecMap () { return EVecMap((Type*)this); }
	#endif

	const Type* ptr() const { return &x; }
	Type* ptr() { return &x; }
	Type operator()(int r) const { return (&x)[r]; }
	Type& operator()(int r) { return (&x)[r]; }
	Point3_ operator - () const {
		return Point3_(
			-x,
			-y,
			-z
		);
	}
	Point3_ operator + (const Point3_& X) const {
		return Point3_(
			x+X.x,
			y+X.y,
			z+X.z
		);
	}
	Point3_ operator - (const Point3_& X) const {
		return Point3_(
			x-X.x,
			y-X.y,
			z-X.z
		);
	}

public:
	Type x, y, z;
};

// simple cv::Matx
template<typename Type, int m, int n>
class Matx
{
public:
	typedef Type value_type;
	enum {
		rows     = m,
		cols     = n,
		channels = rows*cols
	};

	inline Matx() {}
	#ifdef _USE_EIGEN
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF_VECTORIZABLE_FIXED_SIZE(Type,m*n)
	typedef Eigen::Matrix<Type,m,n,(n>1?Eigen::RowMajor:Eigen::Default)> EMat;
	typedef Eigen::Map<const EMat> CEMatMap;
	typedef Eigen::Map<EMat> EMatMap;
	template<typename Derived>
	inline Matx(const Eigen::EigenBase<Derived>& rhs) { operator EMatMap () = rhs; }
	template<typename Derived>
	inline Matx& operator = (const Eigen::EigenBase<Derived>& rhs) { operator EMatMap () = rhs; return *this; }
	inline operator CEMatMap() const { return CEMatMap((const Type*)val); }
	inline operator EMatMap () { return EMatMap((Type*)val); }
	#endif

	Type operator()(int r, int c) const { return val[r*n+c]; }
	Type& operator()(int r, int c) { return val[r*n+c]; }
	Point3_<Type> operator * (const Point3_<Type>& X) const {
		Point3_<Type> R;
		for (int r = 0; r < m; r++) {
			R(r) = Type(0);
			for (int c = 0; c < n; c++)
				R(r) += operator()(r,c)*X(c);
		}
		return R;
	}
	template<int k>
	Matx<Type,m,k> operator * (const Matx<Type,n,k>& M) const {
		Matx<Type,m,k> R;
		for (int r = 0; r < m; r++) {
			for (int l = 0; l < k; l++) {
				R(r,l) = Type(0);
				for (int c = 0; c < n; c++)
					R(r,l) += operator()(r,c)*M(c,l);
			}
		}
		return R;
	}
	Matx<Type,n,m> t() const {
		Matx<Type,n,m> M;
		for (int r = 0; r < m; r++)
			for (int c = 0; c < n; c++)
				M(c,r) = operator()(r,c);
		return M;
	}

	static Matx eye() {
		Matx M;
		memset(M.val, 0, sizeof(Type)*m*n);
		const int shortdim(m < n ? m : n);
		for (int i = 0; i < shortdim; i++)
			M(i,i) = 1;
		return M;
	}

public:
	Type val[m*n];
};

} // namespace cv
#endif
/*----------------------------------------------------------------*/


namespace _INTERFACE_NAMESPACE {

// invalid index
#ifdef _INTERFACE_NO_ID
constexpr uint32_t NO_ID {std::numeric_limits<uint32_t>::max()};
#endif

// custom serialization
namespace ARCHIVE {

// Basic serialization types
struct ArchiveSave {
	std::ostream& stream;
	uint32_t version;
	ArchiveSave(std::ostream& _stream, uint32_t _version)
		: stream(_stream), version(_version) {}
	template<typename _Tp>
	ArchiveSave& operator & (const _Tp& obj);
};
struct ArchiveLoad {
	std::istream& stream;
	uint32_t version;
	ArchiveLoad(std::istream& _stream, uint32_t _version)
		: stream(_stream), version(_version) {}
	template<typename _Tp>
	ArchiveLoad& operator & (_Tp& obj);
};

template<typename _Tp>
bool Save(ArchiveSave& a, const _Tp& obj) {
	const_cast<_Tp&>(obj).serialize(a, a.version);
	return true;
}
template<typename _Tp>
bool Load(ArchiveLoad& a, _Tp& obj) {
	obj.serialize(a, a.version);
	return true;
}

template<typename _Tp>
ArchiveSave& ArchiveSave::operator & (const _Tp& obj) {
	Save(*this, obj);
	return *this;
}
template<typename _Tp>
ArchiveLoad& ArchiveLoad::operator & (_Tp& obj) {
	Load(*this, obj);
	return *this;
}

// Main exporter & importer
template<typename _Tp>
bool SerializeSave(const _Tp& obj, const std::string& fileName, uint32_t version=MVSI_PROJECT_VER) {
	// open the output stream
	std::ofstream stream(fileName, std::ofstream::binary);
	if (!stream.is_open())
		return false;
	// write header
	if (version > 0) {
		// save project ID
		stream.write(MVSI_PROJECT_ID, 4);
		// save project version
		stream.write((const char*)&version, sizeof(uint32_t));
		// reserve some bytes
		const uint32_t reserved(0);
		stream.write((const char*)&reserved, sizeof(uint32_t));
	}
	// serialize out the current state
	ARCHIVE::ArchiveSave serializer(stream, version);
	serializer & obj;
	return true;
}
template<typename _Tp>
bool SerializeLoad(_Tp& obj, const std::string& fileName, uint32_t* pVersion=NULL) {
	// open the input stream
	std::ifstream stream(fileName, std::ifstream::binary);
	if (!stream.is_open())
		return false;
	// read header
	uint32_t version(0);
	// load project header ID
	char szHeader[4];
	stream.read(szHeader, 4);
	if (!stream)
		return false;
	if (strncmp(szHeader, MVSI_PROJECT_ID, 4) != 0) {
		// try to load as the first version that didn't have a header
		const size_t size(fileName.size());
		if (size <= 4)
			return false;
		std::string ext(fileName.substr(size-4));
		std::transform(ext.begin(), ext.end(), ext.begin(), [](char c) { return (char)std::tolower(c); });
		if (ext != ".mvs")
			return false;
		stream.seekg(0, std::ifstream::beg);
	} else {
		// load project version
		stream.read((char*)&version, sizeof(uint32_t));
		if (!stream || version > MVSI_PROJECT_VER)
			return false;
		// skip reserved bytes
		uint32_t reserved;
		stream.read((char*)&reserved, sizeof(uint32_t));
	}
	// serialize in the current state
	ARCHIVE::ArchiveLoad serializer(stream, version);
	serializer & obj;
	if (pVersion)
		*pVersion = version;
	return true;
}


#define ARCHIVE_DEFINE_TYPE(TYPE) \
template<> \
inline bool Save<TYPE>(ArchiveSave& a, const TYPE& v) { \
	a.stream.write((const char*)&v, sizeof(TYPE)); \
	return true; \
} \
template<> \
inline bool Load<TYPE>(ArchiveLoad& a, TYPE& v) { \
	a.stream.read((char*)&v, sizeof(TYPE)); \
	return true; \
}

// Serialization support for basic types
ARCHIVE_DEFINE_TYPE(uint32_t)
ARCHIVE_DEFINE_TYPE(uint64_t)
ARCHIVE_DEFINE_TYPE(float)
ARCHIVE_DEFINE_TYPE(double)

// Serialization support for cv::Matx
template<typename _Tp, int m, int n>
inline bool Save(ArchiveSave& a, const cv::Matx<_Tp,m,n>& _m) {
	a.stream.write((const char*)_m.val, sizeof(_Tp)*m*n);
	return true;
}
template<typename _Tp, int m, int n>
inline bool Load(ArchiveLoad& a, cv::Matx<_Tp,m,n>& _m) {
	a.stream.read((char*)_m.val, sizeof(_Tp)*m*n);
	return true;
}

// Serialization support for cv::Point3_
template<typename _Tp>
inline bool Save(ArchiveSave& a, const cv::Point3_<_Tp>& pt) {
	a.stream.write((const char*)&pt.x, sizeof(_Tp)*3);
	return true;
}
template<typename _Tp>
inline bool Load(ArchiveLoad& a, cv::Point3_<_Tp>& pt) {
	a.stream.read((char*)&pt.x, sizeof(_Tp)*3);
	return true;
}

// Serialization support for std::string
template<>
inline bool Save<std::string>(ArchiveSave& a, const std::string& s) {
	const uint64_t size(s.size());
	Save(a, size);
	if (size > 0)
		a.stream.write(&s[0], sizeof(char)*size);
	return true;
}
template<>
inline bool Load<std::string>(ArchiveLoad& a, std::string& s) {
	uint64_t size;
	Load(a, size);
	if (size > 0) {
		s.resize(size);
		a.stream.read(&s[0], sizeof(char)*size);
	}
	return true;
}

// Serialization support for std::vector
template<typename _Tp>
inline bool Save(ArchiveSave& a, const std::vector<_Tp>& v) {
	const uint64_t size(v.size());
	Save(a, size);
	for (uint64_t i=0; i<size; ++i)
		Save(a, v[i]);
	return true;
}
template<typename _Tp>
inline bool Load(ArchiveLoad& a, std::vector<_Tp>& v) {
	uint64_t size;
	Load(a, size);
	if (size > 0) {
		v.resize(size);
		for (uint64_t i=0; i<size; ++i)
			Load(a, v[i]);
	}
	return true;
}

} // namespace ARCHIVE
/*----------------------------------------------------------------*/


// interface used to export/import MVS input data;
//  - MAX(width,height) is used for normalization
//  - row-major order is used for storing the matrices
struct Interface
{
	typedef cv::Point3_<float> Pos3f;
	typedef cv::Point3_<double> Pos3d;
	typedef cv::Matx<double,3,3> Mat33d;
	typedef cv::Matx<double,4,4> Mat44d;
	typedef cv::Point3_<uint8_t> Col3; // x=B, y=G, z=R
	/*----------------------------------------------------------------*/

	// structure describing a mobile platform with cameras attached to it
	struct Platform {
		// structure describing a camera mounted on a platform
		struct Camera {
			std::string name; // camera's name
			std::string bandName; // camera's band name, ex: RGB, BLUE, GREEN, RED, NIR, THERMAL, etc (optional)
			uint32_t width, height; // image resolution in pixels for all images sharing this camera (optional)
			Mat33d K; // camera's intrinsics matrix (normalized if image resolution not specified), where integer coordinates is by convention the pixel center
			Mat33d R; // camera's rotation matrix relative to the platform
			Pos3d C; // camera's translation vector relative to the platform

			Camera() : width(0), height(0) {}
			bool HasResolution() const { return width > 0 && height > 0; }
			bool IsNormalized() const { return !HasResolution(); }
			static uint32_t GetNormalizationScale(uint32_t width, uint32_t height) { return std::max(width, height); }
			uint32_t GetNormalizationScale() const { return GetNormalizationScale(width, height); }

			// project point: camera to image (homogeneous) coordinates
			inline Pos3d operator * (const Pos3d& X) const {
				return Pos3d(
					K(0,2)+K(0,0)*X.x/X.z,
					K(1,2)+K(1,1)*X.y/X.z,
					1.0);
			}
			// back-project point: image (z is the depth) to camera coordinates
			inline Pos3d operator / (const Pos3d& x) const {
				return Pos3d(
					(x.x-K(0,2))*x.z/K(0,0),
					(x.y-K(1,2))*x.z/K(1,1),
					1.0);
			}

			template <class Archive>
			void serialize(Archive& ar, const unsigned int version) {
				ar & name;
				if (version > 3) {
					ar & bandName;
				}
				if (version > 0) {
					ar & width;
					ar & height;
				}
				ar & K;
				ar & R;
				ar & C;
			}
		};
		typedef std::vector<Camera> CameraArr;

		// structure describing a pose along the trajectory of a platform
		struct Pose {
			Mat33d R; // platform's rotation matrix that rotates a point from world to camera coordinate system
			Pos3d C; // platform's translation vector (position) in world coordinate system

			Pose() {}
			template <typename MAT, typename POS>
			Pose(const MAT& _R, const POS& _C) : R(_R), C(_C) {}

			// translation vector t = -RC
			inline Pos3d GetTranslation() const { return R*(-C); }
			inline void SetTranslation(const Pos3d& T) { C = R.t()*(-T); }

			// combine poses
			inline Pose operator * (const Pose& P) const {
				return Pose(R*P.R, P.R.t()*C+P.C);
			}
			inline Pose& operator *= (const Pose& P) {
				R = R*P.R; C = P.R.t()*C+P.C; return *this;
			}

			// project point: world to local coordinates
			inline Pos3d operator * (const Pos3d& X) const {
				return R * (X - C);
			}
			// back-project point: local to world coordinates
			inline Pos3d operator / (const Pos3d& X) const {
				return R.t() * X + C;
			}

			template <class Archive>
			void serialize(Archive& ar, const unsigned int /*version*/) {
				ar & R;
				ar & C;
			}
		};
		typedef std::vector<Pose> PoseArr;

		std::string name; // platform's name
		CameraArr cameras; // cameras mounted on the platform
		PoseArr poses; // trajectory of the platform

		const Mat33d& GetK(uint32_t cameraID) const {
			return cameras[cameraID].K;
		}
		static Mat33d ScaleK(const Mat33d& _K, double scale) {
			Mat33d K(_K);
			const bool bNormalized(K(0,2) < 3 && K(1,2) < 3);
			K(0,0) *= scale;
			K(1,1) *= scale;
			K(0,2) = bNormalized ? K(0,2)*scale : (K(0,2)+0.5)*scale-0.5;
			K(1,2) = bNormalized ? K(1,2)*scale : (K(1,2)+0.5)*scale-0.5;
			K(0,1) *= scale;
			return K;
		}
		const Mat33d& SetFullK(uint32_t cameraID, const Mat33d& K, uint32_t width, uint32_t height, bool normalize=false) {
			Camera& camera = cameras[cameraID];
			if (normalize) {
				camera.width = camera.height = 0;
				camera.K = ScaleK(K, 1.0/(double)Camera::GetNormalizationScale(width, height));
			} else {
				camera.width = width; camera.height = height;
				camera.K = K;
			}
			return camera.K;
		}
		Mat33d GetFullK(uint32_t cameraID, uint32_t width, uint32_t height) const {
			const Camera& camera = cameras[cameraID];
			if (!camera.IsNormalized() && camera.width == width && camera.height == height)
				return camera.K;
			return ScaleK(camera.K, (double)Camera::GetNormalizationScale(width, height)/
				(camera.IsNormalized()?1.0:(double)camera.GetNormalizationScale()));
		}

		Pose GetPose(uint32_t cameraID, uint32_t poseID) const {
			const Camera& camera = cameras[cameraID];
			const Pose& pose = poses[poseID];
			// add the relative camera pose to the platform
			return Pose{
				camera.R*pose.R,
				pose.R.t()*camera.C+pose.C
			};
		}

		template <class Archive>
		void serialize(Archive& ar, const unsigned int /*version*/) {
			ar & name;
			ar & cameras;
			ar & poses;
		}
	};
	typedef std::vector<Platform> PlatformArr;
	/*----------------------------------------------------------------*/

	// structure describing an image
	struct Image {
		// structure describing how an other image relates to this image in terms of overlap,
		// i.e. how many 3D points are shared between the two images, base-line and common area,
		// useful for ex. when selecting the best images to densly match with
		struct ViewScore {
			uint32_t ID; // image local-ID, the index in this scene images list
			uint32_t points; // number of 3D points shared with the reference image
			float scale; // image scale relative to the reference image
			float angle; // image angle relative to the reference image (radians)
			float area; // common image area relative to the reference image (ratio)
			float score; // aggregated image score relative to the reference image (larger is better)

			template<class Archive>
			void serialize(Archive& ar, const unsigned int /*version*/) {
				ar & ID;
				ar & points;
				ar & scale;
				ar & angle;
				ar & area;
				ar & score;
			}
		};
		
		std::string name; // image file name
		std::string maskName; // segmentation file name (optional)
		uint32_t platformID; // ID of the associated platform
		uint32_t cameraID; // ID of the associated camera on the associated platform
		uint32_t poseID; // ID of the pose of the associated platform
		uint32_t ID; // image global-ID, ex. the ID given outside the current scene, like the index in the full list of image files (optional)
		float minDepth; // minimum depth of the points seen by this image (optional)
		float avgDepth; // average depth of the points seen by this image (optional)
		float maxDepth; // maximum depth of the points seen by this image (optional)
		std::vector<ViewScore> viewScores; // list of view scores for this image (optional)

		Image() : platformID(NO_ID), cameraID(NO_ID), poseID(NO_ID), ID(NO_ID), minDepth(0), avgDepth(0), maxDepth(0) {}

		bool IsValid() const { return poseID != NO_ID; }

		template <class Archive>
		void serialize(Archive& ar, const unsigned int version) {
			ar & name;
			if (version > 4) {
				ar & maskName;
			}
			ar & platformID;
			ar & cameraID;
			ar & poseID;
			if (version > 2) {
				ar & ID;
			}
			if (version > 6) {
				ar & minDepth;
				ar & avgDepth;
				ar & maxDepth;
				ar & viewScores;
			}
		}
	};
	typedef std::vector<Image> ImageArr;
	/*----------------------------------------------------------------*/

	// structure describing a 3D point
	struct Vertex {
		// structure describing one view for a given 3D feature
		struct View {
			uint32_t imageID; // image ID corresponding to this view
			float confidence; // view's confidence (0 - not available)

			template<class Archive>
			void serialize(Archive& ar, const unsigned int /*version*/) {
				ar & imageID;
				ar & confidence;
			}
		};
		typedef std::vector<View> ViewArr;

		Pos3f X; // 3D point position
		ViewArr views; // list of all available views for this 3D feature

		template <class Archive>
		void serialize(Archive& ar, const unsigned int /*version*/) {
			ar & X;
			ar & views;
		}
	};
	typedef std::vector<Vertex> VertexArr;
	/*----------------------------------------------------------------*/

	// structure describing a 3D line
	struct Line {
		// structure describing one view for a given 3D feature
		struct View {
			uint32_t imageID; // image ID corresponding to this view
			float confidence; // view's confidence (0 - not available)

			template<class Archive>
			void serialize(Archive& ar, const unsigned int /*version*/) {
				ar & imageID;
				ar & confidence;
			}
		};
		typedef std::vector<View> ViewArr;

		Pos3f pt1; // 3D line segment end-point
		Pos3f pt2; // 3D line segment end-point
		ViewArr views; // list of all available views for this 3D feature

		template <class Archive>
		void serialize(Archive& ar, const unsigned int /*version*/) {
			ar & pt1;
			ar & pt2;
			ar & views;
		}
	};
	typedef std::vector<Line> LineArr;
	/*----------------------------------------------------------------*/

	// structure describing a 3D point's normal (optional)
	struct Normal {
		Pos3f n; // 3D feature normal

		template <class Archive>
		void serialize(Archive& ar, const unsigned int /*version*/) {
			ar & n;
		}
	};
	typedef std::vector<Normal> NormalArr;
	/*----------------------------------------------------------------*/

	// structure describing a 3D point's color (optional)
	struct Color {
		Col3 c; // 3D feature color

		template <class Archive>
		void serialize(Archive& ar, const unsigned int /*version*/) {
			ar & c;
		}
	};
	typedef std::vector<Color> ColorArr;
	/*----------------------------------------------------------------*/

	// structure describing a Oriented Bounding-Box (optional)
	struct OBB {
		Mat33d rot; // rotation from scene to OBB coordinate system
		Pos3d ptMin; // minimal point represented in OBB coordinate system
		Pos3d ptMax; // maximal point represented in OBB coordinate system

		OBB() : rot(Mat33d::eye()), ptMin(0, 0, 0), ptMax(0, 0, 0) {}

		bool IsValid() const { return ptMin.x < ptMax.x && ptMin.y < ptMax.y && ptMin.z < ptMax.z; }

		template <class Archive>
		void serialize(Archive& ar, const unsigned int /*version*/) {
			ar & rot;
			ar & ptMin;
			ar & ptMax;
		}
	};
	/*----------------------------------------------------------------*/

	PlatformArr platforms; // array of platforms
	ImageArr images; // array of images
	VertexArr vertices; // array of reconstructed 3D points
	NormalArr verticesNormal; // array of reconstructed 3D points' normal (optional)
	ColorArr verticesColor; // array of reconstructed 3D points' color (optional)
	LineArr lines; // array of reconstructed 3D lines (optional)
	NormalArr linesNormal; // array of reconstructed 3D lines' normal (optional)
	ColorArr linesColor; // array of reconstructed 3D lines' color (optional)
	Mat44d transform; // transformation used to convert from absolute to relative coordinate system (optional)
	OBB obb; // minimum oriented bounding box containing the scene (optional)

	Interface() : transform(Mat44d::eye()) {}

	const Mat33d& GetK(uint32_t imageID) const {
		const Image& image = images[imageID];
		return platforms[image.platformID].GetK(image.cameraID);
	}
	Mat33d GetFullK(uint32_t imageID, uint32_t width, uint32_t height) const {
		const Image& image = images[imageID];
		return platforms[image.platformID].GetFullK(image.cameraID, width, height);
	}

	const Platform::Camera& GetCamera(uint32_t imageID) const {
		const Image& image = images[imageID];
		return platforms[image.platformID].cameras[image.cameraID];
	}

	Platform::Pose GetPose(uint32_t imageID) const {
		const Image& image = images[imageID];
		return platforms[image.platformID].GetPose(image.cameraID, image.poseID);
	}

	template <class Archive>
	void serialize(Archive& ar, const unsigned int version) {
		ar & platforms;
		ar & images;
		ar & vertices;
		ar & verticesNormal;
		ar & verticesColor;
		if (version > 0) {
			ar & lines;
			ar & linesNormal;
			ar & linesColor;
			if (version > 1) {
				ar & transform;
				if (version > 5) {
					ar & obb;
				}
			}
		}
	}
};
/*----------------------------------------------------------------*/

} // namespace _INTERFACE_NAMESPACE

#endif // _INTERFACE_MVS_H_
