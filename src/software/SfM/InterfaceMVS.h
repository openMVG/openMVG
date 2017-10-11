#ifndef _INTERFACE_MVS_H_
#define _INTERFACE_MVS_H_


// I N C L U D E S /////////////////////////////////////////////////

#include <fstream>

// D E F I N E S ///////////////////////////////////////////////////

#define MVSI_PROJECT_ID "MVSI" // identifies the project stream
#define MVSI_PROJECT_VER ((uint32_t)2) // identifies the version of a project stream

// set a default namespace name is none given
#ifndef _INTERFACE_NAMESPACE
#define _INTERFACE_NAMESPACE MVS
#endif

// uncomment to enable custom OpenCV data types
// (should be uncommented if OpenCV is not available)
#if !defined(_USE_OPENCV) && !defined(_USE_CUSTOM_CV)
#define _USE_CUSTOM_CV
#endif

#ifndef NO_ID
#define NO_ID std::numeric_limits<uint32_t>::max()
#endif


// S T R U C T S ///////////////////////////////////////////////////

#ifdef _USE_CUSTOM_CV

namespace cv {

// simple cv::Matx
template<typename Type, int m, int n>
class Matx
{
public:
	typedef Type value_type;
	inline Matx() {}
	#ifdef _USE_EIGEN
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF_VECTORIZABLE_FIXED_SIZE(Type,m*n)
	typedef Eigen::Matrix<Type,m,n,(n>1?Eigen::RowMajor : Eigen::ColMajor)> EMat;
	typedef Eigen::Map<const EMat> CEMatMap;
	typedef Eigen::Map<EMat> EMatMap;
	template<typename Derived>
	inline Matx(const Eigen::DenseBase<Derived>& rhs) { operator EMatMap () = rhs; }
	template<typename Derived>
	inline Matx& operator = (const Eigen::DenseBase<Derived>& rhs) { operator EMatMap () = rhs; return *this; }
	inline operator CEMatMap() const { return CEMatMap((const Type*)val); }
	inline operator EMatMap () { return EMatMap((Type*)val); }
	#endif
	static Matx eye() {
		Matx M;
		memset(M.val, 0, sizeof(Type)*m*n);
		const int shortdim(m < n ? m : n);
		for (int i = 0; i < shortdim; i++)
			M(i,i) = 1;
		return M;
	}
	Type operator()(int r, int c) const { return val[r*n+c]; }
	Type& operator()(int r, int c) { return val[r*n+c]; }
public:
	Type val[m*n];
};

// simple cv::Matx
template<typename Type>
class Point3_
{
public:
	typedef Type value_type;
	inline Point3_() {}
	#ifdef _USE_EIGEN
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF_VECTORIZABLE_FIXED_SIZE(Type,3)
	typedef Eigen::Matrix<Type,3,1> EVec;
	typedef Eigen::Map<EVec> EVecMap;
	template<typename Derived>
	inline Point3_(const Eigen::DenseBase<Derived>& rhs) { operator EVecMap () = rhs; }
	template<typename Derived>
	inline Point3_& operator = (const Eigen::DenseBase<Derived>& rhs) { operator EVecMap () = rhs; return *this; }
	inline operator const EVecMap () const { return EVecMap((Type*)this); }
	inline operator EVecMap () { return EVecMap((Type*)this); }
	#endif
public:
	Type x, y, z;
};

} // namespace cv
#endif
/*----------------------------------------------------------------*/


namespace _INTERFACE_NAMESPACE {

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
		std::transform(ext.begin(), ext.end(), ext.begin(), ::towlower);
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
bool Save<TYPE>(ArchiveSave& a, const TYPE& v) { \
	a.stream.write((const char*)&v, sizeof(TYPE)); \
	return true; \
} \
template<> \
bool Load<TYPE>(ArchiveLoad& a, TYPE& v) { \
	a.stream.read((char*)&v, sizeof(TYPE)); \
	return true; \
}

// Serialization support for basic types
ARCHIVE_DEFINE_TYPE(uint32_t)
ARCHIVE_DEFINE_TYPE(uint64_t)
ARCHIVE_DEFINE_TYPE(float)
ARCHIVE_DEFINE_TYPE(double)
#ifdef __APPLE__
  #ifdef __clang__
    ARCHIVE_DEFINE_TYPE(unsigned long)
  #endif
#endif

// Serialization support for cv::Matx
template<typename _Tp, int m, int n>
bool Save(ArchiveSave& a, const cv::Matx<_Tp,m,n>& _m) {
	a.stream.write((const char*)_m.val, sizeof(_Tp)*m*n);
	return true;
}
template<typename _Tp, int m, int n>
bool Load(ArchiveLoad& a, cv::Matx<_Tp,m,n>& _m) {
	a.stream.read((char*)_m.val, sizeof(_Tp)*m*n);
	return true;
}

// Serialization support for cv::Point3_
template<typename _Tp>
bool Save(ArchiveSave& a, const cv::Point3_<_Tp>& pt) {
	a.stream.write((const char*)&pt.x, sizeof(_Tp)*3);
	return true;
}
template<typename _Tp>
bool Load(ArchiveLoad& a, cv::Point3_<_Tp>& pt) {
	a.stream.read((char*)&pt.x, sizeof(_Tp)*3);
	return true;
}

// Serialization support for std::string
template<>
bool Save<std::string>(ArchiveSave& a, const std::string& s) {
	const size_t size(s.size());
	Save(a, size);
	if (size > 0)
		a.stream.write(&s[0], sizeof(char)*size);
	return true;
}
template<>
bool Load<std::string>(ArchiveLoad& a, std::string& s) {
	size_t size;
	Load(a, size);
	if (size > 0) {
		s.resize(size);
		a.stream.read(&s[0], sizeof(char)*size);
	}
	return true;
}

// Serialization support for std::vector
template<typename _Tp>
bool Save(ArchiveSave& a, const std::vector<_Tp>& v) {
	const size_t size(v.size());
	Save(a, size);
	for (size_t i=0; i<size; ++i)
		Save(a, v[i]);
	return true;
}
template<typename _Tp>
bool Load(ArchiveLoad& a, std::vector<_Tp>& v) {
	size_t size;
	Load(a, size);
	if (size > 0) {
		v.resize(size);
		for (size_t i=0; i<size; ++i)
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
			uint32_t width, height; // image resolution in pixels for all images sharing this camera (optional)
			Mat33d K; // camera's intrinsics matrix (normalized if image resolution not specified)
			Mat33d R; // camera's rotation matrix relative to the platform
			Pos3d C; // camera's translation vector relative to the platform

			Camera() : width(0), height(0) {}
			bool HasResolution() const { return width > 0 && height > 0; }
			bool IsNormalized() const { return !HasResolution(); }

			template <class Archive>
			void serialize(Archive& ar, const unsigned int version) {
				ar & name;
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
			Mat33d R; // platform's rotation matrix
			Pos3d C; // platform's translation vector in the global coordinate system

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
		std::string name; // image file name
		uint32_t platformID; // ID of the associated platform
		uint32_t cameraID; // ID of the associated camera on the associated platform
		uint32_t poseID; // ID of the pose of the associated platform

		template <class Archive>
		void serialize(Archive& ar, const unsigned int /*version*/) {
			ar & name;
			ar & platformID;
			ar & cameraID;
			ar & poseID;
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

	PlatformArr platforms; // array of platforms
	ImageArr images; // array of images
	VertexArr vertices; // array of reconstructed 3D points
	NormalArr verticesNormal; // array of reconstructed 3D points' normal (optional)
	ColorArr verticesColor; // array of reconstructed 3D points' color (optional)
	LineArr lines; // array of reconstructed 3D lines
	NormalArr linesNormal; // array of reconstructed 3D lines' normal (optional)
	ColorArr linesColor; // array of reconstructed 3D lines' color (optional)
	Mat44d transform; // transformation used to convert from absolute to relative coordinate system (optional)

	Interface() : transform(Mat44d::eye()) {}

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
			}
		}
	}
};
/*----------------------------------------------------------------*/

} // namespace _INTERFACE_NAMESPACE

#endif // _INTERFACE_MVS_H_
