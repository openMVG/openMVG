
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"
#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/sfm/sfm_data_triangulation.hpp"
#include "openMVG/geometry/rigid_transformation3D_srt.hpp"
#include "openMVG/geometry/Similarity3.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"
#include "third_party/htmlDoc/htmlDoc.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/graph/connectedComponent.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/multiview/essential.hpp"
#include "third_party/progress/progress.hpp"

#include "opencv2/line_descriptor.hpp"
#include "opencv2/core/utility.hpp"
#include <opencv2/imgproc.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/calib3d.hpp>

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <io.h>

#define MATCHES_DIST_THRESHOLD 25

using namespace openMVG;
using namespace openMVG::sfm;
using namespace cv;
using namespace cv::line_descriptor;
using namespace std;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::features;

IplImage *src;
IplImage *img;
bool drawing = false;
CvRect rect;
CvPoint origin;
bool istuNumOne = true;
string image_path = "C:\\Users\\Ethan\\Desktop\\8m\\";
int firstImgId = 0;
int secondImgId = 1;
void onMouse(int event, int x, int y, int flags, void *param)
{
	if (drawing)
	{
		rect.x = MIN(origin.x, x);
		rect.y = MIN(origin.y, y);
		rect.width = abs(origin.x - x);
		rect.height = abs(origin.y - y);
	}

	if (event == CV_EVENT_LBUTTONDOWN && !CV_EVENT_MOUSEMOVE)
	{
		drawing = true;
		origin = cvPoint(x, y);
		rect = cvRect(x, y, 0, 0);
	}
	else if (event == CV_EVENT_LBUTTONUP)
	{
		drawing = false;
		if (rect.height == 0 || rect.width == 0)
		{
			//cvDestroyWindow("ScreenShot");
			return;
		}
		img = cvCreateImage(cvSize(rect.width, rect.height), src->depth, src->nChannels);
		cout << rect.x << " " << rect.y << " " << rect.height << " " << rect.width << endl;
		cvSetImageROI(src, rect);
		cvCopy(src, img, 0);
		cvResetImageROI(src);

		//cvNamedWindow("ScreenShot", 1);
		//cvShowImage("ScreenShot", img);
		if (istuNumOne)
		{
			string filename1 = image_path + "picSmall01.jpg";
			const char * file1 = filename1.c_str();
			cvSaveImage(file1, img);
			istuNumOne = false;
		}
		else
		{
			string filename2 = image_path + "picSmall02.jpg";
			const char *file2 = filename2.c_str();
			cvSaveImage(file2, img);
			istuNumOne = true;
		}

		//cvWaitKey(0);


		//cvDestroyWindow("ScreenShot");
		return;
	}
}

bool sortdes(const KeyLine &k1, const KeyLine &k2)
{
	return k1.lineLength > k2.lineLength;
}

void RotationMatrixToEulerAnglesXYZ(Eigen::Matrix<double, 3, 3>  R,double* euler) {
	//Z
	euler[0] = atan2(-R(0, 1), R(0, 0));
	//Y
	euler[1] = atan2(R(0, 2), sqrt(R(0, 0)*R(0, 0) + R(0, 1)*R(0, 1)));
	//double Y2 = asin(R(0, 2));
	//X
	euler[2] = atan2(-R(1, 2), R(2, 2));
	euler[0] = euler[0] * 180 / M_PI;
	cout << " angleZ: " << euler[0];
	euler[1] = euler[1] * 180 / M_PI;
	cout << " angleY: " << euler[1];
	//Y2 = Y2 * 180 / M_PI;
	//cout << " angleY2: " << Y2;
	euler[2] = euler[2] * 180 / M_PI;
	cout << " angleX: " << euler[2] << endl;
	// getRotMatrixXYZ(euler[0], euler[1], euler[2]);
}

void RotationMatrixToEulerAnglesZXY(Eigen::Matrix<double, 3, 3> rotMatrix, double *euler) {
	//cout <<"rotMatrix(1, 0)"<<rotMatrix(1, 0) << endl;
	//cout << "rotMatrix(0, 0)" << rotMatrix(0, 0) << endl;
	euler[0] = -atan2(rotMatrix(0, 1), rotMatrix(1, 1));
	euler[0] = (euler[0] * 360) / (2 * M_PI);
	cout << "angleZ:" << euler[0] << "  ";
	//double pitch = atan2(-rotMatrix(2, 0) , sqrt(rotMatrix(2, 1) * rotMatrix(2, 1) + rotMatrix(2, 2) * rotMatrix(2, 2)));
	euler[1] = atan2(-rotMatrix(2, 0), rotMatrix(2, 2));
	euler[1] = (euler[1] * 360) / (2 * M_PI);
	cout << "angleY:" << euler[1] << "  ";
	euler[2] = -asin(-rotMatrix(2, 1));
	euler[2] = (euler[2] * 360) / (2 * M_PI);
	cout << "angleX:" << euler[2] << "  ";
	cout << endl;
}

void RotationMatrixToEulerAnglesZYX(Eigen::Matrix<double, 3, 3>  R, double* euler) {
	//Z
	euler[0] = atan2(R(1, 0), R(0, 0));
	//Y
	euler[1] = atan2(-R(2, 0), sqrt(R(2, 1)*R(2, 1) + R(2, 2)*R(2, 2)));
	//X
	euler[2] = atan2(R(2, 1), R(2, 2));
	euler[0] = euler[0] * 180 / M_PI;
	cout << " angleZ: " << euler[0];
	euler[1] = euler[1] * 180 / M_PI;
	cout << " angleY: " << euler[1];
	euler[2] = euler[2] * 180 / M_PI;
	cout << " angleX: " << euler[2] << endl;

}

Eigen::Matrix<double, 3, 3> getRotMatrixZXY(double angleZ, double angleY, double angleX) {
	//转化为弧度
	angleZ = (angleZ * 2 * M_PI) / 360;
	angleX = (angleX * 2 * M_PI) / 360;
	angleY = (angleY * 2 * M_PI) / 360;

	Eigen::Matrix<double, 3, 3> rotX, rotY, rotZ, rotationMatrix;
	//绕Z轴
	rotZ << cos(angleZ), -sin(angleZ), 0, sin(angleZ), cos(angleZ), 0, 0, 0, 1;
	//绕X轴
	rotX << 1, 0, 0, 0, cos(angleX), -sin(angleX), 0, sin(angleX), cos(angleX);
	//绕Y轴
	rotY << cos(angleY), 0, sin(angleY), 0, 1, 0, -sin(angleY), 0, cos(angleY);

	//cout << "rotTaxis:\n" << rotYaxis << endl;
	rotationMatrix = rotZ*rotX*rotY;

	cout << "rotationMatrix:\n" << rotationMatrix << endl;

	return rotationMatrix;
}

Eigen::Matrix<double, 3, 3> getRotMatrixZYX(double angleZ, double angleY, double angleX) {
	//转化为弧度
	angleZ = (angleZ * 2 * M_PI) / 360;
	angleX = (angleX * 2 * M_PI) / 360;
	angleY = (angleY * 2 * M_PI) / 360;

	Eigen::Matrix<double, 3, 3> rotX, rotY, rotZ, rotationMatrix;
	//绕Z轴
	rotZ << cos(angleZ), -sin(angleZ), 0, sin(angleZ), cos(angleZ), 0, 0, 0, 1;
	//绕X轴
	rotX << 1, 0, 0, 0, cos(angleX), -sin(angleX), 0, sin(angleX), cos(angleX);
	//绕Y轴
	rotY << cos(angleY), 0, sin(angleY), 0, 1, 0, -sin(angleY), 0, cos(angleY);

	//cout << "rotTaxis:\n" << rotYaxis << endl;
	rotationMatrix = rotZ*rotY*rotX;

	//cout << "rotationMatrix:\n" << rotationMatrix << endl;

	return rotationMatrix;
}

Eigen::Matrix<double, 3, 3> getRotMatrixXYZ(double angleZ, double angleY, double angleX) {
	//转化为弧度
	angleZ = (angleZ * 2 * M_PI) / 360;
	angleX = (angleX * 2 * M_PI) / 360;
	angleY = (angleY * 2 * M_PI) / 360;

	Eigen::Matrix<double, 3, 3> rotX, rotY, rotZ, rotationMatrix;
	//绕Z轴
	rotZ << cos(angleZ), -sin(angleZ), 0, sin(angleZ), cos(angleZ), 0, 0, 0, 1;
	//绕X轴
	rotX << 1, 0, 0, 0, cos(angleX), -sin(angleX), 0, sin(angleX), cos(angleX);
	//绕Y轴
	rotY << cos(angleY), 0, sin(angleY), 0, 1, 0, -sin(angleY), 0, cos(angleY);

	//cout << "rotTaxis:\n" << rotYaxis << endl;
	rotationMatrix = rotX*rotY*rotZ;

	cout << "rotationMatrix:\n" << rotationMatrix << endl;

	return rotationMatrix;
}

Eigen::Matrix<double, 3, 3> getMatrixFromTxt(string path, string identity)
{
	Eigen::Matrix<double, 3, 3> matrix;
	std::ifstream fin(path, std::ios::in);
	string line;
	while (fin >> line)
	{
		if (line == identity)
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					fin >> line;
					matrix(i, j) = atof(line.c_str());
				}
			}
	}

	return matrix;
}

 double* getAngleFromTxt(string path, string identity)
{
	double angles[3];
	std::ifstream fin(path, std::ios::in);
	string line;
	while (fin >> line)
	{
		if (line == identity)
			for (int i = 0; i < 3; i++)
			{
				fin >> line;
				angles[i] = atof(line.c_str());
			}
	}

	return angles;
}
// 获得 RotationMatrix
 Eigen::Matrix<double, 3, 3> getRotationMatrix(string path)
{
	return getMatrixFromTxt(path, "RotationMatrixFromVector:");
}
// 获得 RemappedRotationMatrix
 Eigen::Matrix<double, 3, 3> getRemappedRotationMatrix(string path)
{
	return getMatrixFromTxt(path, "RemappedRotationMatrixFromVector:");
}

// 获得一个路径下所有的 txt 文件
vector<string> listTxtFiles(string path)
{
	vector<string> txt_files;

	intptr_t  hFile = 0;
	struct _finddata_t fileInfo;
	string pathName, exdName;
	const char* p = pathName.assign(path).append("\\*").c_str();
	if ((hFile = _findfirst(p, &fileInfo)) == -1) {
		return txt_files;
	}
	do
	{
		const string fileName(fileInfo.name);
		if (fileName.length() > 3)
		{
			const string extension = fileName.substr(fileName.length() - 3, 3);
			if (extension == "txt")
			{
				txt_files.push_back(path + fileName);
			}
		}
		std::cout << errno<<endl;
		//	cout << fileInfo.name << (fileInfo.attrib&_A_SUBDIR ? "[folder]" : "[file]") << endl;
	} while (_findnext(hFile, &fileInfo) == 0);

	_findclose(hFile);

	return txt_files;
}

vector<string> listJpgFiles(string path)
{
	vector<string> txt_files;

	intptr_t  hFile = 0;
	struct _finddata_t fileInfo;
	string pathName, exdName;

	if ((hFile = _findfirst(pathName.assign(path).append("\\*").c_str(), &fileInfo)) == -1) {
		return txt_files;
	}
	do
	{
		const string fileName(fileInfo.name);
		if (fileName.length() > 3)
		{
			const string extension = fileName.substr(fileName.length() - 3, 3);
			if (extension == "jpg")
			{
				txt_files.push_back(path + fileName);
			}
		}

		//	cout << fileInfo.name << (fileInfo.attrib&_A_SUBDIR ? "[folder]" : "[file]") << endl;
	} while (_findnext(hFile, &fileInfo) == 0);

	_findclose(hFile);

	return txt_files;
}

 double getDistance(double tx, double ty, double tz, double rx, double ry, double rz) {
	return (sqrt((tx - rx) * (tx - rx) + (ty - ry) * (ty - ry) + (tz - rz) * (tz - rz)));

}

 void getFuYang(double tx, double ty, double tz, double rx, double ry, double rz) {
	double dis22_ = getDistance(tx, ty, tz, tx, ty, rz);
	double dis24 = getDistance(tx, ty, tz, rx, ry, rz);
	double theta = asin(dis22_ / dis24);
	double fuyang = (theta * 360) / (2 * M_PI);
	cout << "俯仰角:" << fuyang << endl;

}
// void getShuiPing(double tx, double ty, double tz, double rx, double ry, double rz){
//	double dis1_3 = getDistance(tx, ty, rz, rx, ry, rz);
//	double dis1_1__ = getDistance(tx, ty, rz, rx, ty, rz);
//	double alpha = asin(dis1_1__ / dis1_3);
//	double shuiping = (alpha * 360) / (2 * M_PI);
//	cout << "水平角:" << shuiping << endl;
//
//}

 void getShuiPing(double tx, double ty, double tz, double rx, double ry, double rz) {
	double angle;
	double vecX, vecY;
	if (tz >= rz) {
		//tz为终点的纵坐标，rz为始点的纵坐标
		//获取到矢量
		vecX = tx - rx;
		vecY = ty - ry;
		if (rx <= tx) {
			//t在r的右边，终点在始点的右边
			angle = acos((0 * vecX + 1 * vecY) / (1 * sqrt(vecX * vecX + vecY * vecY)));
		}
		else {
			//t在r的左边，终点在始点的左边
			angle = 2 * M_PI - acos((0 * vecX + 1 * vecY) / (1 * sqrt(vecX * vecX + vecY * vecY)));
		}
	}
	else {
		//rz为终点的纵坐标，tz为始点的纵坐标
		//获取到矢量
		vecX = rx - tx;
		vecY = ry - ty;
		if (tx <= rx) {
			//r在t的右边，终点在始点的右边
			angle = acos((0 * vecX + 1 * vecY) / (1 * sqrt(vecX * vecX + vecY * vecY)));
		}
		else {
			//r在t的左边，终点在始点的左边
			angle = 2 * M_PI - acos((0 * vecX + 1 * vecY) / (1 * sqrt(vecX * vecX + vecY * vecY)));
		}

	}
	angle = (angle * 360) / (2 * M_PI);
	cout << "水平角:" << angle << endl;

}
// 获得一个路径下所有的jpg 文件

cv::Mat getFundament(string image_path)
{
	////////////*求出两个图像之间的F的关系，求出极线，求出同名像点*/
	//利用P 和 c 求两图之间的F
	/////testPicture01
	//Mat promatric1 = (Mat_<double>(3, 4) << 3359.61, 856.293, 2319.35, 464.852,
	//	-846.7, 3568.64, 1213.34, 136.077,
	//	-0.0675598, 0.0444401, 0.996725, -0.100325);
	//Mat promatric2 = (Mat_<double>(3, 4) << 3012.23, 2118.43, 1959.17, 1117.13,
	//	-1876.79, 3250.25, 915.311, -28.6339,
	//	-0.0364896, 0.138705, 0.989661, -0.128688);
	////第一个相机中心
	//Mat C = (Mat_<double>(4, 1) << -0.174699, -0.111467, 0.0937833, 1.0);


	////testPicture02
	//Mat promatric1 = (Mat_<double>(3, 4) << 3555.77, -14.1684,   2067.54, -182.9,
	//	53.9075,   3575.71,   1539.02, -299.02,
	//	0.0031093, 0.0168167,  0.999854, -0.134442);
	//Mat promatric2 = (Mat_<double>(3, 4) << 3691.34,   350.491,  1714.64,   397.699,
	//	-89.0851,   3656.74,   1308.67, -234.847,
	//	0.0790208, 0.0783662,  0.993788, -0.125345);
	////第一个相机中心
	//Mat C = (Mat_<double>(4, 1) << -0.026432,0.0263055,	0.134101, 1.0);


	///testPicture03
	fstream inProC(image_path + "ProjMat&C.txt", ios::in);

	cv::Mat promatric1(3, 4, CV_64FC1);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
			inProC >> promatric1.at<double>(i, j);
	}

	cv::Mat promatric2(3, 4, CV_64FC1);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
			inProC >> promatric2.at<double>(i, j);
	}
	//第一个相机中心
	cv::Mat C(4, 1, CV_64FC1);
	for (int i = 0; i < 3; i++)
		inProC >> C.at<double>(i, 0);
	C.at<double>(3, 0) = 1.0;

	inProC.close();
	cout << setprecision(15) << promatric1 << endl << promatric2 << endl << C << endl;

	cv::Mat ee = promatric2*C;
	//cout << ee << endl << endl;
	cv::Mat eInvSym = cv::Mat::zeros(3, 3, DataType<double>::type);
	eInvSym.at<double>(0, 1) = -ee.at<double>(2);
	eInvSym.at<double>(0, 2) = ee.at<double>(1);
	eInvSym.at<double>(1, 0) = ee.at<double>(2);
	eInvSym.at<double>(1, 2) = -ee.at<double>(0);
	eInvSym.at<double>(2, 0) = -ee.at<double>(1);
	eInvSym.at<double>(2, 1) = ee.at<double>(0);
	//cout << eInvSym << endl << endl;

	cv::Mat pro1inv = promatric1.inv(DECOMP_SVD);   //求pro的伪逆矩阵
	cv::Mat FundamentEPP = eInvSym*promatric2*pro1inv;
	cout << FundamentEPP << endl << endl;

	/*
	//利用点取点对，求基本矩阵
	const int PointNum=30;
	vector<Point2f> pointTu1(PointNum);
	vector<Point2f> pointTu2(PointNum);

	ifstream infile;
	infile.open("E:\\天线测量\\照片\\tianxianSmall\\tuQ1.txt");

	for (int i = 0; i < PointNum; i++)
	{
	infile >> pointTu1[i].x >> pointTu1[i].y;
	}
	infile.close();

	infile.open("E:\\天线测量\\照片\\tianxianSmall\\tuQ2.txt");

	for (int i = 0; i < PointNum; i++)
	{
	infile >> pointTu2[i].x >> pointTu2[i].y;
	}
	infile.close();

	Mat FundamentalMat= findFundamentalMat(pointTu1, pointTu2, FM_RANSAC, 3, 0.99);
	cout << FundamentalMat << endl;
	*/


	return FundamentEPP;
}



int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "Global Structure from Motion:\n"
    << "-----------------------------------------------------------\n"
    << "Open Source implementation of the paper:\n"
    << "\"Global Fusion of Relative Motions for "
    << "Robust, Accurate and Scalable Structure from Motion.\"\n"
    << "Pierre Moulon, Pascal Monasse and Renaud Marlet. "
    << " ICCV 2013." << std::endl
    << "------------------------------------------------------------"
    << std::endl;


  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sMatchesDir;
  std::string sOutDir = "";
  int iRotationAveragingMethod = int (ROTATION_AVERAGING_L2);
  int iTranslationAveragingMethod = int (TRANSLATION_AVERAGING_SOFTL1);
  std::string sIntrinsic_refinement_options = "ADJUST_ALL";

  cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('r', iRotationAveragingMethod, "rotationAveraging") );
  cmd.add( make_option('t', iTranslationAveragingMethod, "translationAveraging") );
  cmd.add( make_option('f', sIntrinsic_refinement_options, "refineIntrinsics") );

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file] path to a SfM_Data scene\n"
    << "[-m|--matchdir] path to the matches that corresponds to the provided SfM_Data scene\n"
    << "[-o|--outdir] path where the output data will be stored\n"
    << "\n[Optional]\n"
    << "[-r|--rotationAveraging]\n"
      << "\t 1 -> L1 minimization\n"
      << "\t 2 -> L2 minimization (default)\n"
    << "[-t|--translationAveraging]:\n"
      << "\t 1 -> L1 minimization\n"
      << "\t 2 -> L2 minimization of sum of squared Chordal distances\n"
      << "\t 3 -> SoftL1 minimization (default)\n"
    << "[-f|--refineIntrinsics] Intrinsic parameters refinement option\n"
      << "\t ADJUST_ALL -> refine all existing parameters (default) \n"
      << "\t NONE -> intrinsic parameters are held as constant\n"
      << "\t ADJUST_FOCAL_LENGTH -> refine only the focal length\n"
      << "\t ADJUST_PRINCIPAL_POINT -> refine only the principal point position\n"
      << "\t ADJUST_DISTORTION -> refine only the distortion coefficient(s) (if any)\n"
      << "\t -> NOTE: options can be combined thanks to '|'\n"
      << "\t ADJUST_FOCAL_LENGTH|ADJUST_PRINCIPAL_POINT\n"
      <<      "\t\t-> refine the focal length & the principal point position\n"
      << "\t ADJUST_FOCAL_LENGTH|ADJUST_DISTORTION\n"
      <<      "\t\t-> refine the focal length & the distortion coefficient(s) (if any)\n"
      << "\t   \n"
      <<      "\t\t-> refine the principal point position & the distortion coefficient(s) (if any)\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (iRotationAveragingMethod < ROTATION_AVERAGING_L1 ||
      iRotationAveragingMethod > ROTATION_AVERAGING_L2 )  {
    std::cerr << "\n Rotation averaging method is invalid" << std::endl;
    return EXIT_FAILURE;
  }

  const cameras::Intrinsic_Parameter_Type intrinsic_refinement_options =
    cameras::StringTo_Intrinsic_Parameter_Type(sIntrinsic_refinement_options);

  if (iTranslationAveragingMethod < TRANSLATION_AVERAGING_L1 ||
      iTranslationAveragingMethod > TRANSLATION_AVERAGING_SOFTL1 )  {
    std::cerr << "\n Translation averaging method is invalid" << std::endl;
    return EXIT_FAILURE;
  }

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  const std::string sImage_describer = stlplus::create_filespec(sMatchesDir, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    std::cerr << "Invalid: "
      << sImage_describer << " regions type file." << std::endl;
    return EXIT_FAILURE;
  }

  // Features reading
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
  if (!feats_provider->load(sfm_data, sMatchesDir, regions_type)) {
    std::cerr << std::endl
      << "Invalid features." << std::endl;
    return EXIT_FAILURE;
  }
  // Matches reading
  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
  if // Try to read the two matches file formats
  (
    !(matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.e.txt")) ||
      matches_provider->load(sfm_data, stlplus::create_filespec(sMatchesDir, "matches.e.bin")))
  )
  {
    std::cerr << std::endl
      << "Invalid matches file." << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
  {
    if (!stlplus::folder_create(sOutDir))
    {
      std::cerr << "\nCannot create the output directory" << std::endl;
    }
  }

  //---------------------------------------
  // Global SfM reconstruction process
  //---------------------------------------

  openMVG::system::Timer timer;
  GlobalSfMReconstructionEngine_RelativeMotions sfmEngine(
    sfm_data,
    sOutDir,
    stlplus::create_filespec(sOutDir, "Reconstruction_Report.html"));

  // Configure the features_provider & the matches_provider
  sfmEngine.SetFeaturesProvider(feats_provider.get());
  sfmEngine.SetMatchesProvider(matches_provider.get());

  // Configure reconstruction parameters
  sfmEngine.Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);

  // Configure motion averaging method
  sfmEngine.SetRotationAveragingMethod(
    ERotationAveragingMethod(iRotationAveragingMethod));
  sfmEngine.SetTranslationAveragingMethod(
    ETranslationAveragingMethod(iTranslationAveragingMethod));

  if (sfmEngine.Process())
  {
	  std::cout << std::endl << " Total Ac-Global-Sfm took (s): " << timer.elapsed() << std::endl;

	  std::cout << "...Generating SfM_Report.html" << std::endl;
	  Generate_SfM_Report(sfmEngine.Get_SfM_Data(),
		  stlplus::create_filespec(sOutDir, "SfMReconstruction_Report.html"));

	  //获取SFM DATA
	  SfM_Data my_sfm_data;
	  my_sfm_data = sfmEngine.Get_SfM_Data();


	  //旋转转化矩阵的计算
	  //获取虚拟空间的旋转矩阵
	  int cameraNum = my_sfm_data.poses.size();
	  vector<Eigen::Matrix<double, 3, 3>> rotations;
	  for (int i = 0; i < my_sfm_data.poses.size(); i++) {
		  rotations.push_back(my_sfm_data.poses[i].rotation());
		  cout << rotations[i] << endl;
	  }
	  //获取现实空间中的旋转矩阵
	  vector<string> txtNames;

	  //三角测量的两个图像序号

	  txtNames = listTxtFiles(image_path);
	  vector<Eigen::Matrix<double, 3, 3>> rotationsAndroid;
	  Eigen::Matrix<double, 3, 3> tempMat1, tempMat2;
	  double* tempAngle;
	  for (int i = 0; i < cameraNum; i++) {
		  tempAngle = getAngleFromTxt(txtNames[i], "Orientation_NEW_API:");
		  //cout << tempAngle[0] << " " << tempAngle[1] << " "<<tempAngle[2] << " "<<endl;
		  double angleZ = -tempAngle[0];
		  double angleY = tempAngle[2];
		  double angleX = -tempAngle[1];
		  tempMat1 = getRotMatrixZXY(angleZ, angleY, angleX);
		  tempMat2 = getRotationMatrix(txtNames[i]);
		  //cout << tempMat2 << endl;
		  rotationsAndroid.push_back(tempMat1);
	  }


	  Eigen::Matrix<double, 3, 3> rotX, rotZ;
	  double angleX, angleZ;
	  angleX = (180 * 2 * M_PI) / 360;
	  angleZ = (-90 * 2 * M_PI) / 360;
	  rotX << 1, 0, 0, 0, cos(angleX), -sin(angleX), 0, sin(angleX), cos(angleX);
	  rotZ << cos(angleZ), -sin(angleZ), 0, sin(angleZ), cos(angleZ), 0, 0, 0, 1;
	  Eigen::Matrix<double, 3, 3> tempRot;
	  for (int i = 0; i < cameraNum; i++) {
		  //将R取反得到Rc
		  tempRot = rotations[i].inverse();
		  rotations[i] = tempRot;
		  //获取世界空间中的相机旋转矩阵
		  rotationsAndroid[i] = rotationsAndroid[i] * rotX* rotZ;
		  cout << "Remap:" << rotationsAndroid[i] << endl << endl;

	  }
	  cout << "从虚拟空间旋转至现实空间的旋转矩阵:" << endl;
	  double eulerT[3];
	  double eulerTotal[3] = { 0.0, 0.0, 0.0 };
	  std::vector<double> eulerVector;
	  std::vector<Eigen::Matrix<double, 3, 3>> rt(cameraNum);

	  cout << endl << "XYZ" << endl;
	  for (unsigned int i = 0; i < cameraNum; i++) {
		  rt[i] = rotationsAndroid[i] * (rotations[i].inverse());
		  //cout << rt[i] << endl;
		  RotationMatrixToEulerAnglesXYZ(rt[i], eulerT);

	  }
	  cout << endl << "ZXY" << endl;
	  for (unsigned int i = 0; i < cameraNum; i++) {

		  RotationMatrixToEulerAnglesZXY(rt[i], eulerT);

	  }
	  cout << endl << "ZYX" << endl;
	  for (unsigned int i = 0; i < cameraNum; i++) {

		  RotationMatrixToEulerAnglesZYX(rt[i], eulerT);
		  eulerVector.push_back(eulerT[0]);
		  eulerVector.push_back(eulerT[1]);
		  eulerVector.push_back(eulerT[2]);
	  }
	  //判定是否存在180左右的数值,当一个数值存在大于175，需要检测
	  if (abs(eulerVector[0]) >175) {
		  int negativeNum = 0;
		  int positiveNum = 0;
		  for (unsigned int i = 0; i < eulerVector.size() / 3; i++) {
			  if (eulerVector[i * 3] > 0) {
				  positiveNum++;
			  }
			  else {
				  negativeNum++;
			  }
		  }
		  if ((positiveNum != eulerVector.size() / 3) || (negativeNum != eulerVector.size() / 3)) {
			  //有数据需要取反
			  /*cout << endl <<"positiveNum"<< positiveNum << endl;
			  cout << endl << "negativeNum" << positiveNum << endl;*/
			  if (positiveNum >= negativeNum) {
				  for (unsigned int i = 0; i < eulerVector.size() / 3; i++) {
					  eulerVector[i * 3] = abs(eulerVector[i * 3]);
				  }
			  }
			  else {
				  for (unsigned int i = 0; i < eulerVector.size() / 3; i++) {
					  eulerVector[i * 3] = -abs(eulerVector[i * 3]);
				  }
			  }
		  }

	  }

	  if (abs(eulerVector[1]) >175) {
		  int negativeNum = 0;
		  int positiveNum = 0;
		  for (unsigned int i = 0; i < eulerVector.size() / 3; i++) {
			  if (eulerVector[i * 3 + 1] > 0) {
				  positiveNum++;
			  }
			  else {
				  negativeNum++;
			  }
		  }
		  if ((positiveNum != eulerVector.size() / 3) || (negativeNum != eulerVector.size() / 3)) {
			  //有数据需要取反
			  /*cout << endl <<"positiveNum"<< positiveNum << endl;
			  cout << endl << "negativeNum" << positiveNum << endl;*/
			  if (positiveNum >= negativeNum) {
				  for (unsigned int i = 0; i < eulerVector.size() / 3; i++) {
					  eulerVector[i * 3 + 1] = abs(eulerVector[i * 3 + 1]);
				  }
			  }
			  else {
				  for (unsigned int i = 0; i < eulerVector.size() / 3; i++) {
					  eulerVector[i * 3 + 1] = -abs(eulerVector[i * 3 + 1]);
				  }
			  }
		  }

	  }

	  if (abs(eulerVector[2]) >175) {
		  int negativeNum = 0;
		  int positiveNum = 0;
		  for (unsigned int i = 0; i < eulerVector.size() / 3; i++) {
			  if (eulerVector[i * 3 + 2] > 0) {
				  positiveNum++;
			  }
			  else {
				  negativeNum++;
			  }
		  }
		  if ((positiveNum != eulerVector.size() / 3) || (negativeNum != eulerVector.size() / 3)) {
			  //有数据需要取反
			  /*cout << endl <<"positiveNum"<< positiveNum << endl;
			  cout << endl << "negativeNum" << positiveNum << endl;*/
			  if (positiveNum >= negativeNum) {
				  for (unsigned int i = 0; i < eulerVector.size() / 3; i++) {
					  eulerVector[i * 3 + 2] = abs(eulerVector[i * 3 + 2]);
				  }
			  }
			  else {
				  for (unsigned int i = 0; i < eulerVector.size() / 3; i++) {
					  eulerVector[i * 3 + 2] = -abs(eulerVector[i * 3 + 2]);
				  }
			  }
		  }

	  }

	  for (unsigned int i = 0; i < eulerVector.size() / 3; i++) {
		  eulerTotal[0] += eulerVector[i * 3];
		  eulerTotal[1] += eulerVector[i * 3 + 1];
		  eulerTotal[2] += eulerVector[i * 3 + 2];
	  }

	  eulerTotal[0] = eulerTotal[0] / (eulerVector.size() / 3);
	  eulerTotal[1] = eulerTotal[1] / (eulerVector.size() / 3);
	  eulerTotal[2] = eulerTotal[2] / (eulerVector.size() / 3);
	  Eigen::Matrix<double, 3, 3> finalrt;
	  finalrt = getRotMatrixZYX(eulerTotal[0], eulerTotal[1], eulerTotal[2]);
	  //计算均方值
	  double RSME = 0.0;
	  for (unsigned int i = 0; i < eulerVector.size() / 3; i++) {
		  RSME += sqrt((eulerTotal[0] - eulerVector[i * 3 + 0])*(eulerTotal[0] - eulerVector[i * 3 + 0])) +
			  +sqrt((eulerTotal[1] - eulerVector[i * 3 + 1])*(eulerTotal[1] - eulerVector[i * 3 + 1]))
			  + sqrt((eulerTotal[2] - eulerVector[i * 3 + 2])*(eulerTotal[2] - eulerVector[i * 3 + 2]));
	  }
	  cout << "Angle Z: " << eulerTotal[0] << " ";
	  cout << "Angle Y: " << eulerTotal[1] << " ";
	  cout << "Angle X: " << eulerTotal[2] << endl;
	  cout << "RSME: " << RSME << endl;
	  cout << "Average RSME:" << RSME / eulerVector.size() << endl;


	  fstream tr0(image_path + "transformRot.txt", ios::out);
	  for (unsigned int i = 0; i < cameraNum; i++) {
		  tr0 << "rt" << i << endl << rt[i] << endl;
	  }
	  //tr0 << "rt4_2:" << rt4_2 << endl;
	  tr0 << "final_rt:" << finalrt << endl;
	  tr0.close();

	  //通过直线提取与极限约束获取顶点对以进行三角测量

	  Triangulation trianObj;

	  //获取view，包含了图像，内参，外参
	  View *view = my_sfm_data.views.at(firstImgId).get();
	  //获取内外参
	  IntrinsicBase * cam = my_sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
	  //Eigen::Matrix<double, 3, 3> intrinsic1 = my_sfm_data.GetIntrinsics().at(view->id_intrinsic).get(;
	  Pose3 pose = my_sfm_data.GetPoseOrDie(view);


	  //同理获取view，加入trianObj
	  View *view2 = my_sfm_data.views.at(secondImgId).get();
	  //获取内外参
	  IntrinsicBase * cam2 = my_sfm_data.GetIntrinsics().at(view2->id_intrinsic).get();
	  Pose3 pose2 = my_sfm_data.GetPoseOrDie(view2);

	  /*fstream outProj(image_path + "ProjMat&C.txt", ios::out);
	  outProj << setprecision(15) << (cam->get_projective_equivalent(pose)) << endl;
	  outProj << setprecision(15) << (cam2->get_projective_equivalent(pose2)) << endl;
	  outProj << setprecision(15) << pose.center();
	  outProj.close();

	  cv::Mat myFundamentEPP = getFundament(image_path);
	  cout << endl << "Fundamental Matrix:" << myFundamentEPP << endl;*/

	  vector<string> jpgNames;
	  jpgNames = listJpgFiles(image_path);
	  string img1 = jpgNames[firstImgId];
	  src = cvLoadImage(img1.c_str());
	  if (!src) {
		  printf("Could not load image file: %s\n");
		  exit(0);
	  }
	  cvNamedWindow("screenshot", 0);
	  //缩小提取的区域
	  cvSetMouseCallback("screenshot", onMouse, NULL);//捕捉鼠标
	  cvShowImage("screenshot", src);
	  cvWaitKey(0);
	  Point2d pointOne = CvPoint(rect.x, rect.y);
	  cvReleaseImage(&src);
	  cvReleaseImage(&img);



	  //****第二幅图像
	  string img2 = jpgNames[secondImgId];
	  src = cvLoadImage(img2.c_str());
	  if (!src) {
		  printf("Could not load image file: %s\n");
		  exit(0);
	  }
	  cvNamedWindow("screenshot", 0);
	  //缩小提取的区域
	  cvSetMouseCallback("screenshot", onMouse, NULL);//捕捉鼠标
	  cvShowImage("screenshot", src);
	  cvWaitKey(0);
	  Point2d pointTwo = CvPoint(rect.x, rect.y);
	  cvReleaseImage(&src);
	  cvReleaseImage(&img);
	  cvDestroyAllWindows();



	  cout << pointOne << endl;
	  cout << pointTwo << endl;

	  String image_path1 = image_path + "picSmall01.jpg"; //parser.get<String>( 0 );
	  String image_path2 = image_path + "picSmall02.jpg"; //parser.get<String>( 1 );
	  if (image_path1.empty() || image_path2.empty())
	  {
		  //help();
		  return -1;
	  }

	  /* load image */
	  cv::Mat imageMat1 = imread(image_path1, 1);
	  cv::Mat imageMat2 = imread(image_path2, 1);

	  if (imageMat1.data == NULL || imageMat2.data == NULL)
	  {
		  cout << "Error, images could not be loaded. Please, check their path" << endl;
	  }

	  /* create binary masks */
	  cv::Mat mask1 = cv::Mat::ones(imageMat1.size(), CV_8UC1);
	  cv::Mat mask2 = cv::Mat::ones(imageMat2.size(), CV_8UC1);

	  /* create a pointer to a BinaryDescriptor object with default parameters */
	  Ptr<BinaryDescriptor> bd = BinaryDescriptor::createBinaryDescriptor();

	  /* compute lines and descriptors */
	  vector<KeyLine> keylines1, keylines2;
	  cv::Mat descr1, descr2;

	  (*bd)(imageMat1, mask1, keylines1, descr1, false, false);
	  (*bd)(imageMat2, mask2, keylines2, descr2, false, false);

	  //cout << "**"<<keylines1.size() << endl;
	  //cout <<"**"<< keylines2.size() << endl;



	  ////////////*求出两个图像之间的F的关系，求出极线，求出同名像点*/
	  //cv::Mat FundamentEPP = getFundament(image_path);


	  Eigen::Matrix<double, 3, 4> proj1, proj2;
	  Vec3 c;
	  proj1 = cam->get_projective_equivalent(pose);
	  proj2 = cam2->get_projective_equivalent(pose2);
	  c = my_sfm_data.poses[firstImgId].center();
	  //cout << proj1<<endl;
	  //cout << proj2<<endl;
	  //cout << c << endl;
	  //投影矩阵
	  cv::Mat promatric1(3, 4, CV_64FC1);
	  cv::Mat promatric2(3, 4, CV_64FC1);
	  for (int i = 0; i < 3; i++)
	  {
		  for (int j = 0; j < 4; j++) {
			  promatric1.at<double>(i, j) = proj1(i, j);
			  promatric2.at<double>(i, j) = proj2(i, j);
		  }
	  }
	  //第一个相机中心
	  cv::Mat C(4, 1, CV_64FC1);
	  for (int i = 0; i < 3; i++) {
		  C.at<double>(i, 0) = c(i, 0);
	  }
	  C.at<double>(3, 0) = 1.0;
	  //cout << promatric1 << endl;
	  //cout << promatric2 << endl;
	  //cout << C << endl;

	  cv::Mat ee = promatric2*C;
	  //cout << ee << endl << endl;
	  cv::Mat eInvSym = cv::Mat::zeros(3, 3, DataType<double>::type);
	  eInvSym.at<double>(0, 1) = -ee.at<double>(2);
	  eInvSym.at<double>(0, 2) = ee.at<double>(1);
	  eInvSym.at<double>(1, 0) = ee.at<double>(2);
	  eInvSym.at<double>(1, 2) = -ee.at<double>(0);
	  eInvSym.at<double>(2, 0) = -ee.at<double>(1);
	  eInvSym.at<double>(2, 1) = ee.at<double>(0);
	  //cout << eInvSym << endl << endl;

	  cv::Mat pro1inv = promatric1.inv(DECOMP_SVD);   //求pro的伪逆矩阵
	  cv::Mat FundamentEPP = eInvSym*promatric2*pro1inv;
	  //cout << FundamentEPP << endl << endl;

	  Eigen::Matrix<double, 3, 3> Fundament = F_from_P(cam->get_projective_equivalent(pose), cam->get_projective_equivalent(pose2));
	  cv::Mat FundamentFuc(3, 3, CV_64FC1);
	  for (int i = 0; i < 3; i++)
	  {
		  for (int j = 0; j < 3; j++) {
			  FundamentFuc.at<double>(i, j) = Fundament(i, j);
		  }
	  }

	  cout << "@@@@@@@@@@" << endl;
	  cout << FundamentEPP << endl << FundamentFuc << endl << endl;
	  double rular = FundamentEPP.at<double>(0, 0) / FundamentFuc.at<double>(0, 0);
	  for (int i = 0; i < 3; i++)
	  {
		  for (int j = 0; j < 3; j++) {
			  cout << setprecision(15) << FundamentFuc.at<double>(i, j)*rular << endl;
		  }
	  }
	  cout << endl << endl;

	  cout << "分别输入图像1与图像2中对应的两个线段编号:";


	  ///////*选择所需线并显示*/
	  sort(keylines1.begin(), keylines1.end(), sortdes);
	  int limitMax = 10 < keylines1.size() ? 10 : keylines1.size();
	  for (int i = 0; i < limitMax; i++)
	  {
		  line(imageMat1, CvPoint(keylines1[i].startPointX, keylines1[i].startPointY), CvPoint(keylines1[i].endPointX, keylines1[i].endPointY), Scalar(0, 0, 255), 1);
		  string putNum = "1";
		  putNum[0] = '0' + i;
		  Point2d midPoint = CvPoint((keylines1[i].startPointX + keylines1[i].endPointX) / 2, (keylines1[i].startPointY + keylines1[i].endPointY) / 2);
		  putText(imageMat1, putNum, midPoint, CV_FONT_HERSHEY_SIMPLEX, 1, Scalar(0, 0, 255));
	  }
	  imshow("choosePic01", imageMat1);


	  limitMax = 10 < keylines2.size() ? 10 : keylines2.size();
	  sort(keylines2.begin(), keylines2.end(), sortdes);
	  for (int i = 0; i < limitMax; i++)
	  {
		  line(imageMat2, CvPoint(keylines2[i].startPointX, keylines2[i].startPointY), CvPoint(keylines2[i].endPointX, keylines2[i].endPointY), Scalar(0, 0, 255), 1);
		  string putNum = "1";
		  putNum[0] = '0' + i;
		  Point2d midPoint = CvPoint((keylines2[i].startPointX + keylines2[i].endPointX) / 2, (keylines2[i].startPointY + keylines2[i].endPointY) / 2);
		  putText(imageMat2, putNum, midPoint, CV_FONT_HERSHEY_SIMPLEX, 1, Scalar(0, 0, 255));
	  }
	  imshow("choosePic02", imageMat2);

	  waitKey(0);

	  int lineIdx1 = 0; //线1的代号
	  int lineIdx2[10]; //线2的代号
	  int lineTp = 0;
	  vector<Vec2> points_1, points_2;
	  fstream out(image_path + "point.txt", ios::out);
	  vector<Point2f> point1;
	  //char aa;
	  while (std::cin >> lineIdx1)
	  {
		  if (lineIdx1 == -1)
			  break;
		  cin >> lineIdx2[lineTp++];
		  //lineIdx2[lineTp + 1] = lineIdx2[lineTp];
		  //lineTp += 2;

		  point1.push_back(cvPoint(keylines1[lineIdx1].startPointX + pointOne.x, keylines1[lineIdx1].startPointY + pointOne.y));
		  point1.push_back(cvPoint(keylines1[lineIdx1].endPointX + pointOne.x, keylines1[lineIdx1].endPointY + pointOne.y));
		  cout << setprecision(15) << "##" << showpoint << keylines1[lineIdx1].startPointX + pointOne.x << " " << keylines1[lineIdx1].startPointY + pointOne.y << endl;
		  cout << setprecision(15) << "##" << showpoint << keylines1[lineIdx1].endPointX + pointOne.x << " " << keylines1[lineIdx1].endPointY + pointOne.y << endl;

		  out << setprecision(15) << keylines1[lineIdx1].startPointX + pointOne.x << " " << keylines1[lineIdx1].startPointY + pointOne.y << endl;
		  points_1.push_back(Vec2(keylines1[lineIdx1].startPointX + pointOne.x, keylines1[lineIdx1].startPointY + pointOne.y));
		  out << setprecision(15) << keylines1[lineIdx1].endPointX + pointOne.x << " " << keylines1[lineIdx1].endPointY + pointOne.y << endl;
		  points_1.push_back(Vec2(keylines1[lineIdx1].endPointX + pointOne.x, keylines1[lineIdx1].endPointY + pointOne.y));

		  //if (aa == '\n')
		  //	break;
	  }

	  destroyAllWindows();

	  //画图像1中的两点
	  cv::Mat imageMat4 = imread(img1, 1);
	  for (int ii = 0; ii < point1.size(); ii++)
	  {
		  circle(imageMat4, point1[ii], 1.5, Scalar(0, 0, 255), 3, 8, 0);
		  if (ii % 2)
			  line(imageMat4, point1[ii - 1], point1[ii], Scalar(255, 0, 0), 1);
	  }
	  imwrite(image_path + "picGai01.jpg", imageMat4);


	  //cout << point1 << endl;
	  //计算极线
	  vector<cv::Vec3f> corresEpilines;
	  //computeCorrespondEpilines(point1, 1, FundamentEPP, corresEpilines);
	  computeCorrespondEpilines(point1, 1, FundamentFuc, corresEpilines);

	  //画图2中的两点对应的两条极线
	  cv::Mat imageMat3 = imread(img2, 1);
	  lineTp = 0;
	  for (vector<cv::Vec3f>::const_iterator it = corresEpilines.begin(); it != corresEpilines.end(); ++it)
	  {
		  float a1 = (*it)[0];
		  float b1 = (*it)[1];
		  float c1 = (*it)[2];
		  // draw the epipolar line between first and last column
		  line(imageMat3, Point(0.0, -c1 / b1), Point(imageMat3.cols, -(c1 + a1 * imageMat3.cols) / b1), Scalar(0, 255, 0), 1.5);

		  //另一条线的a,b,c
		  int idx2 = lineIdx2[lineTp / 2];
		  lineTp++;
		  float a2 = keylines2[idx2].endPointY - keylines2[idx2].startPointY;    //y2-y1
		  float b2 = keylines2[idx2].startPointX - keylines2[idx2].endPointX;    //x1-x2
		  float c2 = (keylines2[idx2].endPointX + pointTwo.x)*(keylines2[idx2].startPointY + pointTwo.y) - (keylines2[idx2].startPointX + pointTwo.x)*(keylines2[idx2].endPointY + pointTwo.y);    //x2y1-x1y2

																																																   //画出竖直线
		  if (lineTp % 2 == 0)
		  {
			  line(imageMat3, Point(0.0, -c2 / b2), Point(imageMat3.cols, -(c2 + a2 * imageMat3.cols) / b2), Scalar(255, 0, 0), 0.5);
		  }

		  Point2d ans;
		  ans.x = (b1*c2 - b2*c1) / (a1*b2 - a2*b1);
		  ans.y = (a2*c1 - a1*c2) / (a1*b2 - a2*b1);
		  cout << setprecision(15) << "!!" << showpoint << ans << endl;
		  out << setprecision(15) << ans.x << " " << ans.y << endl;
		  points_2.push_back(Vec2(ans.x, ans.y));
		  circle(imageMat3, ans, 0.8, Scalar(0, 0, 255), 3, 8, 0);

	  }

	  waitKey(0);
	  imwrite(image_path + "picGai02.jpg", imageMat3);




	  std::vector<Vec3> points3D;
	  Vec3 Xtemp;
	  for (int i = 0; i < points_1.size(); i++) {
		  cout << "first camera's undistorted pixel coordinate:" << endl << cam->get_ud_pixel(points_1[i]) << endl;
		  cout << "second camera's undistorted pixel coordinate:" << endl << cam2->get_ud_pixel(points_2[i]) << endl;

		  //cout << "@@@@@@@@@@@@@@@@@@" << endl;
		  trianObj.add(cam->get_projective_equivalent(pose), cam->get_ud_pixel(points_1[i]));
		  trianObj.add(cam2->get_projective_equivalent(pose2), cam2->get_ud_pixel(points_2[i]));
		  //trianObj.add(cam->get_projective_equivalent(pose), points_1[i]);
		  //trianObj.add(cam2->get_projective_equivalent(pose2), points_2[i]);
		  Xtemp = trianObj.compute();
		  Xtemp = finalrt*Xtemp;
		  points3D.push_back(Xtemp);
		  trianObj.clear();
	  }

	  fstream outPoints(image_path + "points3D.txt", ios::out);
	  for (int i = 0; i < points3D.size(); i++) {
		  outPoints << points3D[i].x() << " " << points3D[i].y() << " " << points3D[i].z() << " " << 255 << " " << 0 << " " << 0 << endl;
	  }
	  outPoints.close();


	  //计算姿态角
	  for (unsigned int i = 0; i < points3D.size(); i = i + 2) {
		  getFuYang(points3D[i].x(), points3D[i].y(), points3D[i].z(), points3D[i + 1].x(), points3D[i + 1].y(), points3D[i + 1].z());
		  getShuiPing(points3D[i].x(), points3D[i].y(), points3D[i].z(), points3D[i + 1].x(), points3D[i + 1].y(), points3D[i + 1].z());
	  }

	  for (Landmarks::iterator iterL = my_sfm_data.structure.begin();
		  iterL != my_sfm_data.structure.end(); ++iterL)
	  {
		  // iterL->second.X = iterL->second.X - firstVec;
		  iterL->second.X = finalrt*(iterL->second.X);
		  //iterL->second.X.z() = -iterL->second.X.z();
		  //iterL->second.X.x() = -iterL->second.X.x();
	  }

	  for (Poses::iterator iterP = my_sfm_data.poses.begin();
		  iterP != my_sfm_data.poses.end(); ++iterP)
	  {
		  geometry::Pose3 & pose = iterP->second;
		  //pose = sim(pose);
		  //iterP->second.center().x() = -iterP->second.center().x();
		  //iterP->second.center().z() = -iterP->second.center().z();
		  // iterP->second.center() = iterP->second.center() - firstVec;
		  iterP->second.center() = finalrt*iterP->second.center();
	  }
	  getchar();
	  //-- Export to disk computed scene (data & visualizable results)
	  std::cout << "...Export SfM_Data to disk." << std::endl;
	  Save(my_sfm_data,
		  stlplus::create_filespec(sOutDir, "sfm_data", ".json"),
		  ESfM_Data(ALL));

	  Save(my_sfm_data,
		  stlplus::create_filespec(sOutDir, "cloud_and_poses", ".ply"),
		  ESfM_Data(ALL));
	  
	  getchar();
	  return EXIT_SUCCESS;
  }
  getchar();
  return EXIT_SUCCESS;
}
