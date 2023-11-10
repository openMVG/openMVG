// Copyright (c) 2010 libmv authors.
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

// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/test_data_sets.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/projection.hpp"

#include <cmath>
#include <fstream>

namespace openMVG {


nViewDatasetConfigurator::nViewDatasetConfigurator(int fx, int fy,
  int cx, int cy, double distance, double jitter_amount):
  _fx(fx), _fy(fy), _cx(cx), _cy(cy), _dist(distance),
  _jitter_amount(jitter_amount)
{}

NViewDataSet NRealisticCamerasRing(size_t nviews, size_t npoints,
                                   const nViewDatasetConfigurator & config)
{
  //-- Setup a camera circle rig.
  NViewDataSet d;
  d._n = nviews;
  d._K.resize(nviews);
  d._R.resize(nviews);
  d._t.resize(nviews);
  d._C.resize(nviews);
  d._x.resize(nviews);
  d._x_ids.resize(nviews);

  d._X.resize(3, npoints);
  d._X.setRandom();
  d._X *= 0.6;

  Vecu all_point_ids(npoints);
  for (size_t j = 0; j < npoints; ++j)
    all_point_ids[j] = j;

  for (size_t i = 0; i < nviews; ++i) {
    Vec3 camera_center, t, jitter, lookdir;

    const double theta = i * 2 * M_PI / nviews;
    //-- Circle equation
    camera_center << sin(theta), 0.0, cos(theta); // Y axis UP
    camera_center *= config._dist;
    d._C[i] = camera_center;

    jitter.setRandom();
    jitter *= config._jitter_amount / camera_center.norm();
    lookdir = -camera_center + jitter;

    d._K[i] << config._fx,           0, config._cx,
                        0,  config._fy, config._cy,
                        0,           0,          1;
    d._R[i] = LookAt(lookdir);  // Y axis UP
    d._t[i] = -d._R[i] * camera_center; // [t]=[-RC] Cf HZ.
    d._x[i] = Project(d.P(i), d._X);
    d._x_ids[i] = all_point_ids;
  }
  return d;
}

Mat34 NViewDataSet::P(size_t i)const {
  assert(i < _n);
  Mat34 P;
  P_From_KRt(_K[i], _R[i], _t[i], &P);
  return P;
}

Mat34 NViewDataSet::Rt(size_t i)const {
  assert(i < _n);
  Mat34 P;
  return HStack(_R[i],_t[i]);
}

void NViewDataSet::ExportToPLY(
  const std::string & out_file_name)const {
  std::ofstream outfile(out_file_name.c_str(), std::ios_base::out);
  if (outfile) {
    outfile << "ply"
     << std::endl << "format ascii 1.0"
     << std::endl << "comment NViewDataSet export"
     << std::endl << "comment It shows 3D point structure and cameras"
                  << "+ camera looking direction"
     << std::endl << "element vertex " << _X.cols() + _t.size()*2
     << std::endl << "property float x"
     << std::endl << "property float y"
     << std::endl << "property float z"
     << std::endl << "property uchar red"
     << std::endl << "property uchar green"
     << std::endl << "property uchar blue"
     << std::endl << "end_header" << std::endl;

    //-- Export 3D point cloud
    for (Mat3X::Index i = 0; i < _X.cols(); ++i) {
      // Exports the point position and point color
      outfile << _X.col(i).transpose()
        << " " << "255 255 255" << std::endl;
    }

    //-- Export 3D camera position t = -RC
    for (size_t i = 0; i < _t.size(); ++i) {
      // Exports the camera position and camera color
      outfile << (-_R[i].transpose()*_t[i]).transpose()
        << " " << "0 255 0" << std::endl;
    }
    for (size_t i = 0; i < _t.size(); ++i) {
      Vec3 test;
      test << 0, 0 , 0.4;
      // Exports the camera normal
      outfile << ((-_R[i].transpose()*_t[i])+
        (_R[i].transpose()*test)).transpose()
        << " " << "255 0 0" << std::endl;
    }
    outfile.close();
  }
}

NViewDataSet NRealisticCamerasCardioid(size_t nviews, size_t npoints,
                                        const nViewDatasetConfigurator & config)
{
  //-- Setup a camera circle rig.
  NViewDataSet d;
  d._n = nviews;
  d._K.resize(nviews);
  d._R.resize(nviews);
  d._t.resize(nviews);
  d._C.resize(nviews);
  d._x.resize(nviews);
  d._x_ids.resize(nviews);

  d._X.resize(3, npoints);
  d._X.setRandom();
  d._X *= 0.6;

  Vecu all_point_ids(npoints);
  for (size_t j = 0; j < npoints; ++j)
    all_point_ids[j] = j;

  for (size_t i = 0; i < nviews; ++i) {
    Vec3 camera_center, t, jitter, lookdir;

    const double theta = i * 2 * M_PI / nviews;
    //-- Cardioid
    camera_center <<
      2*sin(theta)-(sin(2*theta)),
      0.0,
      2*cos(theta)-(cos(2*theta)); // Y axis UP
    camera_center *= config._dist;
    d._C[i] = camera_center;

    jitter.setRandom();
    jitter *= config._jitter_amount / camera_center.norm();
    lookdir = -camera_center + jitter;

    d._K[i] << config._fx,           0, config._cx,
      0,  config._fy, config._cy,
      0,           0,          1;
    d._R[i] = LookAt(lookdir);  // Y axis UP
    d._t[i] = -d._R[i] * camera_center; // [t]=[-RC] Cf HZ.
    d._x[i] = Project(d.P(i), d._X);
    d._x_ids[i] = all_point_ids;
  }
  return d;
}

// -----------------------------------------------------------------------------

const double K_[2][3] = {
  {2584.9325098195013197, 0, 249.77137587221417903},
  {0, 2584.9325098195013197, 278.31267937919352562}
  //  0 0 1 
};

static constexpr unsigned synth_nviews_ = 4;
static constexpr unsigned synth_npts_ = 5;

// camera format: just like a 3x4 [R|T] but transposed to better fit row-major:
// | R |
// | - |
// | T'|
//
// In this case, this data is from the synthcurves multiview dataset,
// so that instead of T, C is stored:
// | R |
// | - |
// | C'|
const double
cameras_gt_[synth_nviews_][4][3] = {
// extrinsics for frame 42
{
{0.032488343069021832776, 0.14118885304658673752, 0.98944945062394118462},
{-0.28679990702465507635, 0.94965585180689460199, -0.12609352267095536027},
{-0.95743946069465868387, -0.27967744082122686367, 0.071345695037685272211},
{1071.484198049582119, 320.76850873549409471, -85.986368935484179588}
},
// extrinsics for frame 54
{
{0.47790123270283219048, -0.019299093028001035321, -0.87820154679288175981},
{-0.29465811799137531235, -0.94535498497159320408, -0.13957272617219795841},
{-0.82751858304384573461, 0.32547119288449349872, -0.45747294709026314896},
{925.05253488236451176, -384.00318581806010343, 531.87270782843597772}
},
// extrinsics for frame 62
{
{-0.219610930168113061, -0.33532970351451174551, -0.91614683828061438398},
{-0.50550514340331198504, -0.76406316361831605466, 0.40083915975658834796},
{-0.83440732819378671259, 0.55114559958548192675, -0.0017142677930838123856},
{951.80671923514557875, -619.60363267357240602, 4.2313312789905133116}
},
// extrinsics for frame 07
{
{0.91074869806248703874, 0.30990046893594913602, -0.27294414873882894002},
{0.29070436143147293517, -0.95055540259185944407, -0.10924926017208461126},
{-0.29330493214776437449, 0.020152567000437736355, -0.95580651327613819213},
{342.08616590607340413, -28.455793406432913883, 1067.8273738052921544}
}
};

// Input points and tangents extracted from the synthcurves dataset:
//
// http://github.com/rfabbri/synthcurves-multiview-3d-dataset
//
// sub-dataset:
// spherical-ascii-100_views-perturb-radius_sigma10-normal_sigma0_01rad-minsep_15deg-no_two_cams_colinear_with_object-aspect_ratio1
//      
//
// Frame files: 42, 54, 62, 07
//
// Points and tangents:
// 
// 620
// 1009
// 3011
// 3389
// 0-based ids. +1 to get file line
// 
// Extracting this data from those files: use scripts/getlines from minus
// or just by hand
//
// This is in pixel image coordinates

const double
p_gt_[synth_nviews_][synth_npts_][2] = {
// 2D points for frame 42
{                                                // For debugging: versions with K inverse applied (bearings):
{181.53712861181011817, 382.84265753788542952},  //  -0.0263969163609492629696   0.0404381846572818004493   
{361.39404494701216208, 353.17104859076766843},  //    0.04318204388345608935     0.0289595062645567058457  
{342.08123422244102585, 137.87448982117351193},  //    0.0357107421565418803322  -0.0543295382082630978759  
{243.09524202092040923, 298.35373638828008325},  //   -0.0025827110866271374423   0.0077530291150564589753 
{285.88966212823157775, 251.48973104783391364}   //    0.0139726225419090216429  -0.0103766532508938053025 
},                                                                                                          
// 2D points for frame 54                                                                                   
{                                                                                                           
{320.61083811349851658, 199.23585641086629039},  //    0.027404762782851491143   -0.0305914458764143665226 
{177.83962742245475397, 163.63860158201131867},  //   -0.0278273216714591126175  -0.0443625036095002672765 
{225.74947198803425863, 316.24347634112569949},  //   -0.0092930487712645504228   0.0146738055318050841791 
{210.92593414503781446, 312.1127002295278885},   //   -0.0150276425321018564096   0.0130757846566348229222 
{247.68819285385683315, 263.17278766241219046}   //   -0.0008058945486753843479  -0.0058569775648953104064 
},                                                                                                          
// 2D points for frame 62                                                                                   
{                                                                                                           
{330.38724148135435144, 234.16270784506826885},  //    0.0311868357502182447227  -0.0170797385875301327429 
{165.33058499047140799, 291.56955369118014687},  //   -0.0326665359969645907601   0.0051285185441504316239 
{199.86675126597054941, 393.58510880586692338},  //   -0.0193059681119985282471   0.0445939803026898229366 
{313.99181820108196916, 389.73404770358069982},  //    0.0248441466401581739776   0.0431041692195543724164 
{248.50928647922813752, 333.51852292954941959}   //   -0.0004882484893480448784   0.0213567833359837866425 
},                                                                                                          
// 2D points for frame 07                                                                                   
{                                                                                                           
{195.98098993753490049, 156.78341685173867859},  //   -0.0208092032307780894218  -0.0470144818349400167579 
{164.52669179292479384, 134.9377776538624687},   //   -0.0329775279453008857145  -0.0554656267351995718728 
{197.46123952507306853, 358.3286736009371225},   //   -0.0202365578785628635883   0.0309547711275954817722 
{52.110592924391973213, 218.51730104909975694},  //   -0.076466515932992112914   -0.0231322783488336070068 
{177.00162990735412905, 256.9990394819891435}    //   -0.0281515071238518932439  -0.008245337089552343124  
}
};

const double 
t_gt_[synth_nviews_][synth_npts_][2] = {
// 2D tangents for frame 42
{                                                       // For debugging: versions with K iverse + normalization 
{0.99164213923671906681, -0.12901886563609127334},     //    0.9916421392367189557859  -0.129018865636091273341 
{-0.27407482149753986667, 0.96170837171207557148},     //   -0.274074821497539866666    0.9617083717120755714802
{-0.99085473006575885968, 0.1349329607853933799},      //   -0.990854730065758637636    0.1349329607853933521433
{0.21688457417910370073, -0.97619725541672619507},     //    0.216884574179103700731   -0.9761972554167261950653
{-0.88923826058271226991, 0.45744433094731007383}      //   -0.8892382605827123809306   0.4574443309473101293428
},
// 2D tangents for frame 54
{
{-0.98462682094023079582, -0.1746711867628284176},     //   -0.9846268209402306847977  -0.174671186762828389849 
{-0.80826947668316584394, 0.58881274872604572046},     //   -0.808269476683165732922    0.58881274872604572046  
{0.90119839239834154121, 0.43340680375213874731},      //    0.9011983923983414301873   0.4334068037521386917987
{-0.54078671415378842813, 0.84115975283815680452},     //   -0.5407867141537883171054   0.8411597528381566934996
{0.99935984779679032375, 0.035775614761675351982}      //    0.999359847796790545793    0.0357756147616753589213
},
// 2D tangents for frame 62
{
{-0.91615034676027939931, 0.4008348065363335766},      //   -0.9161503467602793993052   0.4008348065363336321099
{0.44184360558182206313, 0.89709209572175763192},      //    0.4418436055818220076219   0.897092095721757631921 
{0.99085696151653412933, -0.13491657353424593713},      //    0.9908569615165341293306  -0.1349165735342459371271,
{0.17976067146968616184, 0.98371037454769549857},      //    0.1797606714696861618386   0.9837103745476956095928
{0.91927766178284631149, -0.39360968045395255954}      //    0.9192776617828464225113  -0.3936096804539526150535
},
// 2D tangents for frame 07
{
{-0.88483975537184744731, -0.46589548969000454948},    //   -0.8848397553718474473072  -0.4658954896900045494768
{-0.95661607022265571221, -0.29135149594907394643},    //   -0.9566160702226556011851  -0.2913514959490739464343
{-0.67890121232474576196, 0.73422962614157050165},     //   -0.6789012123247458729836   0.7342296261415705016518
{-0.91684145433395725089, 0.39925148417356498554},     //   -0.916841454333957250888    0.3992514841735649855359
{-0.017735947168850744321, -0.99984270571826627805}    //   -0.017735947168850744321   -0.9998427057182662780477
}
};

const double pts3d_gt_[synth_npts_][3] = {
{-39.999960000000001514, 40.000016999999999712, -39.999979999999993652},
{-28.799995999999886465, 40.000010000000003174, 40.000010000000003174},
{16.241229516856364512, -45.185185185185176238, 41.368080573302677294},
{-83.024179089510298013, -7.2456979436932478222, -4.412526863075626693},
{-15.900733289698191442, -6.9926202045388530237, 12.583214033874593696}
};

const double t3d_gt_[synth_npts_][3] = {
{0, 0, 1},
{-1, 0, 0},
{-0.34011103186525448727, -0.10551104075352309153, -0.93444738015720296698},
{-0.71448475613149997621, -0.66437364477423688225, 0.21936087478196092393},
{-0.38185692861567399614, 0.23333310127926898403, -0.89428236587534382096}
};

// number of points is hardcoded and number of views is hardcoded
void 
NOrientedPointsCamerasSphere(size_t nviews, size_t npoints, NViewOrientedDataSet *dp, nViewDatasetConfigurator *conf)
{
  //-- Setup a camera rig
  // based on github.com/rfabbri/synthcurves-multiview-3d-dataset
  // by hardcoding points
  NViewOrientedDataSet &d = *dp;

  assert(nviews <= synth_nviews_ && npoints <= synth_npts_);

  d._n = nviews;
  d._K.resize(nviews);
  d._R.resize(nviews);
  d._t.resize(nviews);
  d._C.resize(nviews);
  d._x.resize(nviews);
  d._tgt2d.resize(nviews);
  d._x_ids.resize(nviews);
  d._X.resize(3, npoints);
  d._Tgt3d.resize(3, npoints);

  Vecu all_point_ids(npoints);
  for (size_t p = 0; p < npoints; ++p) {
    all_point_ids[p] = p;
    d._X.col(p) << pts3d_gt_[p][0], pts3d_gt_[p][1], pts3d_gt_[p][2];
    d._Tgt3d.col(p) << t3d_gt_[p][0], t3d_gt_[p][1], t3d_gt_[p][2];
  }
  d._K_raw = K_;
  d._cameras_gt_raw = cameras_gt_;

  for (size_t v = 0; v < nviews; ++v) {
    d._C[v] << cameras_gt_[v][3][0] , cameras_gt_[v][3][1] , cameras_gt_[v][3][2];

    d._K[v] << K_[0][0],           K_[0][1], K_[0][2],
               K_[1][0],           K_[1][1], K_[1][2],
                      0,           0,          1;
    d._R[v] << cameras_gt_[v][0][0] , cameras_gt_[v][0][1] , cameras_gt_[v][0][2],
               cameras_gt_[v][1][0] , cameras_gt_[v][1][1] , cameras_gt_[v][1][2],
               cameras_gt_[v][2][0] , cameras_gt_[v][2][1] , cameras_gt_[v][2][2];

    d._t[v] = -d._R[v] * d._C[v];
    d._x[v].resize(2,npoints);
    d._tgt2d[v].resize(2,npoints);
    for (unsigned p = 0; p < npoints; ++p) {
      d._x[v].col(p) << p_gt_[v][p][0], p_gt_[v][p][1];
      d._tgt2d[v].col(p) << t_gt_[v][p][0], t_gt_[v][p][1];
    }
    d._x_ids[v] = all_point_ids;
  }
  assert (fabs(K_[0][0] - K_[1][1]) < 1e-6);

  conf->_fx = K_[0][0]; conf->_fy = K_[1][1]; conf->_cx = K_[0][2]; conf->_cy = K_[1][2];
  conf->_dist = conf->_jitter_amount = 0; // not used
}


}  // namespace openMVG
