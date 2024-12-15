
#include <sfm_gsd.hpp>

namespace openMVG {
namespace sfm {


bool AreSame(double a, double b){ return std::fabs(a - b) <= std::numeric_limits<double>::epsilon(); }

double get_ground_sampling_distance(IndexT view_id, const SfM_Data & sfm_data, const Vec2 & sensor_size){
    
    int counter = 0;
    double average_camera_height = 0.0;

    for(auto part: sfm_data.GetPoses()){
        average_camera_height += (part.second.center()[2] - average_camera_height)/(counter+1);
        counter++;
    }

    counter = 0;
    double average_ground_height = 0.0;
    for(auto part: sfm_data.GetLandmarks()){
        average_ground_height += (part.second.X(2) - average_ground_height)/(counter+1);
        counter++;
    }

    double average_distance_to_ground = abs(average_camera_height-average_ground_height);

    openMVG::sfm::View *view = sfm_data.GetViews()->at(view_id);
    std::shared_ptr<cameras::IntrinsicBase> intrinsic = sfm_data.GetIntrinsics().at(view->id_intrinsic);
    std::vector<double> params = intrinsic->getParams();

    unsigned int width = intrinsic->w();
    unsigned int height= intrinsic->h();
    
    double focal_mm = params[0] * sensor_size[0]/width;

    // GSD is the maximum value of the horizontal and vertical GSDs
    double GSD_H = -1.0;
    double GSD_W = -1.0;

    if(!AreSame(0.0,sensor_size[0])){
        GSD_H = (average_distance_to_ground * sensor_size[0]) / (focal_mm * height);
    }
    if(!AreSame(0.0,sensor_size[1])){
        GSD_W = (average_distance_to_ground * sensor_size[1]) / (focal_mm * width);
    }
    
    return std::max(GSD_H,GSD_W);
}

double euclidean_distance(Vec2 v1, Vec2 v2){
    return sqrt(  pow(v1[0]-v2[0],2) + pow(v1[1]-v2[1])  );
}

std::vector<Vec2> getRealGroundCoordinates(IndexT view_id, const SfM_Data & sfm_data, const double & average_height){
    openMVG::sfm::View *view = sfm_data.GetViews()->at(view_id);

    if(!sfm_data.IsPoseAndIntrinsicDefined(view)){
        return std::vector<Vec2>(0);
    }

    geometry::Pose3 pose = sfm_data.poses.at(view->id_pose);

    double camera_height = pose.center()[2];

    std::shared_ptr<cameras::IntrinsicBase> intrinsic = sfm_data.GetIntrinsics().at(view->id_intrinsic);

    unsigned int width = intrinsic->w();
    unsigned int height= intrinsic->h();

    std::vector<Vec2> image_points(4);
    image_points[0] = Vec2(0.0,0.0);
    image_points[1] = Vec2(width-1.0,0.0);
    image_points[2] = Vec2(width-1.0,height-1.0);
    image_points[3] = Vec2(0.0,height-1.0);
    
    std::vector<Vec2> ground_points;

    std::vector<double> params = intrinsic->getParams();

    double cu = params[1];
    double cv = params[2];
    double fu = params[0];
    double fv = params[0];

    for(auto p: image_points){
        // backproject the point
        // getting the ray
        Vec2 kp;
        kp[0] = (kp[0] - cu) / fu;
        kp[1] = (kp[1] - cv) / fv;
        
        Vec2 un_kp= intrinsic->remove_disto(kp);

        const double rho2_d = kp[0] * kp[0] + kp[1] * kp[1];
        const double tmpD = std::max(1 + (1 - fu*fu) * rho2_d, 0.0);

        Vec3 cray;
        cray[0] = kp[0];
        cray[1] = kp[1];
        cray[2] = 1 - fu * (rho2_d + 1) / (fu + sqrt(tmpD));

        Mat3 RT = pose.rotation().transpose() * pose.center();

        // compute the scale
        double scale = -(camera_height-average_height) / RT[2][0];

        // compute ground position
        Vec2 gp = Vec2(pose.center()[0]+scale,pose.center()[1]+scale);

        ground_points.push_back(gp);
    }

    return ground_points;
}

}
}