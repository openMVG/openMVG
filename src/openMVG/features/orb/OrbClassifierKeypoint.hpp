#ifndef ORB_CLASSIFIER_KEYPOINT_H
#define ORB_CLASSIFIER_KEYPOINT_H

struct OrbClassifierKeypoint {
    OrbClassifierKeypoint() {
       x = 0.0f;
       y = 0.0f;
       angle = 0.0f;
       size = 0.0f;
    }
    OrbClassifierKeypoint(float _x, float _y, float _angle, float _size) 
        : x(_x),
          y(_y),
          angle(_angle),
          size(_size) {}
    float x, y;
    float angle;
    float size;
};

#endif
