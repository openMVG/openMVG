#pragma once

#include <cmath>

inline float getCoterminalAngle(float angle)
{
  const float f2PI = 2.f*M_PI;
  while (angle > f2PI) {
    angle -= f2PI;
  }
  while (angle < 0.0f) {
    angle += f2PI;
  }
  return angle;
}


struct DeepClassifierKeypoint {
	DeepClassifierKeypoint() :
		x(0.0f),
		y(0.0f),
		l1(0.0f),
		l2(0.0f),
		phi(0.0f),
		a(0.0f),
		b(0.0f),
		c(0.0f) {
	}
	DeepClassifierKeypoint(
			float _x = 0.0f, float _y = 0.0f,
			float _l1 = 0.0f, float _l2 = 0.0f,
			float _phi = 0.0f,
			float _a = 0.0f, float _b = 0.0f, float _c = 0.0f) :
		x(_x),
		y(_y),
		l1(_l1),
		l2(_l2),
		phi(_phi),
		a(_a),
		b(_b),
		c(_c) {

		}
	
	DeepClassifierKeypoint(
			float _x = 0.0f, float _y = 0.0f,
			float _a = 0.0f, float _b = 0.0f, float _c = 0.0f) :
		x(_x),
		y(_y),
		a(_a),
		b(_b),
		c(_c) {
			l1 = (a + c - std::sqrt(a*a + c*c + 4 * b*b - 2 * a*c)) / 2;
			l2 = (a + c + std::sqrt(a*a + c*c + 4 * b*b - 2 * a*c)) / 2;
			l1 = 1.0 / std::sqrt(l1);
			l2 = 1.0 / std::sqrt(l2);

			phi = 0.0;
			if (b == 0)
			{
			  if (a > c)
				phi = M_PI / 2; // else 0
			}
			else
			{
			  const double t = std::atan(2 * b / (a - c));

			  if (a < c)
				phi = t / 2;
			  else
				phi = t / 2 + ((b > 0) ? -M_PI / 2 : M_PI / 2);
			}

			if (l1 > l2)
			{
			  std::swap(l1, l2);
			  phi = getCoterminalAngle(M_PI / 2 - phi);
			}
		}

	float x, y;
 	float l1, l2;
	float phi; 
	float a, b, c;
};
