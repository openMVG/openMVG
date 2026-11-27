// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/geometry/box.hpp"
#include "testing/testing.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

using namespace openMVG;
using namespace openMVG::geometry;

namespace {

struct OutputFile {
  const std::string path;
  explicit OutputFile(const std::string &name) : path(name) {}
  ~OutputFile() { std::remove(path.c_str()); }
};

Box MakeBox() {
  Box box;
  // Match the public Box vertex order, independently of half-plane solving.
  const Vec3 points[8] = {
    {0, 0, 4}, {0, 3, 4}, {2, 3, 4}, {2, 0, 4},
    {0, 0, 0}, {0, 3, 0}, {2, 3, 0}, {2, 0, 0}
  };
  std::copy(points, points + 8, box.points);
  return box;
}

bool ReadHeader(std::istream &stream, int &vertices, int &faces) {
  vertices = faces = 0;
  std::string line;
  if (!std::getline(stream, line) || line != "ply") return false;
  while (std::getline(stream, line)) {
    if (line == "end_header") return true;
    std::istringstream fields(line);
    std::string kind, name;
    int count;
    if (fields >> kind >> name >> count) {
      if (kind == "element" && name == "vertex") vertices = count;
      if (kind == "element" && name == "face") faces = count;
    }
  }
  return false;
}

}  // namespace

TEST(Box, PlyVertexRecords) {
  const Box box = MakeBox();
  const Vec3 colors[] = {Vec3(12, 34, 56), Vec3(0, 255, 0)};
  for (const Vec3 &color : colors) {
    const OutputFile output("box_vertex_records_test.ply");
    CHECK(Box::export_Ply(box, output.path, color));
    std::ifstream input(output.path);
    int vertices, faces;
    CHECK(ReadHeader(input, vertices, faces));
    EXPECT_EQ(8, vertices);
    for (int i = 0; i < vertices; ++i) {
      std::string line;
      CHECK(std::getline(input, line));
      std::istringstream fields(line);
      Vec3 point;
      std::array<int, 3> rgb;
      CHECK(fields >> point[0] >> point[1] >> point[2]
                   >> rgb[0] >> rgb[1] >> rgb[2]);
      std::string extra;
      CHECK(!(fields >> extra));
      EXPECT_MATRIX_NEAR(box.points[i], point, 0.0);
      for (int channel = 0; channel < 3; ++channel)
        EXPECT_EQ(static_cast<int>(color[channel]), rgb[channel]);
    }
  }
}

TEST(Box, PlyClosedSurface) {
  const Box box = MakeBox();
  const OutputFile output("box_closed_surface_test.ply");
  CHECK(Box::export_Ply(box, output.path));
  std::ifstream input(output.path);
  int vertices, faces;
  CHECK(ReadHeader(input, vertices, faces));
  EXPECT_EQ(8, vertices);
  EXPECT_EQ(12, faces);
  // Token parsing keeps the topology assertion independent of vertex lines.
  std::array<Vec3, 8> points;
  for (Vec3 &point : points) {
    int red, green, blue;
    CHECK(input >> point[0] >> point[1] >> point[2] >> red >> green >> blue);
  }
  std::map<std::pair<int, int>, std::pair<int, int>> edges;
  double volume = 0.0;
  double area = 0.0;
  const Vec3 center(1, 1.5, 2);
  for (int face = 0; face < faces; ++face) {
    int count;
    std::array<int, 3> indices;
    CHECK(input >> count >> indices[0] >> indices[1] >> indices[2]);
    EXPECT_EQ(3, count);
    for (int index : indices) CHECK(index >= 0 && index < vertices);
    const Vec3 &a = points[indices[0]];
    const Vec3 &b = points[indices[1]];
    const Vec3 &c = points[indices[2]];
    const Vec3 normal = (b - a).cross(c - a);
    CHECK(normal.dot((a + b + c) / 3.0 - center) > 0.0);
    volume += a.dot(b.cross(c)) / 6.0;
    area += normal.norm() / 2.0;
    for (int i = 0; i < 3; ++i) {
      const int u = indices[i], v = indices[(i + 1) % 3];
      auto &edge = edges[std::minmax(u, v)];
      ++edge.first;
      edge.second += u < v ? 1 : -1;
    }
  }
  EXPECT_NEAR(24.0, volume, 1e-12);
  EXPECT_NEAR(52.0, area, 1e-12);
  EXPECT_EQ(18, edges.size());
  for (const auto &edge : edges) {
    EXPECT_EQ(2, edge.second.first);
    EXPECT_EQ(0, edge.second.second);
  }
}

int main() { TestResult tr; return TestRegistry::runAllTests(tr); }
