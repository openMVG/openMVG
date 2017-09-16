// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "split.hpp"

#include "testing/testing.h"

#include <string>
#include <vector>

TEST(Split, stringEmpty)
{
  const std::string sInput = "";
  const char sDelimiter = ' ';
  std::vector<std::string> vec_str;
  EXPECT_FALSE( stl::split(sInput, sDelimiter, vec_str));
  EXPECT_EQ( 0, vec_str.size() );
}

TEST(Split, delimiterNotExist)
{
  const std::string sInput = "A string";
  const char sDelimiter = '*';
  std::vector<std::string> vec_str;
  EXPECT_FALSE( stl::split(sInput, sDelimiter, vec_str));
  EXPECT_EQ( 1, vec_str.size() );
}

TEST(Split, delimiterExist)
{
  const std::string sInput = "A string";
  const char sDelimiter = ' ';
  std::vector<std::string> vec_str;
  EXPECT_TRUE( stl::split(sInput, sDelimiter, vec_str));
  EXPECT_EQ( 2, vec_str.size() );
}

TEST(Split, stringSplit3part)
{
  const std::string sInput = "A string useless";
  const char sDelimiter = ' ';
  std::vector<std::string> vec_str;
  EXPECT_TRUE( stl::split(sInput, sDelimiter, vec_str));
  EXPECT_EQ( 3, vec_str.size() );
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
