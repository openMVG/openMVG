#include "testing/testing.h"

#include "split.hpp"

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

TEST(Split, ontheSameString)
{
  const char sDelimiter = ';';
  std::vector<std::string> vec_str = {"foo;"};
  EXPECT_TRUE( stl::split(vec_str[0], sDelimiter, vec_str));
  EXPECT_EQ( 1, vec_str.size() );
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

