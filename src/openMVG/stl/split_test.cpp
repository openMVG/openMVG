#include "testing/testing.h"

#include "split.hpp"

#include <string>
#include <vector>

TEST(Split, stringEmpty)
{
  std::string sInput = "";
  std::string sDelimiter = " ";
  std::vector<std::string> vec_str;
  EXPECT_FALSE( split( sInput, sDelimiter, vec_str ) );
  EXPECT_EQ( 1, vec_str.size() );
}

TEST(Split, delimiterEmpty)
{
  std::string sInput = "A string";
  std::string sDelimiter = "";
  std::vector<std::string> vec_str;
  EXPECT_FALSE( split( sInput, sDelimiter, vec_str ) );
  EXPECT_EQ( 0, vec_str.size() );
}


TEST(Split, delimiterNotExist)
{
  std::string sInput = "A string";
  std::string sDelimiter = "_";
  std::vector<std::string> vec_str;
  EXPECT_FALSE( split( sInput, sDelimiter, vec_str ) );
  EXPECT_EQ( 1, vec_str.size() );
}

TEST(Split, delimiterExist)
{
  std::string sInput = "A string";
  std::string sDelimiter = " ";
  std::vector<std::string> vec_str;
  EXPECT_TRUE( split( sInput, sDelimiter, vec_str ) );
  EXPECT_EQ( 2, vec_str.size() );
}

TEST(Split, stringSplit3part)
{
  std::string sInput = "A string useless";
  std::string sDelimiter = " ";
  std::vector<std::string> vec_str;
  EXPECT_TRUE( split( sInput, sDelimiter, vec_str ) );
  EXPECT_EQ( 3, vec_str.size() );
}

TEST(Split, ontheSameString)
{
	std::string sInput = "";
	std::string sDelimiter = ";";
	std::vector<std::string> vec_str;
	vec_str.push_back("foo;");
	EXPECT_TRUE( split( vec_str[0], sDelimiter, vec_str ) );
	EXPECT_EQ( 2, vec_str.size() );
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

