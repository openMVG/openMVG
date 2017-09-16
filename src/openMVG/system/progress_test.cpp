// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (C) 2017 Bjorn Piltz Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "third_party/progress/progress.hpp"

#include "testing/testing.h"

#include <sstream>

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

static const int Count = 100*1000;
static const int num_threads = 10;

// This funcion counts to Count
int singleThreadedCount(C_Progress * progress=nullptr)
{
  if (!progress)
      progress = &C_Progress::dummy();
  int result = 0;
  for (int i = 0; i<Count; i++)
  {
    if (progress->hasBeenCanceled())
      break;

    ++result;
    ++(*progress);
  }
  return result;
}

// This function counts to Count
int multiThreadedCount(C_Progress * progress=nullptr)
{
  if (!progress)
      progress = &C_Progress::dummy();
  int result = 0;

#ifdef OPENMVG_USE_OPENMP
  omp_set_num_threads(num_threads);
  #pragma omp parallel for
#endif
  for (int i = 0; i<Count; i++)
  {
    if (progress->hasBeenCanceled())
      // We are not allowed to use 'break' in an omp loop, 
      continue;

#ifdef OPENMVG_USE_OPENMP
#pragma omp critical
#endif
    {
      ++result;
      ++(*progress);
    }
  }
  return result;
}

class MockProgress : public C_Progress
{
public:
  virtual void restart(unsigned long ulExpected_count, const std::string& msg=std::string())override
  {
    this->ulExpected_count = ulExpected_count;
    this->currentCount = 0;
  }
  virtual bool hasBeenCanceled()const override
  {
    return false;
  }
  virtual void inc_tic() override
  {
    currentCount++;
  }
  unsigned long ulExpected_count = 0, currentCount = 0;
};

class CancelProgress : public MockProgress
{
  virtual bool hasBeenCanceled()const override
  {
    return true;
  }
};

template<int CancelAtNumber>
class CancelAt : public MockProgress
{
  virtual bool hasBeenCanceled()const override
  {
    return currentCount==CancelAtNumber;
  }
};

TEST(Progress, dummy)
{
  // Make sure the singleThreadedCount() works with the default dummy progress
  EXPECT_EQ(Count, singleThreadedCount());

  // Make sure the multiThreadedCount() works with the default dummy progress
  EXPECT_EQ(Count, multiThreadedCount());

  MockProgress progress;
  EXPECT_EQ(Count, singleThreadedCount(&progress));
  EXPECT_EQ(Count, progress.currentCount);

  progress.restart(Count);
  EXPECT_EQ(Count, multiThreadedCount(&progress));
  EXPECT_EQ(Count, progress.currentCount);
}

TEST(Progress, cancel)
{
  // Make sure the singleThreadedCount() works with MockProgress
  CancelProgress progress;
  EXPECT_EQ(0, singleThreadedCount(&progress));

  progress.restart(Count);
  EXPECT_EQ(0, multiThreadedCount(&progress));

  CancelAt<13> cancelAtThirteen;
  EXPECT_EQ(13, singleThreadedCount(&cancelAtThirteen));
}

// Put this include here, because C_Progress_display is not needed earlier
#include "third_party/progress/progress_display.hpp"

TEST(Progress, terminalProgress)
{
  // Test C_Progress_display:
  const std::string expectedOutPut =
    "\nCounting electric sheep\n"
    "0%   10   20   30   40   50   60   70   80   90   100%\n"
    "|----|----|----|----|----|----|----|----|----|----|\n"
    "***************************************************\n";

  std::ostringstream out;

  // Make sure the singleThreadedCount() works with MockProgress
  C_Progress_display progress(Count, out, "\nCounting electric sheep\n");
  EXPECT_EQ(Count, singleThreadedCount(&progress));

  EXPECT_EQ(expectedOutPut, out.str());
  out.seekp(0);

  progress.restart(Count);
  EXPECT_EQ(Count, multiThreadedCount(&progress));
  EXPECT_EQ(expectedOutPut, out.str());
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr); }
/* ************************************************************************* */
