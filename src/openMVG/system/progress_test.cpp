// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (C) 2017 Bjorn Piltz, Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/system/progressinterface.hpp"
#include "openMVG/system/loggerprogress.hpp"

#include "testing/testing.h"

#include <set>

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

static const int Count = 100*1000;
static const int num_threads = 10;

using openMVG::system::ProgressInterface;

// This function counts to Count
int singleThreadedCount(ProgressInterface * progress=nullptr)
{
  if (!progress)
    progress = &ProgressInterface::dummy();
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
int multiThreadedCount(ProgressInterface * progress=nullptr)
{
  if (!progress)
    progress = &ProgressInterface::dummy();
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

class MockProgress : public ProgressInterface
{
public:
  virtual void Restart(const std::uint32_t expected_count, const std::string& msg = {})override
  {
    ProgressInterface::Restart(expected_count, msg);
    this->ulExpected_count = expected_count;
    this->currentCount = 0;
  }
  virtual bool hasBeenCanceled()const override
  {
    return false;
  }
  virtual std::uint32_t operator+=(const std::uint32_t increment) override
  {
    ProgressInterface::operator+=(increment);
    ++currentCount;
    return currentCount;
  }
  std::uint32_t ulExpected_count = 0, currentCount = 0;
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

  progress.Restart(Count);
  EXPECT_EQ(Count, multiThreadedCount(&progress));
  EXPECT_EQ(Count, progress.currentCount);
}

TEST(Progress, cancel)
{
  // Make sure the singleThreadedCount() works with MockProgress
  CancelProgress progress;
  EXPECT_EQ(0, singleThreadedCount(&progress));

  progress.Restart(Count);
  EXPECT_EQ(0, multiThreadedCount(&progress));

  CancelAt<13> cancelAtThirteen;
  EXPECT_EQ(13, singleThreadedCount(&cancelAtThirteen));
}

TEST(Progress, percent)
{
  ProgressInterface progress(10);
  EXPECT_EQ(10, progress.expected_count());
  EXPECT_EQ(0, progress.count());
  EXPECT_EQ(0, progress.Percent());

  progress += 2;
  EXPECT_EQ(2, progress.count());
  EXPECT_EQ(20, progress.Percent());
  progress += 4;
  EXPECT_EQ(6, progress.count());
  EXPECT_EQ(60, progress.Percent());
  progress += 4;
  EXPECT_EQ(10, progress.count());
  EXPECT_EQ(100, progress.Percent());
}

using openMVG::system::LoggerProgress;

TEST(LogProgress, logging)
{
  LoggerProgress progress(10);
  for (int i = 0; i < 10; ++i)
  {
    ++progress;
  }
}

TEST(LogProgress, logging_if)
{
  LoggerProgress progress(10);
  for (int i = 0; i < 10; ++i)
  {
    const bool b_update_progress_required = progress.Increment(1);
    OPENMVG_LOG_INFO_IF(b_update_progress_required) << progress.PercentString();
  }
}

// Check that every pourcentage is displayed once
TEST(LogProgress, conditional_logging)
{
  LoggerProgress progress(800, "", 2);
  std::set<int> previously_displayed_percentage;
  for (int i = 0; i < progress.count(); ++i)
  {
    if (progress.Increment(1)) // If display is triggered
    {
      EXPECT_EQ(0, previously_displayed_percentage.count(progress.Percent()));
      previously_displayed_percentage.insert(progress.Percent());
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr); }
/* ************************************************************************* */
