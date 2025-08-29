// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2023 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: vitus@google.com (Michael Vitus)

#include "ceres/concurrent_queue.h"

#include <chrono>
#include <thread>

#include "ceres/internal/config.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace ceres::internal {

// A basic test of push and pop.
TEST(ConcurrentQueue, PushPop) {
  ConcurrentQueue<int> queue;

  const int num_to_add = 10;
  for (int i = 0; i < num_to_add; ++i) {
    queue.Push(i);
  }

  for (int i = 0; i < num_to_add; ++i) {
    int value;
    ASSERT_TRUE(queue.Pop(&value));
    EXPECT_EQ(i, value);
  }
}

// Push and pop elements from the queue after StopWaiters has been called.
TEST(ConcurrentQueue, PushPopAfterStopWaiters) {
  ConcurrentQueue<int> queue;

  const int num_to_add = 10;
  int value;

  // Pop should return immediately with false with an empty queue.
  ASSERT_FALSE(queue.Pop(&value));

  for (int i = 0; i < num_to_add; ++i) {
    queue.Push(i);
  }

  // Call stop waiters to ensure we can still Push and Pop from the queue.
  queue.StopWaiters();

  for (int i = 0; i < num_to_add; ++i) {
    ASSERT_TRUE(queue.Pop(&value));
    EXPECT_EQ(i, value);
  }

  // Pop should return immediately with false with an empty queue.
  ASSERT_FALSE(queue.Pop(&value));

  // Ensure we can still push onto the queue after StopWaiters has been called.
  const int offset = 123;
  for (int i = 0; i < num_to_add; ++i) {
    queue.Push(i + offset);
  }

  for (int i = 0; i < num_to_add; ++i) {
    int value;
    ASSERT_TRUE(queue.Pop(&value));
    EXPECT_EQ(i + offset, value);
  }

  // Pop should return immediately with false with an empty queue.
  ASSERT_FALSE(queue.Pop(&value));

  // Try calling StopWaiters again to ensure nothing changes.
  queue.StopWaiters();

  queue.Push(13456);
  ASSERT_TRUE(queue.Pop(&value));
  EXPECT_EQ(13456, value);
}

// Push and pop elements after StopWaiters and EnableWaiters has been called.
TEST(ConcurrentQueue, PushPopStopAndStart) {
  ConcurrentQueue<int> queue;

  int value;

  queue.Push(13456);
  queue.Push(256);

  queue.StopWaiters();

  ASSERT_TRUE(queue.Pop(&value));
  EXPECT_EQ(13456, value);

  queue.EnableWaiters();

  // Try adding another entry after enable has been called.
  queue.Push(989);

  // Ensure we can pop both elements off.
  ASSERT_TRUE(queue.Pop(&value));
  EXPECT_EQ(256, value);

  ASSERT_TRUE(queue.Pop(&value));
  EXPECT_EQ(989, value);

  // Re-enable waiting.
  queue.EnableWaiters();

  // Pop should return immediately with false with an empty queue.
  ASSERT_FALSE(queue.Pop(&value));
}

// A basic test for Wait.
TEST(ConcurrentQueue, Wait) {
  ConcurrentQueue<int> queue;

  int value;

  queue.Push(13456);

  ASSERT_TRUE(queue.Wait(&value));
  EXPECT_EQ(13456, value);

  queue.StopWaiters();

  // Ensure waiting returns immediately after StopWaiters.
  EXPECT_FALSE(queue.Wait(&value));
  EXPECT_FALSE(queue.Wait(&value));

  EXPECT_FALSE(queue.Pop(&value));

  // Calling StopWaiters multiple times does not change anything.
  queue.StopWaiters();

  EXPECT_FALSE(queue.Wait(&value));
  EXPECT_FALSE(queue.Wait(&value));

  queue.Push(989);
  queue.Push(789);

  ASSERT_TRUE(queue.Wait(&value));
  EXPECT_EQ(989, value);

  ASSERT_TRUE(queue.Wait(&value));
  EXPECT_EQ(789, value);
}

// Ensure wait blocks until an element is pushed. Also ensure wait does not
// block after StopWaiters is called and there is no value in the queue.
// Finally, ensures EnableWaiters re-enables waiting.
TEST(ConcurrentQueue, EnsureWaitBlocks) {
  ConcurrentQueue<int> queue;

  int value = 0;
  bool valid_value = false;
  bool waiting = false;
  std::mutex mutex;

  std::thread thread([&]() {
    {
      std::lock_guard<std::mutex> lock(mutex);
      waiting = true;
    }

    int element = 87987;
    bool valid = queue.Wait(&element);

    {
      std::lock_guard<std::mutex> lock(mutex);
      waiting = false;
      value = element;
      valid_value = valid;
    }
  });

  // Give the thread time to start and wait.
  std::this_thread::sleep_for(std::chrono::milliseconds(500));

  // Ensure nothing is has been popped off the queue
  {
    std::lock_guard<std::mutex> lock(mutex);
    EXPECT_TRUE(waiting);
    ASSERT_FALSE(valid_value);
    ASSERT_EQ(0, value);
  }

  queue.Push(13456);

  // Wait for the thread to pop the value.
  thread.join();

  EXPECT_TRUE(valid_value);
  EXPECT_EQ(13456, value);
}

TEST(ConcurrentQueue, StopAndEnableWaiters) {
  ConcurrentQueue<int> queue;

  int value = 0;
  bool valid_value = false;
  bool waiting = false;
  std::mutex mutex;

  auto task = [&]() {
    {
      std::lock_guard<std::mutex> lock(mutex);
      waiting = true;
    }

    int element = 87987;
    bool valid = queue.Wait(&element);

    {
      std::lock_guard<std::mutex> lock(mutex);
      waiting = false;
      value = element;
      valid_value = valid;
    }
  };

  std::thread thread_1(task);

  // Give the thread time to start and wait.
  std::this_thread::sleep_for(std::chrono::milliseconds(500));

  // Ensure the thread is waiting.
  {
    std::lock_guard<std::mutex> lock(mutex);
    EXPECT_TRUE(waiting);
  }

  // Unblock the thread.
  queue.StopWaiters();

  thread_1.join();

  // Ensure nothing has been popped off the queue.
  EXPECT_FALSE(valid_value);
  EXPECT_EQ(87987, value);

  // Ensure another call to Wait returns immediately.
  EXPECT_FALSE(queue.Wait(&value));

  queue.EnableWaiters();

  value = 0;
  valid_value = false;
  waiting = false;

  // Start another task waiting for an element to be pushed.
  std::thread thread_2(task);

  // Give the thread time to start and wait.
  std::this_thread::sleep_for(std::chrono::milliseconds(500));

  // Ensure nothing is popped off the queue.
  {
    std::lock_guard<std::mutex> lock(mutex);
    EXPECT_TRUE(waiting);
    ASSERT_FALSE(valid_value);
    ASSERT_EQ(0, value);
  }

  queue.Push(13456);

  // Wait for the thread to pop the value.
  thread_2.join();

  EXPECT_TRUE(valid_value);
  EXPECT_EQ(13456, value);
}

}  // namespace ceres::internal
