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

#include "ceres/thread_pool.h"

#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>

#include "ceres/internal/config.h"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace ceres::internal {

// Adds a number of tasks to the thread pool and ensures they all run.
TEST(ThreadPool, AddTask) {
  int value = 0;
  const int num_tasks = 100;
  {
    ThreadPool thread_pool(2);

    std::condition_variable condition;
    std::mutex mutex;

    for (int i = 0; i < num_tasks; ++i) {
      thread_pool.AddTask([&]() {
        std::lock_guard<std::mutex> lock(mutex);
        ++value;
        condition.notify_all();
      });
    }

    std::unique_lock<std::mutex> lock(mutex);
    condition.wait(lock, [&]() { return value == num_tasks; });
  }

  EXPECT_EQ(num_tasks, value);
}

// Adds a number of tasks to the queue and resizes the thread pool while the
// threads are executing their work.
TEST(ThreadPool, ResizingDuringExecution) {
  int value = 0;

  const int num_tasks = 100;

  // Run this test in a scope to delete the thread pool and all of the threads
  // are stopped.
  {
    ThreadPool thread_pool(/*num_threads=*/2);

    std::condition_variable condition;
    std::mutex mutex;

    // Acquire a lock on the mutex to prevent the threads from finishing their
    // execution so we can test resizing the thread pool while the workers are
    // executing a task.
    std::unique_lock<std::mutex> lock(mutex);

    // The same task for all of the workers to execute.
    auto task = [&]() {
      // This will block until the mutex is released inside the condition
      // variable.
      std::lock_guard<std::mutex> lock(mutex);
      ++value;
      condition.notify_all();
    };

    // Add the initial set of tasks to run.
    for (int i = 0; i < num_tasks / 2; ++i) {
      thread_pool.AddTask(task);
    }

    // Resize the thread pool while tasks are executing.
    thread_pool.Resize(/*num_threads=*/3);

    // Add more tasks to the thread pool to guarantee these are also completed.
    for (int i = 0; i < num_tasks / 2; ++i) {
      thread_pool.AddTask(task);
    }

    // Unlock the mutex to unblock all of the threads and wait until all of the
    // tasks are completed.
    condition.wait(lock, [&]() { return value == num_tasks; });
  }

  EXPECT_EQ(num_tasks, value);
}

// Tests the destructor will wait until all running tasks are finished before
// destructing the thread pool.
TEST(ThreadPool, Destructor) {
  // Ensure the hardware supports more than 1 thread to ensure the test will
  // pass.
  const int num_hardware_threads = std::thread::hardware_concurrency();
  if (num_hardware_threads <= 1) {
    LOG(ERROR)
        << "Test not supported, the hardware does not support threading.";
    return;
  }

  std::condition_variable condition;
  std::mutex mutex;
  // Lock the mutex to ensure the tasks are blocked.
  std::unique_lock<std::mutex> master_lock(mutex);
  int value = 0;

  // Create a thread that will instantiate and delete the thread pool.  This is
  // required because we need to block on the thread pool being deleted and
  // signal the tasks to finish.
  std::thread thread([&]() {
    ThreadPool thread_pool(/*num_threads=*/2);

    for (int i = 0; i < 100; ++i) {
      thread_pool.AddTask([&]() {
        // This will block until the mutex is released inside the condition
        // variable.
        std::lock_guard<std::mutex> lock(mutex);
        ++value;
        condition.notify_all();
      });
    }
    // The thread pool should be deleted.
  });

  // Give the thread pool time to start, add all the tasks, and then delete
  // itself.
  std::this_thread::sleep_for(std::chrono::milliseconds(500));

  // Unlock the tasks.
  master_lock.unlock();

  // Wait for the thread to complete.
  thread.join();

  EXPECT_EQ(100, value);
}

TEST(ThreadPool, Resize) {
  // Ensure the hardware supports more than 1 thread to ensure the test will
  // pass.
  const int num_hardware_threads = std::thread::hardware_concurrency();
  if (num_hardware_threads <= 1) {
    LOG(ERROR)
        << "Test not supported, the hardware does not support threading.";
    return;
  }

  ThreadPool thread_pool(1);

  EXPECT_EQ(1, thread_pool.Size());

  thread_pool.Resize(2);

  EXPECT_EQ(2, thread_pool.Size());

  // Try reducing the thread pool size and verify it stays the same size.
  thread_pool.Resize(1);
  EXPECT_EQ(2, thread_pool.Size());
}

}  // namespace ceres::internal
