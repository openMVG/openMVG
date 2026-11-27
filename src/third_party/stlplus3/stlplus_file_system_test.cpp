// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "filesystemSimplified/file_system.hpp"
#include "CppUnitLite/TestHarness.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>

#ifdef MSWINDOWS
#include <windows.h>
#include <winioctl.h>
#else
#include <unistd.h>
#endif

namespace
{
// Resize a real temporary file without allocating its logical size on disk.
// The Windows case explicitly requests a sparse file before extending it.
class TemporarySizedFile
{
public:
  explicit TemporarySizedFile(std::uint64_t size)
  {
#ifdef MSWINDOWS
    char directory[MAX_PATH] = {};
    char filename[MAX_PATH] = {};
    const DWORD length = GetTempPathA(MAX_PATH, directory);
    if (length == 0 || length >= MAX_PATH ||
        GetTempFileNameA(directory, "mvg", 0, filename) == 0)
      return;
    path_ = filename;
    HANDLE file = CreateFileA(path_.c_str(), GENERIC_WRITE,
                             FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
                             nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (file == INVALID_HANDLE_VALUE)
      return;
    DWORD returned = 0;
    const bool sparse = size <= 0x7fffffffULL ||
        DeviceIoControl(file, FSCTL_SET_SPARSE, nullptr, 0, nullptr, 0,
                        &returned, nullptr) != 0;
    LARGE_INTEGER end;
    end.QuadPart = static_cast<LONGLONG>(size);
    valid_ = sparse && SetFilePointerEx(file, end, nullptr, FILE_BEGIN) &&
             SetEndOfFile(file);
    CloseHandle(file);
#else
    char filename[] = "openmvg-stlplus-XXXXXX";
    const int file = mkstemp(filename);
    if (file < 0)
      return;
    path_ = filename;
    valid_ = ftruncate(file, static_cast<off_t>(size)) == 0;
    close(file);
#endif
  }

  ~TemporarySizedFile()
  {
    if (!path_.empty())
      std::remove(path_.c_str());
  }

  bool valid() const { return valid_; }
  const std::string& path() const { return path_; }

private:
  std::string path_;
  bool valid_ = false;
};
}

TEST(StlplusFileSystem, LargeFileStatus)
{
  const std::uint64_t size = (std::uint64_t(1) << 31) + 17;
  TemporarySizedFile file(size);
  CHECK(file.valid());
  if (!file.valid()) return;

  CHECK(stlplus::is_present(file.path()));
  CHECK(stlplus::is_file(file.path()));
  CHECK(!stlplus::is_folder(file.path()));
  CHECK(static_cast<std::uint64_t>(stlplus::file_size(file.path())) == size);
  CHECK(stlplus::file_created(file.path()) > 0);
  CHECK(stlplus::file_modified(file.path()) > 0);
  CHECK(stlplus::file_accessed(file.path()) > 0);
  CHECK(stlplus::file_delete(file.path()));
  CHECK(!stlplus::is_present(file.path()));
}

TEST(StlplusFileSystem, SmallFileStatus)
{
  TemporarySizedFile file(17);
  CHECK(file.valid());
  if (!file.valid()) return;

  CHECK(stlplus::is_present(file.path()));
  CHECK(stlplus::is_file(file.path()));
  CHECK(!stlplus::is_folder(file.path()));
  CHECK(stlplus::file_size(file.path()) == 17);
  CHECK(stlplus::file_created(file.path()) > 0);
  CHECK(stlplus::file_modified(file.path()) > 0);
  CHECK(stlplus::file_accessed(file.path()) > 0);
}

TEST(StlplusFileSystem, MissingFileStatus)
{
  TemporarySizedFile file(0);
  CHECK(file.valid());
  if (!file.valid()) return;
  CHECK(std::remove(file.path().c_str()) == 0);

  CHECK(!stlplus::is_present(file.path()));
  CHECK(!stlplus::is_file(file.path()));
  CHECK(!stlplus::is_folder(file.path()));
  CHECK(stlplus::file_size(file.path()) == 0);
  CHECK(stlplus::file_created(file.path()) == 0);
  CHECK(stlplus::file_modified(file.path()) == 0);
  CHECK(stlplus::file_accessed(file.path()) == 0);
}

TEST(StlplusFileSystem, DirectoryStatus)
{
  CHECK(stlplus::is_present("."));
  CHECK(stlplus::is_folder("."));
  CHECK(!stlplus::is_file("."));
  CHECK(stlplus::is_present("./"));
  CHECK(stlplus::is_folder("./"));
}

int main()
{
  TestResult result;
  TestRegistry::runAllTests(result);
  return result.getFailureCount();
}
