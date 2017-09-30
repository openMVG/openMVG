////////////////////////////////////////////////////////////////////////////////

//   Author:    Andy Rushton
//   Copyright: (c) Southampton University 1999-2004
//              (c) Andy Rushton           2004 onwards
//   License:   BSD License, see ../docs/license.html

//   This is a portable interface to the file system.

//   The idea is that you write all file system access code using these functions,
//   which are ported to all platforms that we are interested in. Therefore your
//   code is inherently portable.

//   Native Windows version: switched on by macro MSWINDOWS which is
//     enabled in portability_fixes.hpp if _WIN32 is on, which in turn
//     is defined by default in VC++/Borland/Mingw compilers

//   Unix/Gnu version: default variant, no compiler directives are
//     required but _WIN32 must be absent

//   Cygwin/Gnu version: as Unix version but with additional support
//     for Windows drive letters enabled by the CYGWIN compiler directive

////////////////////////////////////////////////////////////////////////////////
#include "./file_system.hpp"
#include "./wildcard.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <ctype.h>

#ifdef MSWINDOWS
#include <windows.h>
#include <dos.h>
#include <direct.h>
#include <fcntl.h>
#include <io.h>
#include <sys/types.h>
#include <sys/stat.h>
#else
#include <dirent.h>
#include <fcntl.h>
#include <sys/param.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

////////////////////////////////////////////////////////////////////////////////

namespace stlplus
{

////////////////////////////////////////////////////////////////////////////////
// definitions of separators

#ifdef MSWINDOWS
  static const char* separator_set = "\\/";
  static const char preferred_separator = '\\';
#else
  static const char* separator_set = "/";
  static const char preferred_separator = '/';
#endif

  static bool is_separator (char ch)
  {
    for (int i = 0; separator_set[i]; i++)
    {
      if (separator_set[i] == ch)
        return true;
    }
    return false;
  }

////////////////////////////////////////////////////////////////////////////////
// implement string comparison of paths - Unix is case-sensitive, Windoze is case-insensitive

#ifdef MSWINDOWS

  static std::string lowercase(const std::string& val)
  {
    std::string text = val;
    for (unsigned i = 0; i < text.size(); i++)
      text[i] = tolower(text[i]);
    return text;
  }

#endif

  bool path_compare(const std::string& l, const std::string& r)
  {
#ifdef MSWINDOWS
    return lowercase(l) == lowercase(r);
#else
    return l == r;
#endif
  }

////////////////////////////////////////////////////////////////////////////////
// Internal data structure used to hold the different parts of a filespec

  class file_specification
  {
  private:
    bool m_relative;                 // true = relative, false = absolute
    std::string m_drive;             // drive - drive letter (e.g. "c:") or the path for an UNC (e.g. "\\somewhere")
                                     //         empty if not known or on Unix
    std::vector<std::string> m_path; // the subdirectory path to follow from the drive
    std::string m_filename;          // the filename
  public:
    file_specification(void) : m_relative(false) {}
    ~file_specification(void) {}

    bool initialise_folder(const std::string& spec);
    bool initialise_file(const std::string& spec);
    bool simplify(void);
    bool make_absolute(const std::string& root = folder_current_full());
    bool make_absolute(const file_specification& root);
    bool make_relative(const std::string& root = folder_current_full());
    bool make_relative(const file_specification& root);
    bool relative(void) const {return m_relative;}
    bool absolute(void) const {return !relative();}
    void set_relative(void) {m_relative = true;}
    void set_absolute(void) {m_relative = false;}

    const std::string& drive(void) const {return m_drive;}
    std::string& drive(void) {return m_drive;}
    void set_drive(const std::string& drive) {m_drive = drive;}

    const std::vector<std::string>& path(void) const {return m_path;}
    std::vector<std::string>& path(void) {return m_path;}
    void set_path(const std::vector<std::string>& path) {m_path = path;}

    void add_subpath(const std::string& subpath) {m_path.push_back(subpath);}
    unsigned subpath_size(void) const {return static_cast<unsigned>(m_path.size());}
    const std::string& subpath_element(unsigned i) const {return m_path[i];}
    void subpath_erase(unsigned i) {m_path.erase(m_path.begin()+i);}

    const std::string& file(void) const {return m_filename;}
    std::string& file(void) {return m_filename;}
    void set_file(const std::string& file) {m_filename = file;}

    std::string image(void) const;
  };

  bool file_specification::initialise_folder(const std::string& folder_spec)
  {
    std::string spec = folder_spec;
    m_relative = true;
    m_drive.erase();
    m_path.clear();
    m_filename.erase();
    unsigned i = 0;
#ifdef MSWINDOWS
    // first split off the drive letter or UNC prefix on Windows
    if (spec.size() >= 2 && isalpha(spec[0]) && spec[1] == ':')
    {
      // found a drive letter
      i = 2;
      m_drive = spec.substr(0, 2);
      m_relative = false;
      // if there is a drive but no path or a relative path, get the current
      // path for this drive and prepend it to the path
      if (i == spec.size() || !is_separator(spec[i]))
      {
        // getdcwd requires the drive number (1..26) not the letter (A..Z)
        char path [MAX_PATH+1];
        int drivenum = toupper(m_drive[0]) - 'A' + 1;
        if (_getdcwd(drivenum, path, MAX_PATH+1))
        {
          // the path includes the drive so we have the drive info twice
          // need to prepend this absolute path to the spec such that any remaining relative path is still retained
          if (!is_separator(path[strlen(path)-1])) spec.insert(2, 1, preferred_separator);
          spec.insert(2, path+2);
        }
        else
        {
          // non-existent drive - fill in just the root directory
          spec.insert(2, 1, preferred_separator);
        }
      }
    }
    else if (spec.size() >= 2 && is_separator(spec[0]) && is_separator(spec[1]))
    {
      // found an UNC prefix
      i = 2;
      // find the end of the prefix by scanning for the next seperator or the end of the spec
      while (i < spec.size() && !is_separator(spec[i])) i++;
      m_drive = spec.substr(0, i);
      m_relative = false;
    }
#endif
#ifdef CYGWIN
    // first split off the drive letter or UNC prefix on Windows - the Cygwin environment supports these too
    if (spec.size() >= 2 && isalpha(spec[0]) && spec[1] == ':')
    {
      // found a drive letter
      i = 2;
      m_drive = spec.substr(0, 2);
      m_relative = false;
      // if there is a drive but no path or a relative path, get the current
      // path for this drive and prepend it to the path
      if (i == spec.size() || !is_separator(spec[i]))
      {
        // non-existent drive - fill in just the root directory
        spec.insert(2, 1, preferred_separator);
      }
    }
    else if (spec.size() >= 2 && is_separator(spec[0]) && is_separator(spec[1]))
    {
      // found an UNC prefix
      i = 2;
      // find the end of the prefix by scanning for the next seperator or the end of the spec
      while (i < spec.size() && !is_separator(spec[i])) i++;
      m_drive = spec.substr(0, i);
      m_relative = false;
    }
#endif
    // check whether the path is absolute or relative and discard the leading / if absolute
    if (i < spec.size() && is_separator(spec[i]))
    {
      m_relative = false;
      i++;
#ifdef MSWINDOWS
      // if there's no drive, fill it in on Windows since absolute paths must have a drive
      if (m_drive.empty())
      {
        m_drive += (char)(_getdrive() - 1 + 'A');
        m_drive += ':';
      }
#endif
    }
    // now extract the path elements - note that a trailing / is not significant since /a/b/c/ === /a/b/c
    // also note that the leading / has been discarded - all paths are relative
    // if absolute() is set, then paths are relative to the drive, else they are relative to the current path
    unsigned start = i;
    while(i <= spec.size())
    {
      // check for element terminated by either a separator or the end of the string
      if ((i == spec.size()) || is_separator(spec[i]))
      {
        // path element found
        std::string element = spec.substr(start, i-start);
        m_path.push_back(element);
        start = i+1;
      }
      i++;
    }
    // TODO - some error handling?
    return simplify();
  }

  bool file_specification::initialise_file(const std::string& spec)
  {
    m_filename.erase();
    // remove last element as the file and then treat the rest as a folder
    size_t i = spec.size();
    while (--i)
    {
      if (is_separator(spec[i]))
        break;
#ifdef MSWINDOWS
      // on windoze you can say a:fred.txt so the colon separates the path from the filename
      else if (i == 1 && spec[i] == ':')
        break;
#endif
    }
    bool result = initialise_folder(spec.substr(0,i+1));
    m_filename = spec.substr(i+1,spec.size()-i-1);
    // TODO - some error handling?
    return result;
  }

  bool file_specification::simplify(void)
  {
    // simplify the path by removing unnecessary empty, . and .. entries - Note that zero-length entries are treated like .
    for (unsigned i = 0; i < m_path.size(); )
    {
      if (m_path[i].empty() || m_path[i].compare(".") == 0)
      {
        // found . or null
        // these both mean do nothing - so simply delete this element
        m_path.erase(m_path.begin()+i);
      }
      else if (m_path[i].compare("..") == 0)
      {
        // found ..
        if (i == 0 && !m_relative)
        {
          // up from the root does nothing so can be deleted
          m_path.erase(m_path.begin()+i);
        }
        else if (i == 0 || m_path[i-1].compare("..") == 0)
        {
          // the first element of a relative path or the previous element is .. then keep it
          i++;
        }
        else
        {
          // otherwise delete this element and the previous one
          // TODO - this is unsafe if there is a link in the path, this is a naive algorithm
          //      - to do it properly need to check each path element for links
          m_path.erase(m_path.begin()+i);
          m_path.erase(m_path.begin()+i-1);
          i--;
        }
      }
      // keep all other elements
      else
        i++;
    }
    // TODO - error checking?
    return true;
  }

  bool file_specification::make_absolute(const std::string& root)
  {
    // simply call the other version of make_absolute
    file_specification rootspec;
    rootspec.initialise_folder(root);
    return make_absolute(rootspec);
  }

  bool file_specification::make_absolute(const file_specification& rootspec)
  {
    // test whether already an absolute path in which case there's nothing to do
    if (!absolute())
    {
      // initialise the result with the root and make the root absolute
      file_specification result = rootspec;
      result.make_absolute();
      // now append this's relative path and filename to the root's absolute path
      for (unsigned i = 0; i < subpath_size(); i++)
        result.add_subpath(subpath_element(i));
      result.set_file(file());
      // now the result is the absolute path, so transfer it to this
      *this = result;
    }
    // and simplify to get rid of any unwanted .. or . elements
    return simplify();
  }

  bool file_specification::make_relative(const std::string& root)
  {
    // simply call the other version of make_relative
    file_specification rootspec;
    rootspec.initialise_folder(root);
    return make_relative(rootspec);
  }

  bool file_specification::make_relative(const file_specification& rootspec)
  {
    // test whether already an relative path in which case there's nothing to do
    if (!relative())
    {
      // initialise the result with the root and make the root absolute
      file_specification absolute_root = rootspec;
      absolute_root.make_absolute();
      // now compare elements of the absolute root with elements of this to find the common path
      // if the drives are different, no conversion can take place and the result must be absolute, else clear the drive
      if (!path_compare(drive(), absolute_root.drive())) return true;
      set_drive("");
      // first remove leading elements that are identical to the corresponding element in root
      unsigned i = 0;
      while(subpath_size() > 0 && 
            i < absolute_root.subpath_size() && 
            path_compare(subpath_element(0), absolute_root.subpath_element(i)))
      {
        subpath_erase(0);
        i++;
      }
      // now add a .. prefix for every element in root that is different from this
      while (i < absolute_root.subpath_size())
      {
        m_path.insert(m_path.begin(), "..");
        i++;
      }
      set_relative();
    }
    // and simplify to get rid of any unwanted .. or . elements
    return simplify();
  }

  std::string file_specification::image(void) const
  {
    std::string result = m_drive;
    if (absolute())
      result += preferred_separator;
    if (!m_path.empty())
    {
      for (unsigned i = 0; i < m_path.size(); i++)
      {
        if (i != 0) result += std::string(1,preferred_separator);
        result += m_path[i];
      }
    }
    else if (relative())
      result += '.';
    // add a trailing / to the last directory element
    if (result.empty() || !is_separator(result[result.size()-1]))
      result += preferred_separator;
    if (!m_filename.empty())
      result += m_filename;
    return result;
  }

////////////////////////////////////////////////////////////////////////////////
// classifying functions

// Under both Windows and Unix, the stat function is used for classification

// Under Gnu/Linux, the following classifications are defined
// source: Gnu/Linux man page for stat(2) http://linux.die.net/man/2/stat
//   S_IFMT 	0170000	bitmask for the file type bitfields
//   S_IFSOCK 	0140000	socket (Note this overlaps with S_IFDIR)
//   S_IFLNK 	0120000	symbolic link
//   S_IFREG 	0100000	regular file
//   S_IFBLK 	0060000	block device
//   S_IFDIR 	0040000	directory
//   S_IFCHR 	0020000	character device
//   S_IFIFO 	0010000	FIFO
// There are also some Posix-standard macros:
//   S_ISREG(m)        is it a regular file? 
//   S_ISDIR(m)        directory? 
//   S_ISCHR(m)        character device? 
//   S_ISBLK(m)        block device? 
//   S_ISFIFO(m)       FIFO (named pipe)? 
//   S_ISLNK(m)        symbolic link? (Not in POSIX.1-1996.) 
//   S_ISSOCK(m)       socket? (Not in POSIX.1-1996.)
// Under Windows, the following are defined:
// source: Header file sys/stat.h distributed with Visual Studio 10
//   _S_IFMT  (S_IFMT)   0xF000 file type mask
//   _S_IFREG (S_IFREG)  0x8000 regular
//   _S_IFDIR (S_IFDIR)  0x4000 directory
//   _S_IFCHR (S_IFCHR)  0x2000 character special
//   _S_IFIFO            0x1000 pipe

#ifdef MSWINDOWS
// file type tests are not defined for some reason on Windows despite them providing the stat() function!
#define R_OK 4
#define W_OK 2
// Posix-style macros for Windows
#ifndef S_ISREG
#define S_ISREG(mode)  ((mode & _S_IFMT) == _S_IFREG)
#endif
#ifndef S_ISDIR
#define S_ISDIR(mode)  ((mode & _S_IFMT) == _S_IFDIR)
#endif
#ifndef S_ISCHR
#define S_ISCHR(mode)  ((mode & _S_IFMT) == _S_IFCHR)
#endif
#ifndef S_ISBLK
#define S_ISBLK(mode)  (false)
#endif
#ifndef S_ISFIFO
#define S_ISFIFO(mode) ((mode & _S_IFMT) == _S_IFIFO)
#endif
#ifndef S_ISLNK
#define S_ISLNK(mode)  (false)
#endif
#ifndef S_ISSOCK
#define S_ISSOCK(mode) (false)
#endif
#endif

  bool is_present (const std::string& thing)
  {
    // strip off any trailing separator because that will cause the stat function to fail
    std::string path = thing;
    if (!path.empty() && is_separator(path[path.size()-1]))
      path.erase(path.size()-1,1);
    // now test if this thing exists using the built-in stat function
    struct stat buf;
    return stat(path.c_str(), &buf) == 0;
  }

  bool is_folder (const std::string& thing)
  {
    // strip off any trailing separator because that will cause the stat function to fail
    std::string path = thing;
    if (!path.empty() && is_separator(path[path.size()-1]))
      path.erase(path.size()-1,1);
    // now test if this thing exists using the built-in stat function and if so, is it a folder
    struct stat buf;
    if (!(stat(path.c_str(), &buf) == 0))
      return false;
    // If the object is present, see if it is a directory
    // this is the Posix-approved way of testing
    return S_ISDIR(buf.st_mode);
  }

  bool is_file (const std::string& thing)
  {
    // strip off any trailing separator because that will cause the stat function to fail
    std::string path = thing;
    if (!path.empty() && is_separator(path[path.size()-1]))
      path.erase(path.size()-1,1);
    // now test if this thing exists using the built-in stat function and if so, is it a file
    struct stat buf;
    if (!(stat(path.c_str(), &buf) == 0))
      return false;
    // If the object is present, see if it is a file or file-like object
    // Note that devices are neither folders nor files
    // this is the Posix-approved way of testing
    return S_ISREG(buf.st_mode) || S_ISLNK(buf.st_mode) || S_ISSOCK(buf.st_mode) || S_ISFIFO(buf.st_mode);
  }

////////////////////////////////////////////////////////////////////////////////
// file functions

  bool file_exists (const std::string& filespec)
  {
    return is_file(filespec);
  }

  bool file_readable (const std::string& filespec)
  {
    // a file is readable if it exists and can be read
    if (!file_exists(filespec)) return false;
    return access(filespec.c_str(),R_OK)==0;
  }

  bool file_writable (const std::string& filespec)
  {
    // a file is writable if it exists as a file and is writable or if
    // it doesn't exist but could be created and would be writable
    if (is_present(filespec))
    {
      if (!is_file(filespec)) return false;
      return access(filespec.c_str(),W_OK)==0;
    }
    std::string dir = folder_part(filespec);
    if (dir.empty()) dir = ".";
    return folder_writable(dir);
  }

  size_t file_size (const std::string& filespec)
  {
    struct stat buf;
    if (!(stat(filespec.c_str(), &buf) == 0)) return 0;
    return buf.st_size;
  }

  bool file_delete (const std::string& filespec)
  {
    if (!is_file(filespec)) return false;
    return remove(filespec.c_str())==0;
  }

  bool file_rename (const std::string& old_filespec, const std::string& new_filespec)
  {
    if (!is_file(old_filespec)) return false;
    return rename(old_filespec.c_str(), new_filespec.c_str())==0;
  }

  bool file_copy (const std::string& old_filespec, const std::string& new_filespec)
  {
    if (!is_file(old_filespec)) return false;
    // do an exact copy - to do this, use binary mode
    bool result = true;
    FILE* old_file = fopen(old_filespec.c_str(),"rb");
    FILE* new_file = fopen(new_filespec.c_str(),"wb");
    if (!old_file)
      result = false;
    else if (!new_file)
      result = false;
    else
    {
      for (int byte = getc(old_file); byte != EOF; byte = getc(old_file))
        putc(byte,new_file);
    }
    if (old_file) fclose(old_file);
    if (new_file) fclose(new_file);
    return result;
  }

  bool file_move (const std::string& old_filespec, const std::string& new_filespec)
  {
    // try to move the file by renaming - if that fails then do a copy and delete the original
    if (file_rename(old_filespec, new_filespec))
      return true;
    if (!file_copy(old_filespec, new_filespec))
      return false;
    // I'm not sure what to do if the delete fails - is that an error?
    // I've made it an error and then delete the copy so that the original state is recovered
    if (file_delete(old_filespec))
      return true;
    file_delete(new_filespec);
    return false;
  }

  time_t file_created (const std::string& filespec)
  {
    struct stat buf;
    if (!(stat(filespec.c_str(), &buf) == 0)) return 0;
    return buf.st_ctime;
  }

  time_t file_modified (const std::string& filespec)
  {
    struct stat buf;
    if (!(stat(filespec.c_str(), &buf) == 0)) return 0;
    return buf.st_mtime;
  }

  time_t file_accessed (const std::string& filespec)
  {
    struct stat buf;
    if (!(stat(filespec.c_str(), &buf) == 0)) return 0;
    return buf.st_atime;
  }

  std::string create_filespec (const std::string& directory, const std::string& filename)
  {
    std::string result = directory;
    // if directory is empty then no directory part will be added
    // add trailing slash if the directory was specified and does not have a trailing slash
    if (!result.empty() && !is_separator(result[result.size()-1]))
      result += preferred_separator;
    // if filename is null or empty, nothing will be added so the path is then a directory path
    result += filename;
    return result;
  }

  std::string create_filespec (const std::string& directory, const std::string& basename, const std::string& extension)
  {
    return create_filespec(directory, create_filename(basename, extension));
  }

  std::string create_filename(const std::string& basename, const std::string& extension)
  {
    std::string name = basename;
    // extension is optional - so the dot is also optional
    if (!extension.empty())
    {
      if (extension[0] != '.') name += '.';
      name += extension;
    }
    return name;
  }

////////////////////////////////////////////////////////////////////////////////
// folder functions

  bool folder_create (const std::string& directory)
  {
#ifdef MSWINDOWS
    return mkdir(directory.c_str()) == 0;
#else
    return mkdir(directory.c_str(), 0777) == 0;
#endif
  }

  bool folder_exists (const std::string& directory)
  {
    return is_folder(directory);
  }

  bool folder_readable (const std::string& directory)
  {
    // a folder is readable if it exists and has read access
    std::string dir = directory;
    if (dir.empty()) dir = ".";
    if (!folder_exists(dir)) return false;
    return access(dir.c_str(),R_OK)==0;
  }

  bool folder_writable (const std::string& directory)
  {
    // a folder is writable if it exists and has write access
    std::string dir = directory;
    if (dir.empty()) dir = ".";
    if (!folder_exists(dir)) return false;
    return access(dir.c_str(),W_OK)==0;
  }

  bool folder_delete (const std::string& directory, bool recurse)
  {
    std::string dir = directory;
    if (dir.empty()) dir = ".";
    if (!folder_exists(dir)) return false;
    bool result = true;
    // depth-first traversal ensures that directory contents are deleted before trying to delete the directory itself
    if (recurse)
    {
      std::vector<std::string> subdirectories = folder_subdirectories(dir);
      for (std::vector<std::string>::size_type d = 0; d < subdirectories.size(); ++d)
        if (!folder_delete(folder_down(dir,subdirectories[d]),true)) 
          result = false;
      std::vector<std::string> files = folder_files(dir);
      for (std::vector<std::string>::size_type f = 0; f < files.size(); ++f)
        if (!file_delete(create_filespec(dir, files[f]))) 
          result = false;
    }
    if (rmdir(dir.c_str())!=0) result = false;
    return result;
  }

  bool folder_rename (const std::string& old_directory, const std::string& new_directory)
  {
    if (!folder_exists(old_directory)) return false;
    return rename(old_directory.c_str(), new_directory.c_str())==0;
  }

  bool folder_empty(const std::string& directory)
  {
    std::string dir = directory.empty() ? std::string(".") : directory;
    bool result = true;
#ifdef MSWINDOWS
    std::string wildcard = create_filespec(dir, "*.*");
    intptr_t handle = -1;
    _finddata_t fileinfo;
    for (bool OK = (handle = _findfirst((char*)wildcard.c_str(), &fileinfo)) != -1; OK; OK = (_findnext(handle, &fileinfo)==0))
    {
      std::string strentry = fileinfo.name;
      if (strentry.compare(".")!=0 && strentry.compare("..")!=0)
      {
        result = false;
        break;
      }
    }
    _findclose(handle);
#else
    DIR* d = opendir(dir.c_str());
    if (d)
    {
      for (dirent* entry = readdir(d); entry; entry = readdir(d))
      {
        std::string strentry = entry->d_name;
        if (strentry.compare(".")!=0 && strentry.compare("..")!=0)
        {
          result = false;
          break;
        }
      }
      closedir(d);
    }
#endif
    return result;
  }

  bool folder_set_current(const std::string& folder)
  {
    if (!folder_exists(folder))
      return false;
#ifdef MSWINDOWS
    // Windose implementation - this returns non-zero for success
    return (SetCurrentDirectoryA(folder.c_str()) != 0);
#else
    // Unix implementation - this returns zero for success
    return (chdir(folder.c_str()) == 0);
#endif
  }

  std::string folder_current (void)
  {
    return ".";
  }

  std::string folder_current_full(void)
  {
    // It's not clear from the documentation whether the buffer for a path should be one byte longer
    // than the maximum path length to allow for the null termination, so I have made it so anyway
#ifdef MSWINDOWS
    char abspath [MAX_PATH+1];
    return std::string(_fullpath(abspath, ".", MAX_PATH+1));
#else
    char pathname [MAXPATHLEN+1];
    char* result = getcwd(pathname,MAXPATHLEN+1);
    if (!result)
    {
      // should really report the error from errno
      return std::string();
    }
    return std::string(result);
#endif
  }

  std::string folder_down (const std::string& directory, const std::string& subdirectory)
  {
    file_specification spec;
    spec.initialise_folder(directory);
    spec.add_subpath(subdirectory);
    return spec.image();
  }

  std::string folder_up (const std::string& directory, unsigned levels)
  {
    file_specification spec;
    spec.initialise_folder(directory);
    // add .. elements for each level to go up, then use simplify
    // TODO - doesn't work if any of the elements are links
    for (unsigned i = 0; i < levels; i++)
      spec.add_subpath("..");
    spec.simplify();
    return spec.image();
  }

  std::vector<std::string> folder_subdirectories (const std::string& directory)
  {
    return folder_wildcard(directory, "*", true, false);
  }

  std::vector<std::string> folder_files (const std::string& directory)
  {
    return folder_wildcard(directory, "*", false, true);
  }

  std::vector<std::string> folder_all(const std::string& directory)
  {
    return folder_wildcard(directory, "*", true, true);
  }

  std::vector<std::string> folder_wildcard (const std::string& directory, const std::string& wild, bool subdirs, bool files)
  {
    std::string dir = directory.empty() ? std::string(".") : directory;
    std::vector<std::string> results;
#ifdef MSWINDOWS
    std::string wildcard = create_filespec(dir, wild);
    intptr_t handle = -1;
    _finddata_t fileinfo;
    for (bool OK = (handle = _findfirst((char*)wildcard.c_str(), &fileinfo)) != -1; OK; OK = (_findnext(handle, &fileinfo)==0))
    {
      std::string strentry = fileinfo.name;
      if (strentry.compare(".")!=0 && strentry.compare("..")!=0)
        if ((subdirs && (fileinfo.attrib & _A_SUBDIR)) || (files && !(fileinfo.attrib & _A_SUBDIR)))
          results.push_back(strentry);
    }
    _findclose(handle);
#else
    DIR* d = opendir(dir.c_str());
    if (d)
    {
      for (dirent* entry = readdir(d); entry; entry = readdir(d))
      {
        std::string strentry = entry->d_name;
        if (strentry.compare(".")!=0 && strentry.compare("..")!=0)
        {
          std::string subpath = create_filespec(dir, strentry);
          if (((subdirs && is_folder(subpath)) || (files && is_file(subpath))) && (wildcard(wild, strentry)))
            results.push_back(strentry);
        }
      }
      closedir(d);
    }
#endif
    return results;
  }

  std::string folder_home (void)
  {
    if (getenv("HOME"))
      return std::string(getenv("HOME"));
#ifdef MSWINDOWS
    if (getenv("HOMEDRIVE") || getenv("HOMEPATH"))
      return std::string(getenv("HOMEDRIVE")) + std::string(getenv("HOMEPATH"));
    return "C:\\";
#else
    if (getenv("USER"))
      return folder_down("/home", std::string(getenv("USER")));
    if (getenv("USERNAME"))
      return folder_down("/home", std::string(getenv("USERNAME")));
    return "";
#endif
  }

////////////////////////////////////////////////////////////////////////////////
// path functions convert between full and relative paths

  bool is_full_path(const std::string& path)
  {
    file_specification spec;
    spec.initialise_folder(path.empty() ? std::string(".") : path);
    return spec.absolute();
  }

  bool is_relative_path(const std::string& path)
  {
    file_specification spec;
    spec.initialise_folder(path.empty() ? std::string(".") : path);
    return spec.relative();
  }

  static std::string full_path(const std::string& root, const std::string& path)
  {
    // convert path to a full path using root as the start point for relative paths
    // decompose the path and test whether it is already an absolute path, in which case just return it
    file_specification spec;
    spec.initialise_folder(path.empty() ? std::string(".") : path);
    if (spec.relative())
    {
      // okay, so the path is relative after all, so we need to combine it with the root path
      // decompose the root path and check whether it is relative
      file_specification rootspec;
      rootspec.initialise_folder(root.empty() ? std::string(".") : root);
      if (rootspec.relative())
        rootspec.make_absolute();
      // Now do the conversion of the path relative to the root
      spec.make_absolute(rootspec);
    }
    // clean up the path
    spec.simplify();
    return spec.image();
  }

  static std::string relative_path(const std::string& root, const std::string& path)
  {
    // convert path to a relative path, using the root path as its starting point
    // first convert both paths to full paths relative to CWD
    file_specification rootspec;
    rootspec.initialise_folder(root.empty() ? std::string(".") : root);
    if (rootspec.relative())
      rootspec.make_absolute();
    file_specification spec;
    spec.initialise_folder(path.empty() ? std::string(".") : path);
    if (spec.relative())
      spec.make_absolute();
    // now make path spec relative to the root spec
    spec.make_relative(rootspec);
    // clean up the path
    spec.simplify();
    return spec.image();
  }

  std::string folder_to_path (const std::string& path, const std::string& directory)
  {
    return full_path(path, directory);
  }

  std::string filespec_to_path (const std::string& path, const std::string& spec)
  {
    return create_filespec(folder_to_path(path, folder_part(spec)),filename_part(spec));
  }

  std::string folder_to_path(const std::string& folder)
  {
    return folder_to_path(folder_current(), folder);
  }

  std::string filespec_to_path(const std::string& filespec)
  {
    return filespec_to_path(folder_current(), filespec);
  }

  std::string folder_to_relative_path(const std::string& root, const std::string& folder)
  {
    return relative_path(root, folder);
  }

  std::string filespec_to_relative_path(const std::string& root, const std::string& spec)
  {
    return create_filespec(folder_to_relative_path(root, folder_part(spec)),filename_part(spec));
  }

  std::string folder_to_relative_path(const std::string& folder)
  {
    return folder_to_relative_path(folder_current(), folder);
  }

  std::string filespec_to_relative_path(const std::string& filespec)
  {
    return filespec_to_relative_path(folder_current(), filespec);
  }

  std::string folder_append_separator(const std::string& folder)
  {
    std::string result = folder;
    if (result.empty() || !is_separator(result[result.size()-1]))
      result += preferred_separator;
    return result;
  }

  std::string folder_remove_end_separator(const std::string& folder)
  {
    std::string result = folder;
    // remove the end separator if it is not the only element in the path
#ifdef MSWINDOWS
    // allow for the drive letter
    if (result.size() >= 2 && isalpha(result[0]) && result[1] == ':')
      if ((result.size() > 3) && is_separator(result[result.size()-1]))
        result.erase(result.size()-1,1);
#endif
    if ((result.size() > 1) && is_separator(result[result.size()-1]))
      result.erase(result.size()-1,1);
    return result;
  }

////////////////////////////////////////////////////////////////////////////////

  std::string basename_part (const std::string& spec)
  {
    std::string fname = filename_part(spec);
    // scan back through filename until a '.' is found and remove suffix
    // the whole filename is the basename if there is no '.'
    std::string::size_type i = fname.find_last_of('.');
    // observe Unix convention that a dot at the start of a filename is part of the basename, not the extension
    if (i != 0 && i != std::string::npos)
      fname.erase(i, fname.size()-i);
    return fname;
  }

  std::string filename_part (const std::string& spec)
  {
    // scan back through filename until a preferred_separator is found and remove prefix;
    // if there is no preferred_separator then remove nothing, i.e. the whole filespec is filename
    size_t i = spec.size();
    while (i--)
    {
      if (is_separator(spec[i]))
        return spec.substr(i+1,spec.size()-i-1);
    }
    return spec;
  }

  std::string extension_part (const std::string& spec)
  {
    std::string fname = filename_part(spec);
    // scan back through filename until a '.' is found and remove prefix;
    std::string::size_type i = fname.find_last_of('.');
    // observe Unix convention that a dot at the start of a filename is part of the name, not the extension;
    if (i != 0 && i != std::string::npos)
      fname.erase(0, i+1);
    else
      fname.erase();
    return fname;
  }

  std::string folder_part (const std::string& spec)
  {
    // scan back through filename until a separator is found and remove prefix
    // if there is no separator, remove the whole
    size_t i = spec.size();
    while (i--)
    {
      if (is_separator(spec[i]))
        return spec.substr(0,i);
    }
    return std::string();
  }

  std::vector<std::string> filespec_elements(const std::string& filespec)
  {
    file_specification spec;
    spec.initialise_file(filespec);
    std::vector<std::string> result = spec.path();
    if (!spec.drive().empty()) result.insert(result.begin(),spec.drive());
    if (!spec.file().empty()) result.push_back(spec.file());
    return result;
  }

  std::vector<std::string> folder_elements(const std::string& folder)
  {
    file_specification spec;
    spec.initialise_folder(folder);
    std::vector<std::string> result = spec.path();
    if (!spec.drive().empty()) result.insert(result.begin(),spec.drive());
    return result;
  }

////////////////////////////////////////////////////////////////////////////////
// mimic the command lookup used by the shell

// Windows looks at the following locations:
// 1) application root
// 2) current directory
// 3) 32-bit system directory
// 4) 16-bit system directory
// 5) windows system directory
// 6) %path%
// currently only (2) and (6) has been implemented although many system folders are on the path anyway
// also implement the implied .exe extension on commands with no path (see CreateProcess documentation)
// TODO - PATHEXT handling to find non-exe executables

  std::string path_lookup (const std::string& command)
  {
    std::string path = std::string(".") + PATH_SPLITTER + getenv("PATH");
    return lookup(command, path);
  }

  std::string lookup (const std::string& command, const std::string& path, const std::string& splitter)
  {
    // first check whether the command is already a path and check whether it exists
    if (!folder_part(command).empty())
    {
      if (file_exists(command))
        return command;
    }
    else
    {
      // command is just a name - so do path lookup
      // split the path into its elements
      std::vector<std::string> paths;
      if (!path.empty())
      {
        for(std::string::size_type offset = 0;;)
        {
          std::string::size_type found = path.find(splitter, offset);
          if (found != std::string::npos)
          {
            paths.push_back(path.substr(offset, found-offset));
            offset = found + splitter.size();
          }
          else
          {
            paths.push_back(path.substr(offset, path.size()-offset));
            break;
          }
        }
      }
      // now lookup each path to see if it its the matching one
      for (unsigned i = 0; i < paths.size(); i++)
      {
        std::string spec = create_filespec(paths[i], command);
        if (file_exists(spec))
        {
          return spec;
        }
      }
    }
#ifdef MSWINDOWS
    // if there is no extension, try recursing on each possible extension
    // TODO iterate through PATHEXT
    if (extension_part(command).empty())
      return lookup(create_filespec(folder_part(command), basename_part(command), "exe"), path, splitter);
#endif
    // if path lookup failed, return empty string to indicate error
    return std::string();
  }

////////////////////////////////////////////////////////////////////////////////

  std::string install_path(const std::string& argv0)
  {
    std::string bin_directory = folder_part(argv0);
    if (bin_directory.empty())
    {
      // do path lookup to find the executable path
      bin_directory = folder_part(path_lookup(argv0));
    }
    return bin_directory;
  }

////////////////////////////////////////////////////////////////////////////////

} // end namespace stlplus
