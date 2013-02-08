////////////////////////////////////////////////////////////////////////////////

//   Author:    Andy Rushton
//   Copyright: (c) Southampton University 1999-2004
//              (c) Andy Rushton           2004 onwards
//   License:   BSD License, see ../docs/license.html

//   Simple wildcard matching function.

//   WARNING: wheel re-invention follows
//   Given that all shells perform wildcard matching, why don't the library writers put it in the C run-time????????

////////////////////////////////////////////////////////////////////////////////
#include "./wildcard.hpp"

namespace stlplus
{

  // function for testing whether a character matches a set
  // I can't remember the exact rules and I have no definitive references but:
  // a set contains characters, escaped characters (I think) and ranges in the form a-z
  // The character '-' can only appear at the start of the set where it is not interpreted as a range
  // This is a horrible mess - blame the Unix folks for making a hash of wildcards
  // first expand any ranges and remove escape characters to make life more palatable

  static bool match_set (const std::string& set, char match)
  {
    std::string simple_set;
    for (std::string::const_iterator i = set.begin(); i != set.end(); ++i)
    {
      switch(*i)
      {
      case '-':
      {
        if (i == set.begin())
        {
          simple_set += *i;
        }
        else if (i+1 == set.end())
        {
          return false;
        }
        else
        {
          // found a set. The first character is already in the result, so first remove it (the set might be empty)
          simple_set.erase(simple_set.end()-1);
          char last = *++i;
          for (char ch = *(i-2); ch <= last; ch++)
          {
            simple_set += ch;
          }
        }
        break;
      }
      case '\\':
        if (i+1 == set.end()) {return false;}
        simple_set += *++i;
        break;
      default:
        simple_set += *i;
        break;
      }
    }
    std::string::size_type result = simple_set.find(match);
    return result != std::string::npos;
  }

  // the recursive bit - basically whenever a * is found you recursively call this for each candidate substring match
  // until either it succeeds or you run out of string to match
  // for each * in the wildcard another level of recursion is created

  static bool match_remainder (const std::string& wild, std::string::const_iterator wildi, const std::string& match, std::string::const_iterator matchi)
  {
    //cerr << "match_remainder called at " << *matchi << " with wildcard " << *wildi << endl;
    while (wildi != wild.end() && matchi != match.end())
    {
      //cerr << "trying to match " << *matchi << " with wildcard " << *wildi << endl;
      switch(*wildi)
      {
      case '*':
      {
        ++wildi;
        ++matchi;
        for (std::string::const_iterator i = matchi; i != match.end(); ++i)
        {
          // deal with * at the end of the wildcard - there is no remainder then
          if (wildi == wild.end())
          {
            if (i == match.end()-1)
              return true;
          }
          else if (match_remainder(wild, wildi, match, i))
          {
            return true;
          }
        }
        return false;
      }
      case '[':
      {
        // scan for the end of the set using a similar method for avoiding escaped characters
        bool found = false;
        std::string::const_iterator end = wildi + 1;
        for (; !found && end != wild.end(); ++end)
        {
          switch(*end)
          {
          case ']':
          {
            // found the set, now match with its contents excluding the brackets
            if (!match_set(wild.substr(wildi - wild.begin() + 1, end - wildi - 1), *matchi))
              return false;
            found = true;
            break;
          }
          case '\\':
            if (end == wild.end()-1)
              return false;
            ++end;
            break;
          default:
            break;
          }
        }
        if (!found)
          return false;
        ++matchi;
        wildi = end;
        break;
      }
      case '?':
        ++wildi;
        ++matchi;
        break;
      case '\\':
        if (wildi == wild.end()-1)
          return false;
        ++wildi;
        if (*wildi != *matchi)
          return false;
        ++wildi;
        ++matchi;
        break;
      default:
        if (*wildi != *matchi)
          return false;
        ++wildi;
        ++matchi;
        break;
      }
    }
    bool result = wildi == wild.end() && matchi == match.end();
    return result;
  }

  // like all recursions the exported function has a simpler interface than the
  // recursive function and is just a 'seed' to the recursion itself

  bool wildcard(const std::string& wild, const std::string& match)
  {
    return match_remainder(wild, wild.begin(), match, match.begin());
  }

} // end namespace stlplus
