/**
 * @file cmdLine.h
 * @brief Command line option parsing
 * @author Pascal Monasse
 * 
 * Copyright (c) 2012-2013 Pascal Monasse
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CMDLINE_H
#define CMDLINE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <cassert>

/// Base class for option/switch
class Option {
public:
    char c; ///< Option letter (eg 's' for option -s)
    bool used; ///< Does the command line use that option?
    std::string longName; /// Optional long name (eg "switch" for --switch)

    /// Constructor with short name/long name
    Option(char d, std::string name)
    : c(d), used(false), longName(name) {}
    virtual ~Option() {}
    virtual bool check(int& argc, char* argv[])=0; ///< Option found at argv[0]?
    virtual Option* clone() const=0; ///< Copy
};

/// Option on/off is called a switch
class OptionSwitch : public Option {
public:
    /// Constructor with short name/long name (optional)
    OptionSwitch(char c, std::string name="")
    : Option(c,name) {}
    /// Find switch in argv[0]
    bool check(int& argc, char* argv[]) {
        if(std::string("-")+c==argv[0] ||
           (!longName.empty() && std::string("--")+longName==argv[0])) {
            used = true;
            std::rotate(argv, argv+1, argv+argc);
            argc -= 1;
            return true;
        } else if(std::string(argv[0]).find(std::string("-")+c)==0) {
            used = true; // Handle multiple switches in single option
            std::rotate(argv[0]+1, argv[0]+2,
                        argv[0]+std::string(argv[0]).size()+1);
            return true;
        }
        return false;
    }
    /// Copy
    Option* clone() const {
        return new OptionSwitch(c, longName);
    }
};

/// Option with an argument of type T. The type must be readable by operator>>
template <class T>
class OptionField : public Option {
public:
    /// Constructor. The result with be stored in variable @field.
    OptionField(char c, T& field, std::string name="")
    : Option(c,name), _field(field) {}
    /// Find option in argv[0] and argument in argv[1]. Throw an exception
    /// (type std::string) if the argument cannot be read.
    bool check(int& argc, char* argv[]) {
        std::string param; int arg=0;
        if(std::string("-")+c==argv[0] ||
           (!longName.empty() && std::string("--")+longName==argv[0])) {
            if(argc<=1)
                throw std::string("Option ")
                    +argv[0]+" requires argument";
            param=argv[1]; arg=2;
        } else if(std::string(argv[0]).find(std::string("-")+c)==0) {
            param=argv[0]+2; arg=1;
        } else if(!longName.empty() &&
                  std::string(argv[0]).find(std::string("--")+longName+'=')==0){
            size_t size=(std::string("--")+longName+'=').size();
            param=std::string(argv[0]).substr(size); arg=1;
        }
        if(arg>0) {
            if(! read_param(param))
                throw std::string("Unable to interpret ")
                    +param+" as argument of "+argv[0];
            used = true;
            std::rotate(argv, argv+arg, argv+argc);
            argc -= arg;
            return true;
        }
        return false;
    }
    bool read_param(const std::string& param) {
        std::stringstream str(param);
        return !((str >> _field).fail() || !str.eof());
    }
    /// Copy
    Option* clone() const {
        return new OptionField<T>(c, _field, longName);
    }
private:
    T& _field; ///< Reference to variable where to store the value
};

/// Template specialization to be able to take parameter including space.
/// Generic method would do >>_field (stops at space) and test eof (false).
template <>
inline bool OptionField<std::string>::read_param(const std::string& param) {
    _field = param;
    return true;
}

/// New switch option
OptionSwitch make_switch(char c, std::string name="") {
    return OptionSwitch(c, name);
}

/// New option with argument.
template <class T>
OptionField<T> make_option(char c, T& field, std::string name="") {
    return OptionField<T>(c, field, name);
}

/// Command line parsing
class CmdLine {
    std::vector<Option*> opts;
public:
    /// Destructor
    ~CmdLine() {
        std::vector<Option*>::iterator it=opts.begin();
        for(; it != opts.end(); ++it)
            delete *it;
    }
    /// Add an option
    void add(const Option& opt) {
        opts.push_back( opt.clone() );
    }
    /// Parse of command line acting as a filter. All options are virtually
    /// removed from the command line.
    void process(int& argc, char* argv[]) throw(std::string) {
        std::vector<Option*>::iterator it=opts.begin();
        for(; it != opts.end(); ++it)
            (*it)->used = false;
        for(int i=1; i<argc;) {
            if(std::string("--")==argv[i]) { // "--" means stop option parsing
                std::rotate(argv+i, argv+i+1, argv+argc);
                -- argc;
                break;
            }
            bool found=false; // Find option
            for(it=opts.begin(); it != opts.end(); ++it) {
                int n = argc-i;
                found = (*it)->check(n, argv+i);
                if(found) {
                    argc = n+i;
                    break;
                }
            }
            if(! found) { // A single dash and a negative number are not options
                if(std::string(argv[i]).size()>1 && argv[i][0] == '-') {
                    std::istringstream str(argv[i]);
                    float v;
                    if(! (str>>v).eof())
                        throw std::string("Unrecognized option ")+argv[i];
                }
                ++i;
            }
        }
    }
    /// Was the option used in last parsing?
    bool used(char c) const {
        std::vector<Option*>::const_iterator it=opts.begin();
        for(; it != opts.end(); ++it)
            if((*it)->c == c)
                return (*it)->used;
        assert(false); // Called with non-existent option, probably a bug
        return false;
    }
};

#endif
