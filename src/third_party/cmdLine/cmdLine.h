/**
 * @file cmdLine.h
 * @brief Command line option parsing
 * @author Pascal Monasse
 * 
 * Copyright (c) 2012 Pascal Monasse
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
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
        if(std::string("-")+c==argv[0] ||
           (!longName.empty() && std::string("--")+longName==argv[0])) {
            std::stringstream str;
            if(argc<=1)
                throw std::string("Error: Option ")
                    +argv[0]+" requires argument";
            str << argv[1];
            if((str >> _field).fail() || !str.eof())
                throw std::string("Error: Unable to interpret ")
                    +argv[1]+" as argument of "+argv[0];
            used = true;
            std::rotate(argv, argv+2, argv+argc);
            argc -= 2;
            return true;
        }
        return false;
    }
    /// Copy
    Option* clone() const {
        return new OptionField<T>(c, _field, longName);
    }
private:
    T& _field; ///< Reference to variable where to store the value
};

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
            if(! found) {
                if(std::string(argv[i]).size()>1 && argv[i][0] == '-')
                    throw std::string("Unrecognized option ")+argv[i];
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
