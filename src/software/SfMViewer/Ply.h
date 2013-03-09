#ifndef PLY_H
#define PLY_H

#include<cassert>
#include<climits> /* CHAR_BIT */
#include<cstddef>
#include<iostream>
#include<limits>
#include<fstream>
#include<sstream>
#include<string>
#include<utility>
#include<vector>

#if CHAR_BIT != 8
# error 8-bit bytes assumed
#endif /* CHAR_BIT */

/*!
 *  A class to read and write files in PLY (PoLYgon ?) format.
 *
 *  Sample code to read a PLY file:
 *
 *  \code
 *    Ply ply;
 *
 *    // Read PLY header
 *    if (!ply.open(path)) {
 *        // ...
 *    }
 *
 *    // Iterate over elements
 *    for (Ply::ElementsIterator it = ply.elements_begin();
 *         it != ply.elements_end(); ++it) {
 *        const Ply::Element& element = *it;
 *
 *        if (element.name() == "some_element_to_skip") {
 *            if (!ply.skip(element))
 *                // ...
 *            else
 *                continue;
 *        } else {
 *            const size_type& num_elements = element.count();
 *
 *            for (size_type i = 0; i != num_elements; ++i) {
 *                for (Ply::PropertiesIterator it2 =
 *                         element.properties_begin();
 *                    it2 != element.properties_end(); ++it2) {
 *                    const Ply::Property& property = *it2;
 *
 *                    if (it2->name() == "some_property_to_read") {
 *                        double d;
 *                        ply.read(*it2, d);
 *                        // ...
 *                    } else {
 *                        if (!ply.skip(*it2))
 *                             // ...
 *                         else
 *                             continue;
 *                    }
 *                }
 *            }
 *        }
 *    }
 *
 *    ply.close();
 *  \endcode
 *
 *  Sample code to write a PLY file:
 *
 *  \code
 *      Ply ply;
 *
 *      Ply::ElementList elements;
 *
 *      // ...
 *
 *      // Write PLY header
 *      if (!ply.open(path, elements, Ply::PLY_FORMAT_ASC)) {
 *          // ...
 *      }
 *
 *      // Iterate over elements
 *      for (Ply::ElementsIterator it = ply.elements_begin();
 *           it != ply.elements_end(); ++it) {
 *          const Ply::Element& element = *it;
 *
 *          const Ply::string_type name = element.name();
 *          const size_type& num_elements = element.count();
 *
 *          for (size_type i = 0; i != num_elements; ++i) {
 *              for (Ply::PropertiesIterator it2 =
 *                       element.properties_begin();
 *                  it2 = element.properties_end(); ++it2) {
 *                  const Ply::Property& property = *it2;
 *
 *                  if (it2->name() == "some_property_to_read") {
 *                      double d;
 *                      ply.write(*it2, d);
 *                      // ...
 *                  } else
 *                      // ...
 *              }
 *          }
 *      }
 *
 *      ply.close();
 *  \endcode
 */
class Ply {
    public:
        typedef std::size_t size_type;
        typedef std::string string_type;

        /* Comments ***********************************************/

        typedef std::vector<string_type> CommentList;
        typedef CommentList::const_iterator CommentsIterator;

        /* Object informations ************************************/

        typedef std::vector<string_type> ObjInfoList;
        typedef ObjInfoList::const_iterator ObjInfosIterator;

        /* Format *************************************************/

        enum Format {
            PLY_FORMAT_ASC     = 0,
            PLY_FORMAT_BIN_BE  = 1,
            PLY_FORMAT_BIN_LE  = 2,
            PLY_FORMAT_UNKNOWN = 3,
        };

        /* Type ***************************************************/

        enum Type {
            PLY_TYPE_CHAR    = 0,
            PLY_TYPE_INT8    = PLY_TYPE_CHAR,
            PLY_TYPE_UCHAR   = 1,
            PLY_TYPE_UINT8   = PLY_TYPE_UCHAR,
            PLY_TYPE_SHORT   = 2,
            PLY_TYPE_INT16   = PLY_TYPE_SHORT,
            PLY_TYPE_USHORT  = 3,
            PLY_TYPE_UINT16  = PLY_TYPE_USHORT,
            PLY_TYPE_INT     = 4,
            PLY_TYPE_INT32   = PLY_TYPE_INT,
            PLY_TYPE_UINT    = 5,
            PLY_TYPE_UINT32  = PLY_TYPE_UINT,
            PLY_TYPE_FLOAT   = 6,
            PLY_TYPE_FLOAT32 = PLY_TYPE_FLOAT,
            PLY_TYPE_DOUBLE  = 7,
            PLY_TYPE_FLOAT64 = PLY_TYPE_DOUBLE,
            PLY_TYPE_LIST    = 8,
            PLY_TYPE_UNKNOWN = 9,
        };

        typedef char int8_type;
        typedef unsigned char uint8_type;
        typedef short int16_type;
        typedef unsigned short uint16_type;
        typedef int int32_type;
        typedef unsigned int uint32_type;
        typedef float float32_type;
        typedef double float64_type;

        /* Property ***********************************************/

        class Property {
            public:
                /*! Constructs a scalar property. */
                Property(const string_type& name, Type type) :
                    _name(name),
                    _type(type),
                    _count_type(PLY_TYPE_UNKNOWN),
                    _value_type(PLY_TYPE_UNKNOWN)
                {
                    assert(!name.empty());
                    assert(type != PLY_TYPE_LIST);
                    assert(type != PLY_TYPE_UNKNOWN);
                }

                /*! Constructs a list property. */
                Property(const string_type& name, Type count_type,
                                                  Type value_type) :
                    _name(name),
                    _type(PLY_TYPE_LIST),
                    _count_type(count_type),
                    _value_type(value_type)
                {
                    assert(!name.empty());
                    assert(count_type != PLY_TYPE_UNKNOWN);
                    assert(count_type != PLY_TYPE_LIST);
                    assert(value_type != PLY_TYPE_UNKNOWN);
                    assert(value_type != PLY_TYPE_LIST);
                }

                /*! Returns the name of the property. */
                const string_type& name() const { return _name; }

                /*! Returns the type of the property. */
                Type type() const { return _type; }

                /*! Returns the count type of the property (For lists only)). */
                Type count_type() const
                {
                    assert(type() == PLY_TYPE_LIST);
                    return _count_type;
                }

                /*! Returns the value type of the property (For lists only)). */
                Type value_type() const
                {
                    assert(type() == PLY_TYPE_LIST);
                    return _value_type;
                }

            protected:
                /*const*/ string_type _name;
                /*const*/ Type _type;
                /*const*/ Type _count_type; /* For lists only */
                /*const*/ Type _value_type; /* For lists only */
        };

        /* Element ************************************************/

        class Element {
            public:
                typedef std::vector<Property> PropertyList;
                typedef PropertyList::const_iterator PropertiesIterator;

            public:
                /*! Constructs an element. */
                Element(const string_type& name, size_type count) :
                    _name(name),
                    _count(count) { }

                /*! Returns the name of the element. */
                const string_type& name() const { return _name; }

                /*! Returns the count of the element. */
                const size_type& count() const { return _count; }

                /*!
                 *  Returns the size required for this element
                 *       or 0 if the size is variable.
                 */
                size_type size() const
                {
                    size_type s = 0;

                    for (PropertiesIterator it = properties_begin();
                         it != properties_end(); ++it) {
                        const size_type t = type_size(it->type());

                        if (t == 0) {
                            s = 0;
                            break;
                        }

                        s += t;
                    }

                    return s;
                }

                /* Properties *************************/

                const PropertyList& properties() const { return _properties; }

                /*! Returns the number of properties of this element. */
                size_type num_properties() const { return _properties.size(); }

                PropertiesIterator properties_begin() const { return _properties.begin(); }

                PropertiesIterator properties_end() const { return _properties.end(); }

                /*! Adds a property to the element. */
                bool add_property(const Property& property)
                {
                    assert(property.type() != PLY_TYPE_UNKNOWN);
                    assert(!property.name().empty());

                    if (has_property(property.name())) {
                        std::cerr << "Ply: element \"" << this->name()
                                  << "\" already a property \""
                                  << property.name() << "\"" << std::endl;
                        return false;
                    }

                    _properties.push_back(property);

                    return true;
                }

                const Property& property(size_type i) const { return _properties[i]; }

                PropertiesIterator find_property(const string_type& name) const
                {
                    PropertiesIterator it;

                    for (it = properties_begin(); it != properties_end(); ++it)
                        if (it->name() == name)
                            return it;

                    return it; /* = properties_end() */
                }

                bool has_property(const string_type& name) const { return (find_property(name) != properties_end()); }

            protected:
                /* const */ string_type _name;
                /* const */ size_type _count;
                PropertyList _properties;
        };

        typedef Element::PropertiesIterator PropertiesIterator;
        typedef std::vector<Element> ElementList;
        typedef ElementList::const_iterator ElementsIterator;

    public:
        /*! Constructs a new Ply object. */
        Ply() :
            _format(PLY_FORMAT_UNKNOWN),
            _version(0.0)
        {
            assert(type_size(PLY_TYPE_INT8) == 1);
            assert(type_size(PLY_TYPE_UINT8) == 1);
            assert(type_size(PLY_TYPE_INT16) == 2);
            assert(type_size(PLY_TYPE_UINT16) == 2);
            assert(type_size(PLY_TYPE_INT32) == 4);
            assert(type_size(PLY_TYPE_UINT32) == 4);
            assert(type_size(PLY_TYPE_FLOAT32) == 4);
            assert(type_size(PLY_TYPE_FLOAT64) == 8);

            _find_native_format();
        }

        /* Open/Close *********************************************/

        /*! Opens the given file and write header information. */
        bool open(string_type path, const ElementList& elements,
                                    Format format = PLY_FORMAT_UNKNOWN)
        {
            ObjInfoList obj_infos;
            CommentList comments;

            return open(path, elements, obj_infos, comments, format);
        }

        /*! Opens the given file and write header information. */
        bool open(const string_type& path, const ElementList& elements,
                                           const ObjInfoList& obj_infos,
                                           const CommentList& comments,
                                           Format format = PLY_FORMAT_UNKNOWN);

        /*! Opens the given file and read header information. */
        bool open(const string_type& path);

        /*! Sets the floating-point output precision. */
        void precision(size_type num_digits10) { _os.precision(num_digits10); }

        /*! Sets the floating-point output format to scientific. */
        void scientific() { _os.setf(std::ios::scientific,
                                     std::ios::floatfield); }

        /*! Sets the floating-point output format to fixed. */
        void fixed() { _os.setf(std::ios::scientific,
                                std::ios::floatfield); }

        void close();

        /* Skip ***************************************************/

        /*! Skips the given element. */
        bool skip(const Element& element)
        {
            if (!_element_begin(element))
                return false;

            if (format() == PLY_FORMAT_ASC) {
                if (!_skip_ASC(element))
                    return false;
            } else {
                if (!_skip_BIN(element))
                    return false;
            }

            return _element_end();
        }

        /*! Skips the given property. */
        bool skip(const Property& property)
        {
            if (!_property_begin(property))
                return false;

            switch (format()) {
                case PLY_FORMAT_ASC:
                    if (!_skip_ASC(property))
                        return false;
                    break;
                case PLY_FORMAT_BIN_BE:
                case PLY_FORMAT_BIN_LE:
                    if (!_skip_BIN(property))
                        return false;
                    break;
                default: /* PLY_FORMAT_UNKNOWN */
                    assert(0);
                    return false;
            }

            return _property_end();
        }

        /* Read ***************************************************/

        bool read_begin(const Element& element) { return _element_begin(element); }

        template<class T>
        bool read(const Property& property, T& value)
        {
            if (!_property_begin(property))
                return false;

            const Type type = property.type();

            if (!_type_read(format(), native_format(), type, _is, value))
                return false;

            return _property_end();
        }

        template<class T>
        bool read_count(const Property& property, T& count)
        {
            assert(property.type() == PLY_TYPE_LIST);

            if (!_property_begin(property))
                return false;

            const Type count_type = property.count_type();

            if (!_type_read(format(), native_format(), count_type, _is, count))
                return false;

            _current_property_list_count = count;
            _current_property_list_index = 0;

            return true;
        }

        template<class T>
        bool read_value(const Property& property, T& value)
        {
            assert(property.type() == PLY_TYPE_LIST);

            if (_current_property_list_index >= _current_property_list_count) {
                std::cerr << "Ply: reading too many list property values"
                          << std::endl;
                return false;
            }

            const Type value_type = property.value_type();

            if (!_type_read(format(), native_format(), value_type, _is, value))
                return false;

            _current_property_list_index++;

            return true;
        }

        bool read_end() { return _element_end(); }

        /* Write **************************************************/

        bool write_begin(const Element& element) { return _element_begin(element); }

        template<class T>
        bool write(const Property& property, T value)
        {
            if (!_property_begin(property))
                return false;

            if (format() == PLY_FORMAT_ASC && _current_property_num != 0)
                if ((_os << ' ') == 0) {
                    std::cerr << "Ply: I/O error" << std::endl;
                    return false;
                }

            const Type type = property.type();

            if (!_type_write(format(), native_format(), type, _os, value))
                return false;

            return _property_end();
        }

        template<class T>
        bool write_count(const Property& property, T count)
        {
            assert(property.type() == PLY_TYPE_LIST);

            if (!_property_begin(property))
                return false;

            if (format() == PLY_FORMAT_ASC && _current_property_num != 0)
                if ((_os << ' ') == 0) {
                    std::cerr << "Ply: I/O error" << std::endl;
                    return false;
                }

            const Type count_type = property.count_type();

            _current_property_list_count = count;
            _current_property_list_index = 0;

            return _type_write(format(), native_format(), count_type, _os, count);
        }

        template<class T>
        bool write_value(const Property& property, T value)
        {
            assert(property.type() == PLY_TYPE_LIST);

            if (format() == PLY_FORMAT_ASC)
                if ((_os << ' ') == 0) {
                    std::cerr << "Ply: I/O error" << std::endl;
                    return false;
                }

            if (_current_property_list_index >= _current_property_list_count) {
                std::cerr << "Ply: writing too many list property values"
                          << std::endl;
                return false;
            }

            const Type value_type = property.value_type();

            if (!_type_write(format(), native_format(), value_type, _os, value))
                return false;

            _current_property_list_index++;

            return true;
        }

        bool write_end()
        {
            if (format() == PLY_FORMAT_ASC)
                if ((_os << std::endl) == 0) {
                    std::cerr << "Ply: I/O error" << std::endl;
                    return false;
                }

            return _element_end();
        }

        /* Format *************************************************/

        Format native_format() const { return _native_format; }

        Format format() const { return _format; }

        float version() const { return _version; }

        /* Comments ***********************************************/

        size_type num_comments() const { return _comments.size(); }

        CommentsIterator comments_begin() const { return _comments.begin(); }

        CommentsIterator comments_end() const { return _comments.end(); }

        /* Object information *************************************/

        size_type num_obj_infos() const { return _obj_infos.size(); }

        ObjInfosIterator obj_infos_begin() const { return _obj_infos.begin(); }

        ObjInfosIterator obj_infos_end() const { return _obj_infos.end(); }

        /* Elements ***********************************************/

        /*! Returns a for the elements. */
        const ElementList& elements() const { return _elements; }

        size_type num_elements() const { return _elements.size(); }

        /*! Returns an iterator for the first element. */
        ElementsIterator elements_begin() const { return _elements.begin(); }

        /*! Returns a past-the-end iterator for the elements. */
        ElementsIterator elements_end() const { return _elements.end(); }

        const Element& element(size_type i) const { return _elements[i]; }

        /*! Finds an element with the given name. */
        ElementsIterator find_element(const string_type& name) const
        {
            ElementsIterator it;

            for (it = elements_begin(); it != elements_end(); ++it)
                if (it->name() == name)
                    return it;

            return it; /* = elements_end() */
        }

        /*! Checks whether an element with the given name exists. */
        bool has_element(const string_type& name) const { return (find_element(name) != elements_end()); }

        /*! Returns the element being read. */
        const Element& current_element() const { return element(_current_element_num); }

        /*! Returns the property being read. */
        const Property& current_property() const { return current_element().property(_current_property_num); }

        /* Miscellaneous ******************************************/

        /*! Converts the given format enum to the appropriate string. */
        static string_type format_enum_to_string(Format format)
        {
            switch (format) {
                case PLY_FORMAT_ASC:
                    return "ascii";
                case PLY_FORMAT_BIN_BE:
                    return "binary_big_endian";
                case PLY_FORMAT_BIN_LE:
                    return "binary_little_endian";
                default: /* PLY_FORMAT_UNKNOWN */
                    return "unknown";
            }
        }

        /*! Converts the given format string to the appropriate enum. */
        static Format format_string_to_enum(const string_type& s)
        {
            typedef std::pair<Format,string_type> FormatStringPair;

            static const FormatStringPair pairs[] = {
                FormatStringPair(PLY_FORMAT_ASC   , "ascii"               ),
                FormatStringPair(PLY_FORMAT_BIN_BE, "binary_big_endian"   ),
                FormatStringPair(PLY_FORMAT_BIN_LE, "binary_little_endian"),
            };
            static size_type num_pairs = sizeof(pairs)
                                       / sizeof(FormatStringPair);

            for (size_type i = 0; i != num_pairs; ++i)
                if (pairs[i].second == s)
                    return pairs[i].first;

            return PLY_FORMAT_UNKNOWN;
        }

        /*! Returns the size (in bytes) of the type or 0 if it a list. */
        static size_type type_size(Type type)
        {
            switch (type) {
                case PLY_TYPE_INT8: /* PLY_TYPE_CHAR */
                    return sizeof(int8_type);
                case PLY_TYPE_UINT8: /* PLY_TYPE_UCHAR */
                    return sizeof(uint8_type);
                case PLY_TYPE_INT16: /* PLY_TYPE_SHORT */
                    return sizeof(int16_type);
                case PLY_TYPE_UINT16: /* PLY_TYPE_USHORT */
                    return sizeof(uint16_type);
                case PLY_TYPE_INT32: /* PLY_TYPE_INT */
                    return sizeof(int32_type);
                case PLY_TYPE_UINT32: /* PLY_TYPE_UINT */
                    return sizeof(uint32_type);
                case PLY_TYPE_FLOAT32: /* PLY_TYPE_FLOAT */
                    return sizeof(float32_type);
                case PLY_TYPE_FLOAT64: /* PLY_TYPE_DOUBLE */
                    return sizeof(float64_type);
                default: /* PLY_TYPE_LIST, PLY_TYPE_UNKNOWN */
                    return 0;
            }
        }

        /*! Converts the given type enum to its string representation. */
        static string_type type_enum_to_string(Type type)
        {
            switch (type) {
                case PLY_TYPE_CHAR: /* PLY_TYPE_INT8 */
                    return "char";
                case PLY_TYPE_UCHAR: /* PLY_TYPE_UINT8 */
                    return "uchar";
                case PLY_TYPE_SHORT: /* PLY_TYPE_INT16 */
                    return "short";
                case PLY_TYPE_USHORT: /* PLY_TYPE_UINT16 */
                    return "ushort";
                case PLY_TYPE_INT: /* PLY_TYPE_INT32 */
                    return "int";
                case PLY_TYPE_UINT: /* PLY_TYPE_UINT32 */
                    return "uint";
                case PLY_TYPE_FLOAT: /* PLY_TYPE_FLOAT32 */
                    return "float";
                case PLY_TYPE_DOUBLE: /* PLY_TYPE_FLOAT64 */
                    return "double";
                case PLY_TYPE_LIST:
                    return "list";
                default: /* PLY_TYPE_UNKNOWN */
                    return "unknown";
            }
        }

        /*! Converts the given type string to the appropriate enum. */
        static Type type_string_to_enum(const string_type& s)
        {
            typedef std::pair<Type,string_type> TypeStringPair;

            static const TypeStringPair pairs[] = {
                TypeStringPair(PLY_TYPE_CHAR   , "char"   ),
                TypeStringPair(PLY_TYPE_INT8   , "int8"   ),
                TypeStringPair(PLY_TYPE_UCHAR  , "uchar"  ),
                TypeStringPair(PLY_TYPE_UINT8  , "uint8"  ),
                TypeStringPair(PLY_TYPE_SHORT  , "short"  ),
                TypeStringPair(PLY_TYPE_INT16  , "int16"  ),
                TypeStringPair(PLY_TYPE_USHORT , "ushort" ),
                TypeStringPair(PLY_TYPE_UINT16 , "uint16" ),
                TypeStringPair(PLY_TYPE_INT    , "int"    ),
                TypeStringPair(PLY_TYPE_INT32  , "int32"  ),
                TypeStringPair(PLY_TYPE_UINT   , "uint"   ),
                TypeStringPair(PLY_TYPE_UINT32 , "uint32" ),
                TypeStringPair(PLY_TYPE_FLOAT  , "float"  ),
                TypeStringPair(PLY_TYPE_FLOAT32, "float32"),
                TypeStringPair(PLY_TYPE_DOUBLE , "double" ),
                TypeStringPair(PLY_TYPE_FLOAT64, "float64"),
                TypeStringPair(PLY_TYPE_LIST   , "list"   ),
            };
            static size_type num_pairs = sizeof(pairs)
                                       / sizeof(TypeStringPair);

            for (size_type i = 0; i != num_pairs; ++i)
                if (pairs[i].second == s)
                    return pairs[i].first;

            return PLY_TYPE_UNKNOWN;
        }

    protected:
        bool _element_begin(const Element& element)
        {
            if (&current_element() != &element) {
                _current_element_num++;

                if (_current_element_num == num_elements()) {
                    std::cout << "Ply: no more element to read/write"
                              << std::endl;
                    return false;
                }

                if (&current_element() != &element) {
                    std::cout << "Ply: reading/writing wrong element: \""
                              << element.name() << "\" instead of \""
                              << current_element().name() << '"' << std::endl;
                    return false;
                }
            }

            _current_property_num = 0;

            return true;
        }

        bool _element_end()
        {
            const Element& current_element = this->current_element();

            const size_type num_properties = current_element.num_properties();

            if (_current_property_num < num_properties) {
                if (current_property().type() == PLY_TYPE_LIST)
                    return _property_end();
                else {
                    std::cout << "Ply: reading/writing incomplete \""
                              << current_element.name() << "\" element"
                              << std::endl;
                    return false;
                }
            }

            return true;
        }

        bool _property_begin(const Property& property)
        {
            const Property& current_property = this->current_property();

            if (&current_property != &property) {
                if (current_property.type() == PLY_TYPE_LIST) {
                    if (!_property_end())
                        return false;
                } else {
                    std::cout << "Ply: reading/writing wrong property: \""
                              << property.name() << "\" instead of \""
                              << current_property.name() << '"'
                              << std::endl;
                    return false;
                }
            }

            _current_property_list_index = 0;
            _current_property_list_count = 0;

            return true;
        }

        bool _property_end()
        {
            const Property& current_property = this->current_property();

            if (current_property.type() == PLY_TYPE_LIST
             && _current_property_list_index < _current_property_list_count) {
                std::cout << "Ply: reading/writing incomplete \""
                          << current_property.name() << "\" list property"
                          << std::endl;
                return false;
            }

            _current_property_num++;

            return true;
        }

        bool _skip_ASC(const Element& element);

        bool _skip_BIN(const Element& element);

        bool _skip_ASC(const Property& property);

        bool _skip_BIN(const Property& property);

        bool _check_last_element() const;

        bool _check_end_of_header() const;

        bool _parse_format(std::istringstream& iss, size_type line_num);

        void _parse_comment(std::istringstream& iss, size_type line_num);

        void _parse_obj_info(std::istringstream& iss, size_type line_num);

        bool _parse_element(std::istringstream& iss, size_type line_num);

        bool _parse_property(std::istringstream& iss, size_type line_num);

        /*! Returns the first element. */
        const Element& _first_element() const { return *_elements.begin(); }

        /*! Returns the last element. */
        const Element& _last_element() const { return *_elements.rbegin(); }

        /*! Returns the last element. */
        Element& _last_element() { return *_elements.rbegin(); }

        /*! Finds the byte order used by the machine. */
        void _find_native_format();

        /*!
         *  Reads one item of the given type and converts it to the
         *  supplied type.
         */
        template<class T>
        static bool _type_read(Format format, Format native_format, Type type,
                               std::istream& is, T& t)
        {
            switch (format) {
                case PLY_FORMAT_ASC:
                    return _type_read_ASC(type, is, t);
                case PLY_FORMAT_BIN_BE:
                    return _type_read_BIN_BE(native_format, type, is, t);
                case PLY_FORMAT_BIN_LE:
                    return _type_read_BIN_LE(native_format, type, is, t);
                default: /* PLY_FORMAT_UNKNOWN */
                    assert(0);
                    return false;
            }
        }

        /*! Converts one item to the supplied type and writes it. */
        template<class T>
        static bool _type_write(Format format, Format native_format, Type type,
                                std::ostream& os, T t)
        {
            switch (format) {
                case PLY_FORMAT_ASC:
                    return _type_write_ASC(type, os, t);
                case PLY_FORMAT_BIN_BE:
                    return _type_write_BIN_BE(native_format, type, os, t);
                case PLY_FORMAT_BIN_LE:
                    return _type_write_BIN_LE(native_format, type, os, t);
                default: /* PLY_FORMAT_UNKNOWN */
                    assert(0);
                    return false;
            }
        }

        /* ASCII **************************************************/

        template<class T>
        static bool _type_read_ASC(std::istream& is, T& t)
        {
            if (!(is >> t))
                return false;
            else
                return true;
        }

        template<class T>
        static bool _type_read_ASC(Type type, std::istream& is, T& t);

        template<class T>
        static bool _type_write_ASC(std::ostream& os, T t)
        {
            if (!(os << t))
                return false;
            else
                return true;
        }

        template<class T>
        static bool _type_write_ASC(Type type, std::ostream& os, T t);

        /* Binary *************************************************/

        /* Byte-swapping **********************/

        template<class T>
        static void _swab(T& t)
        {
            const size_type size = sizeof(T);
            const size_type half_size = size / 2;

            uint8_type* p = reinterpret_cast<uint8_type*>(&t);

            for (size_type i = 0; i != half_size; ++i)
                std::swap(p[i], p[size - i - 1]);
        }

        template<class T>
        static void _swab(T* p, size_type n)
        {
            const size_type size = sizeof(T);
            const size_type half_size = size / 2;

            uint8_type* q = reinterpret_cast<uint8_type*>(p);

            for (size_type i = 0; i != n; ++i)
                for (size_type j = 0; j != half_size; ++j)
                    std::swap(q[j], q[size - j - 1]);
        }

        /* Binary big-endian ******************/

        template<class T>
        static bool _type_read_BIN_BE(Format native_format,
                                      std::istream& is, T& t)
        {
            const size_type size = sizeof(T);

            char* p = reinterpret_cast<char*>(&t);

            if (!is.read(p, size))
                return false;

            if (native_format != PLY_FORMAT_BIN_BE)
                _swab(t);

            return true;
        }

        template<class T>
        static bool _type_read_BIN_BE(Format native_format, Type type,
                                     std::istream& is, T& t);

        template<class T>
        static bool _type_write_BIN_BE(Format native_format,
                                       std::ostream& os, T t)
        {
            if (native_format != PLY_FORMAT_BIN_BE)
                _swab(t);

            const size_type size = sizeof(T);

            const char* p = reinterpret_cast<const char*>(&t);

            if (!os.write(p, size))
                return false;

            return true;
        }

        template<class T>
        static bool _type_write_BIN_BE(Format native_format, Type type,
                                       std::ostream& os, T t);

        /* Binary little-endian ***************/

        template<class T>
        static bool _type_read_BIN_LE(Format native_format,
                                      std::istream& is, T& t)
        {
            const size_type size = sizeof(T);

            char* p = reinterpret_cast<char*>(&t);

            if (!is.read(p, size))
                return false;

            if (native_format != PLY_FORMAT_BIN_LE)
                _swab(t);

            return true;
        }

        template<class T>
        static bool _type_read_BIN_LE(Format native_format, Type type,
                                      std::istream& is, T& t);

        template<class T>
        static bool _type_write_BIN_LE(Format native_format,
                                       std::ostream& os, T t)
        {
            if (native_format != PLY_FORMAT_BIN_LE)
                _swab(t);

            const size_type size = sizeof(T);

            const char* p = reinterpret_cast<const char*>(&t);

            if (!os.write(p, size))
                return false;

            return true;
        }

        template<class T>
        static bool _type_write_BIN_LE(Format native_format, Type type,
                                       std::ostream& os, T t);

    protected:
        Format _format;
        Format _native_format;
        float _version;
        CommentList _comments;
        ObjInfoList _obj_infos;
        ElementList _elements;

        std::ifstream _is;
        std::ofstream _os;

        size_type _current_element_num;
        size_type _current_property_num;
        size_type _current_property_list_index;
        size_type _current_property_list_count;
};

template<>
inline bool Ply::_type_read_ASC<Ply::int8_type>(std::istream& is,
                                                int8_type& t)
{
#undef min
#undef max
    const int8_type min = std::numeric_limits<int8_type>::min();
    const int8_type max = std::numeric_limits<int8_type>::max();

    int i;

    if (!(is >> i))
        return false;

    if (i < min || i > max)
        return false;

    t = int8_type(i);

    return true;
}

template<>
inline bool Ply::_type_write_ASC<Ply::int8_type>(std::ostream& os,
                                                 int8_type t)
{
    return ((os << int32_type(t)) != 0);
}

template<>
inline bool Ply::_type_read_ASC<Ply::uint8_type>(std::istream& is,
                                                 uint8_type& t)
{
 // const uint8_type min = std::numeric_limits<uint8_type>::min();
    const uint8_type max = std::numeric_limits<uint8_type>::max();

    unsigned int u;

    if (!(is >> u))
        return false;

    if (/* u < min || */ u > max)
        return false;

    t = uint8_type(u);

    return true;
}

template<>
inline bool Ply::_type_write_ASC<Ply::uint8_type>(std::ostream& os,
                                                  uint8_type t)
{
    return ((os << uint32_type(t)) != 0);
}

template<class T>
inline bool Ply::_type_read_ASC(Type type, std::istream& is, T& t)
{
    bool b;

    int8_type i8 = 0;
    uint8_type u8 = 0;
    int16_type i16 = 0;
    uint16_type u16 = 0;
    int32_type i32 = 0;
    uint32_type u32 = 0;
    float32_type f32 = 0;
    float64_type f64 = 0;

    switch (type) {
        case PLY_TYPE_INT8: /* PLY_TYPE_CHAR */
            b = _type_read_ASC<int8_type>(is, i8);
            t = T(i8);
            break;
        case PLY_TYPE_UINT8: /* PLY_TYPE_UCHAR */
            b = _type_read_ASC<uint8_type>(is, u8);
            t = T(u8);
            break;
        case PLY_TYPE_INT16: /* PLY_TYPE_SHORT */
            b = _type_read_ASC<int16_type>(is, i16);
            t = T(i16);
            break;
        case PLY_TYPE_UINT16: /* PLY_TYPE_USHORT */
            b = _type_read_ASC<uint16_type>(is, u16);
            t = T(u16);
            break;
        case PLY_TYPE_INT32: /* PLY_TYPE_INT */
            b = _type_read_ASC<int32_type>(is, i32);
            t = T(i32);
            break;
        case PLY_TYPE_UINT32: /* PLY_TYPE_UINT */
            b = _type_read_ASC<uint32_type>(is, u32);
            t = T(u32);
            break;
        case PLY_TYPE_FLOAT32: /* PLY_TYPE_FLOAT */
            b = _type_read_ASC<float32_type>(is, f32);
            t = T(f32);
            break;
        case PLY_TYPE_FLOAT64: /* PLY_TYPE_DOUBLE */
            b = _type_read_ASC<float64_type>(is, f64);
            t = T(f64);
            break;
        default: /* PLY_TYPE_LIST, PLY_TYPE_UNKNOWN */
            assert(0);
            return false;
            break;
    }

    return b;
}

template<class T>
inline bool Ply::_type_write_ASC(Type type, std::ostream& os, T t)
{
    bool b;

    int8_type i8 = 0;
    uint8_type u8 = 0;
    int16_type i16 = 0;
    uint16_type u16 = 0;
    int32_type i32 = 0;
    uint32_type u32 = 0;
    float32_type f32 = 0;
    float64_type f64 = 0;

    switch (type) {
        case PLY_TYPE_INT8: /* PLY_TYPE_CHAR */
            i8 = int8_type(t);
            b = _type_write_ASC<int8_type>(os, i8);
            break;
        case PLY_TYPE_UINT8: /* PLY_TYPE_UCHAR */
            u8 = uint8_type(t);
            b = _type_write_ASC<uint8_type>(os, u8);
            break;
        case PLY_TYPE_INT16: /* PLY_TYPE_SHORT */
            i16 = int16_type(t);
            b = _type_write_ASC<int16_type>(os, i16);
            break;
        case PLY_TYPE_UINT16: /* PLY_TYPE_USHORT */
            u16 = uint16_type(t);
            b = _type_write_ASC<uint16_type>(os, u16);
            break;
        case PLY_TYPE_INT32: /* PLY_TYPE_INT */
            i32 = int32_type(t);
            b = _type_write_ASC<int32_type>(os, i32);
            break;
        case PLY_TYPE_UINT32: /* PLY_TYPE_UINT */
            u32 = uint32_type(t);
            b = _type_write_ASC<uint32_type>(os, u32);
            break;
        case PLY_TYPE_FLOAT32: /* PLY_TYPE_FLOAT */
            f32 = float32_type(t);
            b = _type_write_ASC<float32_type>(os, f32);
            break;
        case PLY_TYPE_FLOAT64: /* PLY_TYPE_DOUBLE */
            f64 = float64_type(t);
            b = _type_write_ASC<float64_type>(os, f64);
            break;
        default: /* PLY_TYPE_LIST, PLY_TYPE_UNKNOWN */
            assert(0);
            return false;
            break;
    }

    return b;
}

template<class T>
inline bool Ply::_type_read_BIN_BE(Format native_format, Type type,
                                   std::istream& is, T& t)
{
    bool b;

    int8_type i8 = 0;
    uint8_type u8 = 0;
    int16_type i16 = 0;
    uint16_type u16 = 0;
    int32_type i32 = 0;
    uint32_type u32 = 0;
    float32_type f32 = 0;
    float64_type f64 = 0;

    switch (type) {
        case PLY_TYPE_INT8: /* PLY_TYPE_CHAR */
            b = _type_read_BIN_BE<int8_type>(native_format, is, i8);
            t = T(i8);
            break;
        case PLY_TYPE_UINT8: /* PLY_TYPE_UCHAR */
            b = _type_read_BIN_BE<uint8_type>(native_format, is, u8);
            t = T(u8);
            break;
        case PLY_TYPE_INT16: /* PLY_TYPE_SHORT */
            b = _type_read_BIN_BE<int16_type>(native_format, is, i16);
            t = T(i16);
            break;
        case PLY_TYPE_UINT16: /* PLY_TYPE_USHORT */
            b = _type_read_BIN_BE<uint16_type>(native_format, is, u16);
            t = T(u16);
            break;
        case PLY_TYPE_INT32: /* PLY_TYPE_INT */
            b = _type_read_BIN_BE<int32_type>(native_format, is, i32);
            t = T(i32);
            break;
        case PLY_TYPE_UINT32: /* PLY_TYPE_UINT */
            b = _type_read_BIN_BE<uint32_type>(native_format, is, u32);
            t = T(u32);
            break;
        case PLY_TYPE_FLOAT32: /* PLY_TYPE_FLOAT */
            b = _type_read_BIN_BE<float32_type>(native_format, is, f32);
            t = T(f32);
            break;
        case PLY_TYPE_FLOAT64: /* PLY_TYPE_DOUBLE */
            b = _type_read_BIN_BE<float64_type>(native_format, is, f64);
            t = T(f64);
            break;
        default: /* PLY_TYPE_LIST, PLY_TYPE_UNKNOWN */
            assert(0);
            return false;
            break;
    }

    return b;
}

template<class T>
inline bool Ply::_type_write_BIN_BE(Format native_format, Type type,
                                    std::ostream& os, T t)
{
    bool b;

    int8_type i8 = 0;
    uint8_type u8 = 0;
    int16_type i16 = 0;
    uint16_type u16 = 0;
    int32_type i32 = 0;
    uint32_type u32 = 0;
    float32_type f32 = 0;
    float64_type f64 = 0;

    switch (type) {
        case PLY_TYPE_INT8: /* PLY_TYPE_CHAR */
            i8 = int8_type(t);
            b = _type_write_BIN_BE<int8_type>(native_format, os, i8);
            break;
        case PLY_TYPE_UINT8: /* PLY_TYPE_UCHAR */
            u8 = uint8_type(t);
            b = _type_write_BIN_BE<uint8_type>(native_format, os, u8);
            break;
        case PLY_TYPE_INT16: /* PLY_TYPE_SHORT */
            i16 = int16_type(t);
            b = _type_write_BIN_BE<int16_type>(native_format, os, i16);
            break;
        case PLY_TYPE_UINT16: /* PLY_TYPE_USHORT */
            u16 = uint16_type(t);
            b = _type_write_BIN_BE<uint16_type>(native_format, os, u16);
            break;
        case PLY_TYPE_INT32: /* PLY_TYPE_INT */
            i32 = int32_type(t);
            b = _type_write_BIN_BE<int32_type>(native_format, os, i32);
            break;
        case PLY_TYPE_UINT32: /* PLY_TYPE_UINT */
            u32 = uint32_type(t);
            b = _type_write_BIN_BE<uint32_type>(native_format, os, u32);
            break;
        case PLY_TYPE_FLOAT32: /* PLY_TYPE_FLOAT */
            f32 = float32_type(t);
            b = _type_write_BIN_BE<float32_type>(native_format, os, f32);
            break;
        case PLY_TYPE_FLOAT64: /* PLY_TYPE_DOUBLE */
            f64 = float64_type(t);
            b = _type_write_BIN_BE<float64_type>(native_format, os, f64);
            break;
        default: /* PLY_TYPE_LIST, PLY_TYPE_UNKNOWN */
            assert(0);
            return false;
            break;
    }

    return b;
}

template<class T>
inline bool Ply::_type_read_BIN_LE(Format native_format, Type type,
                                   std::istream& is, T& t)
{
    bool b;

    int8_type i8 = 0;
    uint8_type u8 = 0;
    int16_type i16 = 0;
    uint16_type u16 = 0;
    int32_type i32 = 0;
    uint32_type u32 = 0;
    float32_type f32 = 0;
    float64_type f64 = 0;

    switch (type) {
        case PLY_TYPE_INT8: /* PLY_TYPE_CHAR */
            b = _type_read_BIN_LE<int8_type>(native_format, is, i8);
            t = T(i8);
            break;
        case PLY_TYPE_UINT8: /* PLY_TYPE_UCHAR */
            b = _type_read_BIN_LE<uint8_type>(native_format, is, u8);
            t = T(u8);
            break;
        case PLY_TYPE_INT16: /* PLY_TYPE_SHORT */
            b = _type_read_BIN_LE<int16_type>(native_format, is, i16);
            t = T(i16);
            break;
        case PLY_TYPE_UINT16: /* PLY_TYPE_USHORT */
            b = _type_read_BIN_LE<uint16_type>(native_format, is, u16);
            t = T(u16);
            break;
        case PLY_TYPE_INT32: /* PLY_TYPE_INT */
            b = _type_read_BIN_LE<int32_type>(native_format, is, i32);
            t = T(i32);
            break;
        case PLY_TYPE_UINT32: /* PLY_TYPE_UINT */
            b = _type_read_BIN_LE<uint32_type>(native_format, is, u32);
            t = T(u32);
            break;
        case PLY_TYPE_FLOAT32: /* PLY_TYPE_FLOAT */
            b = _type_read_BIN_LE<float32_type>(native_format, is, f32);
            t = T(f32);
            break;
        case PLY_TYPE_FLOAT64: /* PLY_TYPE_DOUBLE */
            b = _type_read_BIN_LE<float64_type>(native_format, is, f64);
            t = T(f64);
            break;
        default: /* PLY_TYPE_LIST, PLY_TYPE_UNKNOWN */
            assert(0);
            return false;
            break;
    }

    return b;
}

template<class T>
inline bool Ply::_type_write_BIN_LE(Format native_format, Type type,
                                    std::ostream& os, T t)
{
    bool b;

    int8_type i8 = 0;
    uint8_type u8 = 0;
    int16_type i16 = 0;
    uint16_type u16 = 0;
    int32_type i32 = 0;
    uint32_type u32 = 0;
    float32_type f32 = 0;
    float64_type f64 = 0;

    switch (type) {
        case PLY_TYPE_INT8: /* PLY_TYPE_CHAR */
            i8 = int8_type(t);
            b = _type_write_BIN_LE<int8_type>(native_format, os, i8);
            break;
        case PLY_TYPE_UINT8: /* PLY_TYPE_UCHAR */
            u8 = uint8_type(t);
            b = _type_write_BIN_LE<uint8_type>(native_format, os, u8);
            break;
        case PLY_TYPE_INT16: /* PLY_TYPE_SHORT */
            i16 = int16_type(t);
            b = _type_write_BIN_LE<int16_type>(native_format, os, i16);
            break;
        case PLY_TYPE_UINT16: /* PLY_TYPE_USHORT */
            u16 = uint16_type(t);
            b = _type_write_BIN_LE<uint16_type>(native_format, os, u16);
            break;
        case PLY_TYPE_INT32: /* PLY_TYPE_INT */
            i32 = int32_type(t);
            b = _type_write_BIN_LE<int32_type>(native_format, os, i32);
            break;
        case PLY_TYPE_UINT32: /* PLY_TYPE_UINT */
            u32 = uint32_type(t);
            b = _type_write_BIN_LE<uint32_type>(native_format, os, u32);
            break;
        case PLY_TYPE_FLOAT32: /* PLY_TYPE_FLOAT */
            f32 = float32_type(t);
            b = _type_write_BIN_LE<float32_type>(native_format, os, f32);
            t = T(f32);
            break;
        case PLY_TYPE_FLOAT64: /* PLY_TYPE_DOUBLE */
            f64 = float64_type(t);
            b = _type_write_BIN_LE<float64_type>(native_format, os, f64);
            break;
        default: /* PLY_TYPE_LIST, PLY_TYPE_UNKNOWN */
            assert(0);
            return false;
            break;
    }

    return b;
}

#endif /* PLY_HPP */
