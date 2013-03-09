#include "Ply.h"

bool Ply::open(const string_type& path, const ElementList& elements,
                                        const ObjInfoList& obj_infos,
                                        const CommentList& comments,
                                        Format format)
{
    _comments.clear();
    _obj_infos.clear();
    _elements.clear();

    if (format != PLY_FORMAT_ASC
     && format != PLY_FORMAT_BIN_BE
     && format != PLY_FORMAT_BIN_LE
     && format != PLY_FORMAT_UNKNOWN) {
        std::cerr << "Ply: unknown format" << std::endl;
        return false;
    }

    if (format == PLY_FORMAT_UNKNOWN)
        format = native_format();

    if (elements.size() == 0) {
        std::cerr << "Ply: no element given" << std::endl;
        return false;
    }

    for (ElementsIterator it = elements.begin();
         it != elements.end(); ++it)
        if (it->num_properties() == 0) {
            std::cerr << "Ply: element \"" << it->name()
                      << "\" has no property" << std::endl;
            return false;
        }

    _format = format;
    _version = 1.0f;

    _comments = comments;
    _obj_infos = obj_infos;
    _elements = elements;

    _current_element_num = 0;
    _current_property_num = 0;
    _current_property_list_index = 0;
    _current_property_list_count = 0;

    _os.open(path.c_str(), std::ios::out | std::ios::binary);

    if (!_os.is_open()) {
        std::cerr << "Ply: can't open file \"" << path << "\""
                  << std::endl;
        return false;
    }

    const int major = int(_version);
    const int minor = int(10 * (_version - major));

    _os << "ply" << std::endl
        << "format " << format_enum_to_string(_format)
        << ' ' << major << '.' << minor << std::endl;

    for (CommentsIterator it = comments_begin();
         it != comments_end(); ++it)
        _os << "comment " << *it << std::endl;

    for (ObjInfosIterator it = obj_infos_begin();
         it != obj_infos_end(); ++it)
        _os << "obj_info " << *it << std::endl;

    for (ElementsIterator it = elements_begin();
         it != elements_end(); ++it) {
        _os << "element " << it->name() << ' '
                          << it->count() << std::endl;

        for (PropertiesIterator it2 = it->properties_begin();
             it2 != it->properties_end(); ++it2) {
            const Type type = it2->type();

            _os << "property " << type_enum_to_string(type) << ' ';

            if (type == PLY_TYPE_LIST)
                _os << type_enum_to_string(it2->count_type()) << ' '
                    << type_enum_to_string(it2->value_type()) << ' ';

            _os << it2->name() << std::endl;
        }
    }

    _os << "end_header" << std::endl;

    return (_os != 0);
}

bool Ply::open(const string_type& path)
{
    _comments.clear();
    _obj_infos.clear();
    _elements.clear();

    _is.open(path.c_str(), std::ios::in | std::ios::binary);

    if (!_is.is_open()) {
        std::cerr << "Ply: can't open file \"" << path << "\""
                  << std::endl;
        return false;
    }

    string_type line_str;
    size_type line_num = 0;

    while (getline(_is, line_str)) {
        /* Extract tokens */
        std::istringstream iss(line_str);
        string_type token;

        if (!(iss >> token)) {
            std::cerr << "Ply: can't extract first token on line "
                      << (line_num + 1) << std::endl;
            _is.close();
            return false;
        }

        /* The first line must be "ply" */
        if (line_num == 0) {
            if (token != "ply") {
                std::cerr << "Ply: invalid first line in header"
                          << std::endl;
                _is.close();
                return false;
            }

            if (iss >> token)
                std::cerr << "Ply: skipping garbage in header"
                          << std::endl;

            line_num++;
            continue;
        }

        /* The end of the header */
        if (token == "end_header") {
            if (iss >> token)
                std::cerr << "Ply: skipping garbage in header"
                          << std::endl;

            if (!_check_end_of_header()) {
                _is.close();
                return false;
            } else
                break;
        }

        if (token == "format") {
            if (!_parse_format(iss, line_num)) {
                _is.close();
                return false;
            }
        } else if (token == "comment") {
            _parse_comment(iss, line_num);
        } else if (token == "obj_info") {
            _parse_obj_info(iss, line_num);
        } else if (token == "element") {
            if (!_parse_element(iss, line_num)) {
                _is.close();
                return false;
            }
        } else if (token == "property") {
            if (!_parse_property(iss, line_num)) {
                _is.close();
                return false;
            }
        } else {
            std::cerr << "Ply: unsupported token \"" << token
                      << "\" on line " << (line_num + 1) << std::endl;
            _is.close();
            return false;
        }

        line_num++;
    }

    _current_element_num = 0;
    _current_property_num = 0;
    _current_property_list_index = 0;
    _current_property_list_count = 0;

    return true;
}

void Ply::close()
{
    if (_os.is_open())
        _os.close();

    if (_is.is_open())
        _is.close();
}

bool Ply::_skip_ASC(const Element& element)
{
    const size_type count = element.count();

    for (size_type i = 0; i != count; ++i)
        for (PropertiesIterator it = element.properties_begin();
             it != element.properties_end(); ++it) {
            if (!_property_begin(*it))
                return false;

            if (!_skip_ASC(*it)) {
                std::cerr << "Ply: error while skipping property \""
                          << it->name() << "\" of element \""
                          << element.name() << "\" " << (i + 1)
                          << std::endl;
                return false;
            }

            if (!_property_end())
                return false;
        }

    return true;
}

bool Ply::_skip_BIN(const Element& element)
{
    const size_type count = element.count();
    const size_type size = element.size();

    if (size != 0) {
        /* Avoid calling _property_begin()/_property_end() for nothing */
        _current_property_num += element.num_properties();

        return (_is.seekg(count * size, std::ios::cur) != 0);
    }

    for (size_type i = 0; i != count; ++i)
        for (PropertiesIterator it = element.properties_begin();
             it != element.properties_end(); ++it) {
            if (!_property_begin(*it))
                return false;

            if (!_skip_BIN(*it)) {
                std::cerr << "Ply: error while skipping property \""
                          << it->name() << "\" of element \""
                          << element.name() << "\" " << (i + 1)
                          << std::endl;
                return false;
            }

            if (!_property_end())
                return false;
        }

    return true;
}

bool Ply::_skip_ASC(const Property& property)
{
    const Type type = property.type();

    assert(type != PLY_TYPE_UNKNOWN);

    if (type == PLY_TYPE_LIST) {
        const Type count_type = property.count_type();

        size_type count;

        if (!_type_read_ASC(count_type, _is, count))
            return false;

        if (count == 0)
            return true;

        const Type value_type = property.value_type();

        int dummy;

        for (size_type i = 0; i != count; ++i)
            if (!_type_read_ASC(value_type, _is, dummy))
                return false;
    } else {
        int dummy;

        if (!_type_read_ASC(type, _is, dummy))
            return false;
    }

    return true;
}

bool Ply::_skip_BIN(const Property& property)
{
    const Type type = property.type();

    assert(type != PLY_TYPE_UNKNOWN);

    if (type == PLY_TYPE_LIST) {
        const Type count_type = property.count_type();
        const Type value_type = property.value_type();

        size_type count = 0;

        _type_read(format(), native_format(), count_type, _is, count);

        const size_type value_type_size = type_size(value_type);

        return (_is.seekg(count * value_type_size, std::ios::cur) != 0);
    } else {
        const size_type size = type_size(type);

        return (_is.seekg(size, std::ios::cur) != 0);
    }
}

bool Ply::_check_last_element() const
{
    const Element& element = _last_element();

    if (element.num_properties() == 0) {
        std::cerr << "Ply: element \"" << element.name()
                  << "\" has no property" << std::endl;
        return false;
    } else
        return true;
}

bool Ply::_check_end_of_header() const
{
    if (format() == PLY_FORMAT_UNKNOWN) {
        std::cerr << "Ply: unspecified format" << std::endl;
        return false;
    }

    if (num_elements() == 0) {
        std::cerr << "Ply: no element found" << std::endl;
        return false;
    }

    if (!_check_last_element())
        return false;

    return true;
}

bool Ply::_parse_format(std::istringstream& iss, size_type line_num)
{
    if (line_num != 1) {
        std::cerr << "Ply: unexpected format token on line "
                  << (line_num + 1) << std::endl;
        return false;
    }

    string_type token;

    if (!(iss >> token)) {
        std::cerr << "Ply: can't extract format token on line "
                  << (line_num + 1) << std::endl;
        return false;
    }

    Format format = format_string_to_enum(token);

    if (format == PLY_FORMAT_UNKNOWN) {
        std::cerr << "Ply: unsupported format on line "
                  << (line_num + 1) << std::endl;
        return false;
    }

    if ((format == PLY_FORMAT_BIN_BE
      || format == PLY_FORMAT_BIN_LE)
     && native_format() == PLY_FORMAT_UNKNOWN) {
        std::cerr << "Ply: unknown native format" << std::endl;
        return false;
    }

    float version;

    if (!(iss >> version)) {
        std::cerr << "Ply: can't extract version token on line "
                  << (line_num + 1) << std::endl;
        return false;
    }

    if (iss >> token) {
        std::cerr << "Ply: unexpected token \"" << token
                  << "\" on line " << (line_num + 1) << std::endl;
        return false;
    }

    if (version != 1.0f) {
        std::cerr << "Ply: unsupported version on line "
                  << (line_num + 1) << std::endl;
        return false;
    }

    _format = format;
    _version = version;

    return true;
}

void Ply::_parse_comment(std::istringstream& iss, size_type line_num)
{
    iss >> std::noskipws;
    char c;
    iss >> c; /* Skip first char */

    string_type comment;

    while (iss >> c)
        comment += c;

    _comments.push_back(comment);
}

void Ply::_parse_obj_info(std::istringstream& iss, size_type line_num)
{
    iss >> std::noskipws;
    char c;
    iss >> c; /* Skip first char */

    string_type obj_info;

    while (iss >> c)
        obj_info += c;

    _obj_infos.push_back(obj_info);
}

bool Ply::_parse_element(std::istringstream& iss, size_type line_num)
{
    if (num_elements() != 0)
        if (!_check_last_element())
            return false;

    string_type name;

    if (!(iss >> name)) {
        std::cerr << "Ply: can't extract name token on line "
                  << (line_num + 1) << std::endl;
        return false;
    }

    if (has_element(name)) {
        std::cerr << "Ply: already specified element on line "
                  << (line_num + 1) << std::endl;
        return false;
    }

    size_type count;

    if (!(iss >> count)) {
        std::cerr << "Ply: can't extract count token on line "
                  << (line_num + 1) << std::endl;
        return false;
    }

    string_type token;

    if (iss >> token) {
        std::cout << "Ply: unexpected token \"" << token
                  << "\" on line " << (line_num + 1) << std::endl;
        return false;
    }

    _elements.push_back(Element(name, count));

    return true;
}

bool Ply::_parse_property(std::istringstream& iss, size_type line_num)
{
    if (num_elements() == 0) {
        std::cerr << "Ply: property specified w/o element on line "
                  << (line_num + 1) << std::endl;
        return false;
    }

    string_type token;

    if (!(iss >> token)) {
        std::cerr << "Ply: can't extract type token on line "
                  << (line_num + 1) << std::endl;
        return false;
    }

    Type type = type_string_to_enum(token);

    if (type == PLY_TYPE_UNKNOWN) {
        std::cerr << "Ply: unsupported type on line "
                  << (line_num + 1) << std::endl;
        return false;
    }

    Type count_type = PLY_TYPE_UNKNOWN;
    Type value_type = PLY_TYPE_UNKNOWN;

    if (type == PLY_TYPE_LIST) {
        if (!(iss >> token)) {
            std::cerr << "Ply: can't extract count type token"
                                 " on line " << (line_num + 1)
                      << std::endl;
            return false;
        }

        count_type = type_string_to_enum(token);

        if (count_type == PLY_TYPE_UNKNOWN) {
            std::cerr << "Ply: unsupported count type on line "
                      << (line_num + 1) << std::endl;
            return false;
        } else if (count_type == PLY_TYPE_CHAR
                || count_type == PLY_TYPE_SHORT
                || count_type == PLY_TYPE_INT
                || count_type == PLY_TYPE_FLOAT
                || count_type == PLY_TYPE_DOUBLE
                || count_type == PLY_TYPE_LIST) {
            std::cerr << "Ply: invalid count type ("
                      << type_enum_to_string(count_type)
                      << ") on line " << (line_num + 1) << std::endl;
            return false;
        }

        if (!(iss >> token)) {
            std::cerr << "Ply: can't extract value type token"
                             " on line " << (line_num + 1)
                      << std::endl;
            return false;
        }

        value_type = type_string_to_enum(token);

        if (value_type == PLY_TYPE_UNKNOWN) {
            std::cerr << "Ply: unsupported value type on line "
                      << (line_num + 1) << std::endl;
            return false;
        } else if (value_type == PLY_TYPE_LIST) {
            std::cerr << "Ply: invalid value type (list) on line "
                      << (line_num + 1) << std::endl;
            return false;
        }
    }

    if (!(iss >> token)) {
        std::cerr << "Ply: can't extract name token on line "
                  << (line_num + 1) << std::endl;
        return false;
    }

    /* Get the last element */
    Element& element = _last_element();

    if (element.has_property(token)) {
        std::cerr << "Ply: element \"" << element.name()
                  << "\" already has a property \"" << token << '"' 
                  << (line_num + 1) << std::endl;
        return false;
    }

    bool b;

    if (type == PLY_TYPE_LIST)
        b = element.add_property(Property(token, count_type,
                                                 value_type));
    else
        b = element.add_property(Property(token, type));

    return b;
}

void Ply::_find_native_format()
{
    const uint32_type u32 = 0x12345678;
    const uint8_type* ptr = reinterpret_cast<const uint8_type*>(&u32);

    if (ptr[0] == 0x12 && ptr[1] == 0x34
     && ptr[2] == 0x56 && ptr[3] == 0x78)
        _native_format = PLY_FORMAT_BIN_BE;
    else if (ptr[0] == 0x78 && ptr[1] == 0x56
          && ptr[2] == 0x34 && ptr[3] == 0x12)
        _native_format = PLY_FORMAT_BIN_LE;
    else
        _native_format = PLY_FORMAT_UNKNOWN;
}
