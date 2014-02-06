#ifndef JB_GRAMMAR_H
#define JB_GRAMMAR_H

#include <string>

namespace JsonBox {
	// Structural characters.
	namespace Structural {
		const char BEGIN_ARRAY = '[';
		const char BEGIN_OBJECT = '{';
		const char END_ARRAY = ']';
		const char END_OBJECT = '}';
		const char NAME_SEPARATOR = ':';
		const char VALUE_SEPARATOR = ',';
		const char BEGIN_END_STRING = '"';
	}
	
	// Whitespace characters.
	namespace Whitespace {
		const char SPACE = ' ';
		const char HORIZONTAL_TAB = '\t';
		const char NEW_LINE = '\n';
		const char CARRIAGE_RETURN = '\r';
	}
	
	// Literals.
	namespace Literals {
		const std::string FALSE_STRING = "false";
		const std::string TRUE_STRING = "true";
		const std::string NULL_STRING = "null";
	}
	
	// Numbers.
	namespace Numbers {
		const std::string DIGITS = "0123456789ABCDEFabcdef";
		const char DECIMAL_POINT = '.';
		const char LOWER_EXP = 'e';
		const char UPPER_EXP = 'E';
		const char MINUS = '-';
		const char PLUS = '+';
	}
	
	// Strings.
	namespace Strings {
		// C++ string characters.
		namespace Std {
			const char QUOTATION_MARK = '"';
			const char REVERSE_SOLIDUS = '\\';
			const char SOLIDUS = '/';
			const char BACKSPACE = '\b';
			const char FORM_FEED = '\f';
			const char LINE_FEED = '\n';
			const char CARRIAGE_RETURN = '\r';
			const char TAB = '\t';
		}
		// Json escape strings.
		namespace Json {
			const std::string QUOTATION_MARK = "\\\"";
			const std::string REVERSE_SOLIDUS = "\\\\";
			const std::string SOLIDUS = "\\/";
			const std::string BACKSPACE = "\\b";
			const std::string FORM_FEED = "\\f";
			const std::string LINE_FEED = "\\n";
			const std::string CARRIAGE_RETURN = "\\r";
			const std::string TAB = "\\t";
			const std::string BEGIN_UNICODE = "\\u";
			namespace Escape {
				const char BEGIN_ESCAPE = '\\';
				const char QUOTATION_MARK = '"';
				const char REVERSE_SOLIDUS = '\\';
				const char SOLIDUS = '/';
				const char BACKSPACE = 'b';
				const char FORM_FEED = 'f';
				const char LINE_FEED = 'n';
				const char CARRIAGE_RETURN = 'r';
				const char TAB = 't';
				const char BEGIN_UNICODE = 'u';
			}
		}
	}
}

#endif