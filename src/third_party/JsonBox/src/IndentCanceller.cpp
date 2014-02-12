#include <JsonBox/IndentCanceller.h>

#include <JsonBox/Grammar.h>

namespace JsonBox {
	IndentCanceller::IndentCanceller() : afterBackSlash(false), inString(false) {
	}

	std::streambuf::int_type IndentCanceller::operator()(std::streambuf &destination,
	                                                     std::streambuf::int_type character) {
		std::streambuf::char_type tmpChar = std::streambuf::traits_type::to_char_type(character);

		// If we encounter a quotation mark.
		if (tmpChar == Structural::BEGIN_END_STRING) {
			// If we're not in a string, we change that. If we're in a string,
			// we change that only if we're not after an escape back slash.
			inString = !inString || (afterBackSlash);
		}

		// We determine if we start a backslash escape or not.
		afterBackSlash = inString && !afterBackSlash && (tmpChar == Strings::Json::Escape::BEGIN_ESCAPE);

		return (tmpChar != Whitespace::NEW_LINE && tmpChar != Whitespace::HORIZONTAL_TAB && tmpChar != Whitespace::CARRIAGE_RETURN && (inString || tmpChar != Whitespace::SPACE)) ? (destination.sputc(tmpChar)) : (0);
	}
}
