#include <JsonBox/Escaper.h>

#include <JsonBox/Grammar.h>
#include <JsonBox/Value.h>

namespace JsonBox {
	Escaper::Escaper() : afterBackSlash(false), inString(false) {
	}

	std::streambuf::int_type Escaper::operator()(std::streambuf &destination,
	                                             std::streambuf::int_type character) {
		bool notEscaped = true;
		std::streambuf::char_type tmpChar = std::streambuf::traits_type::to_char_type(character);

		// If we encounter a quotation mark.
		if (tmpChar == Structural::BEGIN_END_STRING) {
			// If we're not in a string, we change that. If we're in a string,
			// we change that only if we're not after an escape back slash.
			inString = !inString || (afterBackSlash);

		} else if (inString && !afterBackSlash) {
			// If we are in a string definition and we're not after a backslash
			// escape.
			if (tmpChar == Strings::Std::REVERSE_SOLIDUS) {
				destination.sputn(Strings::Json::REVERSE_SOLIDUS.c_str(), Strings::Json::REVERSE_SOLIDUS.size());
				notEscaped = false;

			} else if (tmpChar == Strings::Std::BACKSPACE) {
				destination.sputn(Strings::Json::BACKSPACE.c_str(), Strings::Json::BACKSPACE.size());
				notEscaped = false;

			} else if (tmpChar == Strings::Std::FORM_FEED) {
				destination.sputn(Strings::Json::FORM_FEED.c_str(), Strings::Json::FORM_FEED.size());
				notEscaped = false;

			} else if (tmpChar == Strings::Std::LINE_FEED) {
				destination.sputn(Strings::Json::LINE_FEED.c_str(), Strings::Json::LINE_FEED.size());
				notEscaped = false;

			} else if (tmpChar == Strings::Std::TAB) {
				destination.sputn(Strings::Json::TAB.c_str(), Strings::Json::TAB.size());
				notEscaped = false;

			} else if (tmpChar >= '\0' && tmpChar <= '\x1f') {
				std::string tmp(Value::escapeToUnicode(tmpChar));
				destination.sputn(tmp.c_str(), tmp.size());
				notEscaped = false;
			}

		}

		// We determine if we start a backslash escape or not.
		afterBackSlash = inString && !afterBackSlash && (tmpChar == Strings::Json::Escape::BEGIN_ESCAPE);
		return (notEscaped) ? (destination.sputc(tmpChar)) : (0);
	}
}
