#include <JsonBox/Indenter.h>

#include <JsonBox/Grammar.h>

namespace JsonBox {
	Indenter::Indenter() : atStartOfLine(true) {
	}

	std::streambuf::int_type Indenter::operator()(std::streambuf &destination,
	                                              std::streambuf::int_type character) {
		std::streambuf::char_type tmpChar = std::streambuf::traits_type::to_char_type(character);

		if (atStartOfLine && tmpChar != Whitespace::NEW_LINE) {
			destination.sputc(Whitespace::HORIZONTAL_TAB);
		}

		atStartOfLine = (tmpChar == Whitespace::NEW_LINE);
		return destination.sputc(tmpChar);
	}
}
