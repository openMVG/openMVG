#ifndef JB_INDENT_CANCELLER_H
#define JB_INDENT_CANCELLER_H

#include <streambuf>

namespace JsonBox {
	/**
	 * Cancels indentations to a streambuf.
	 * @see JsonBox::OutputFilter
	 */
	class IndentCanceller {
	public:
		/**
		 * Default constructor.
		 */
		IndentCanceller();

		/**
		 * Inserts a tab character at the start of each line.
		 * @param destination Streambuf in which to insert the tab character.
		 * @param character Character to insert in the streambuf.
		 * @return Unspecified value not equal to traits::eof() on success,
		 * traits::eof() on failure.
		 */
		std::streambuf::int_type operator()(std::streambuf &destination,
		                                    std::streambuf::int_type character);
	private:
		bool afterBackSlash;
		bool inString;
	};
}

#endif
