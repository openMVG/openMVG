#ifndef JB_INDENTER_H
#define JB_INDENTER_H

#include <streambuf>

namespace JsonBox {
	/**
	 * Adds a level of indentation to a streambuf.
	 * @see JsonBox::OutputFilter
	 */
	class Indenter {
	public:
		/**
		 * Default constructor.
		 */
		Indenter();

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
		/// Used to indicate if we are at the start of a new line.
		bool atStartOfLine;
	};
}

#endif
