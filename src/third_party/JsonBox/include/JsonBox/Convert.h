/**
 * @file
 */
#ifndef JB_CONVERTER_H
#define JB_CONVERTER_H

#include <string>
#include <stdint.h>

namespace JsonBox {
	
	typedef std::basic_string<int32_t> String32;

	/**
	 * This class is used to encode/decode/transcode UTF8, 16 and 32.
	 */
	class Convert {
	public:
		/**
		 * Encode the given UTF32 string to a 8bit UTF8 one.
		 * @param utf32String UTF32 string to convert to UTF8.
		 * @return UTF8 string resulting from the conversion.
		 */
		static std::string encodeToUTF8(const String32& utf32String);

		/**
		 * Decode the given 8bit UTF8 string to an UTF32 string.
		 * @param utf8String UTF8 string to convert to UTF32.
		 * @return UTF32 string resulting from the conversion.
		 */
		static String32 decodeUTF8(const std::string& utf8String);
	};
}

#endif
