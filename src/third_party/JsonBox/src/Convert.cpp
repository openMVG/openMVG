#include <JsonBox/Convert.h>

#include <sstream>

#define         MASKBITS                0x3F //00111111
#define         MASK1BYTE               0x80 //10000000
#define         MASK2BYTES              0xC0 //11000000
#define         MASK3BYTES              0xE0 //11100000
#define         MASK4BYTES              0xF0 //11110000
#define         MASK5BYTES              0xF8 //11111000
#define         MASK6BYTES              0xFC //11111100

namespace JsonBox {
	std::string Convert::encodeToUTF8(const String32& utf32String) {
		std::stringstream result;

		for(String32::const_iterator i = utf32String.begin() ; i < utf32String.end(); ++i) {
			// 0xxxxxxx
			if(*i < 0x80) {
				result << static_cast<char>(*i);
			}
			// 110xxxxx 10xxxxxx
			else if(*i < 0x800) {
				result << static_cast<char>(MASK2BYTES | (*i >> 6));
				result << static_cast<char>(MASK1BYTE | (*i & MASKBITS));
			}
			// 1110xxxx 10xxxxxx 10xxxxxx
			else if(*i < 0x10000) {
				result << static_cast<char>(MASK3BYTES | (*i >> 12));
				result << static_cast<char>(MASK1BYTE | (*i >> 6 & MASKBITS));
				result << static_cast<char>(MASK1BYTE | (*i & MASKBITS));
			}
			// 11110xxx 10xxxxxx 10xxxxxx 10xxxxxx
			else if(*i < 0x200000) {
				result << static_cast<char>(MASK4BYTES | (*i >> 18));
				result << static_cast<char>(MASK1BYTE | (*i >> 12 & MASKBITS));
				result << static_cast<char>(MASK1BYTE | (*i >> 6 & MASKBITS));
				result << static_cast<char>(MASK1BYTE | (*i & MASKBITS));
			}
			// 111110xx 10xxxxxx 10xxxxxx 10xxxxxx 10xxxxxx
			else if(*i < 0x4000000) {
				result << static_cast<char>(MASK5BYTES | (*i >> 24));
				result << static_cast<char>(MASK1BYTE | (*i >> 18 & MASKBITS));
				result << static_cast<char>(MASK1BYTE | (*i >> 12 & MASKBITS));
				result << static_cast<char>(MASK1BYTE | (*i >> 6 & MASKBITS));
				result << static_cast<char>(MASK1BYTE | (*i & MASKBITS));
			}
			// 1111110x 10xxxxxx 10xxxxxx 10xxxxxx 10xxxxxx 10xxxxxx
			else if(*i < 0x8000000) {
				result << static_cast<char>(MASK6BYTES | (*i >> 30));
				result << static_cast<char>(MASK1BYTE | (*i >> 18 & MASKBITS));
				result << static_cast<char>(MASK1BYTE | (*i >> 12 & MASKBITS));
				result << static_cast<char>(MASK1BYTE | (*i >> 6 & MASKBITS));
				result << static_cast<char>(MASK1BYTE | (*i & MASKBITS));
			}
		}

		return result.str();
	}

	String32 Convert::decodeUTF8(const std::string& utf8String) {
		String32 result;

		String32::value_type tmpUnicodeChar;

		for(std::string::const_iterator i = utf8String.begin() ; i < utf8String.end();) {

			// 1111110x 10xxxxxx 10xxxxxx 10xxxxxx 10xxxxxx 10xxxxxx
			if((*i & MASK6BYTES) == MASK6BYTES) {
				tmpUnicodeChar = ((*i & 0x01) << 30) | ((*(i + 1) & MASKBITS) << 24)
				                 | ((*(i + 2) & MASKBITS) << 18) | ((*(i + 3)
				                         & MASKBITS) << 12)
				                 | ((*(i + 4) & MASKBITS) << 6) | (*(i + 5) & MASKBITS);
				i += 6;
			}
			// 111110xx 10xxxxxx 10xxxxxx 10xxxxxx 10xxxxxx
			else if((*i & MASK5BYTES) == MASK5BYTES) {
				tmpUnicodeChar = ((*i & 0x03) << 24) | ((*(i + 1)
				                                        & MASKBITS) << 18)
				                 | ((*(i + 2) & MASKBITS) << 12) | ((*(i + 3)
				                         & MASKBITS) << 6)
				                 | (*(i + 4) & MASKBITS);
				i += 5;
			}
			// 11110xxx 10xxxxxx 10xxxxxx 10xxxxxx
			else if((*i & MASK4BYTES) == MASK4BYTES) {
				tmpUnicodeChar = ((*i & 0x07) << 18) | ((*(i + 1)
				                                        & MASKBITS) << 12)
				                 | ((*(i + 2) & MASKBITS) << 6) | (*(i + 3) & MASKBITS);
				i += 4;
			}
			// 1110xxxx 10xxxxxx 10xxxxxx
			else if((*i & MASK3BYTES) == MASK3BYTES) {
				tmpUnicodeChar = ((*i & 0x0F) << 12) | ((*(i + 1) & MASKBITS) << 6)
				                 | (*(i + 2) & MASKBITS);
				i += 3;
			}
			// 110xxxxx 10xxxxxx
			else if((*i & MASK2BYTES) == MASK2BYTES) {
				tmpUnicodeChar = ((*i & 0x1F) << 6) | (*(i + 1) & MASKBITS);
				i += 2;
			}
			// 0xxxxxxx
			else { // if(*i < MASK1BYTE)
				tmpUnicodeChar = *i;
				i += 1;
			}

			result.push_back(tmpUnicodeChar);
		}

		return result;
	}
}
