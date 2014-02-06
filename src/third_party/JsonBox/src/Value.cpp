#include <JsonBox/Value.h>

#include <cassert>
#include <stack>
#include <sstream>
#include <list>
#include <iomanip>
#include <fstream>

#include <JsonBox/Grammar.h>
#include <JsonBox/Convert.h>
#include <JsonBox/OutputFilter.h>
#include <JsonBox/IndentCanceller.h>
#include <JsonBox/SolidusEscaper.h>
#include <JsonBox/Array.h>
#include <JsonBox/Object.h>

namespace JsonBox {

	const std::string Value::EMPTY_STRING = std::string();
	const double Value::EMPTY_DOUBLE = 0.0;
	const Object Value::EMPTY_OBJECT = Object();
	const Array Value::EMPTY_ARRAY = Array();

	std::string Value::escapeMinimumCharacters(const std::string &str) {
		std::stringstream result;

		// For each character in the string.
		for (std::string::const_iterator i = str.begin(); i != str.end(); ++i) {
			if (*i == Strings::Std::QUOTATION_MARK) {
				result << Strings::Json::QUOTATION_MARK;

			} else if (*i == Strings::Std::REVERSE_SOLIDUS) {
				result << Strings::Json::REVERSE_SOLIDUS;

			} else if (*i == Strings::Std::BACKSPACE) {
				result << Strings::Json::BACKSPACE;

			} else if (*i == Strings::Std::FORM_FEED) {
				result << Strings::Json::FORM_FEED;

			} else if (*i == Strings::Std::LINE_FEED) {
				result << Strings::Json::LINE_FEED;

			} else if (*i == Strings::Std::CARRIAGE_RETURN) {
				result << Strings::Json::CARRIAGE_RETURN;

			} else if (*i == Strings::Std::TAB) {
				result << Strings::Json::TAB;

			} else if (*i >= '\0' && *i <= '\x1f') {
				result << Value::escapeToUnicode(*i);

			} else {
				result << *i;
			}
		}

		return result.str();
	}

	std::string Value::escapeAllCharacters(const std::string &str) {
		std::stringstream result;

		// For each character in the string.
		for (std::string::const_iterator i = str.begin(); i != str.end(); ++i) {
			if (*i == Strings::Std::QUOTATION_MARK) {
				result << Strings::Json::QUOTATION_MARK;

			} else if (*i == Strings::Std::REVERSE_SOLIDUS) {
				result << Strings::Json::REVERSE_SOLIDUS;

			} else if (*i == Strings::Std::SOLIDUS) {
				result << Strings::Json::SOLIDUS;

			} else if (*i == Strings::Std::BACKSPACE) {
				result << Strings::Json::BACKSPACE;

			} else if (*i == Strings::Std::FORM_FEED) {
				result << Strings::Json::FORM_FEED;

			} else if (*i == Strings::Std::LINE_FEED) {
				result << Strings::Json::LINE_FEED;

			} else if (*i == Strings::Std::CARRIAGE_RETURN) {
				result << Strings::Json::CARRIAGE_RETURN;

			} else if (*i == Strings::Std::TAB) {
				result << Strings::Json::TAB;

			} else if (*i >= '\0' && *i <= '\x1f') {
				result << Value::escapeToUnicode(*i);

			} else {
				result << *i;
			}

		}

		return result.str();
	}


	const std::string Value::escapeToUnicode(char charToEscape) {
		std::stringstream result;

		if (charToEscape >= '\0' && charToEscape <= '\x1f') {
			result << "\\u00";
			result << std::hex << std::setfill('0') << std::setw(2) << static_cast<int>(charToEscape);
		}

		return result.str();
	}

	Value::Value() : type(NULL_VALUE), data() {
	}

	Value::Value(std::istream &input) : type(NULL_VALUE), data() {
		loadFromStream(input);
	}

	Value::Value(const std::string &newString) : type(STRING),
		data(new std::string(newString)) {
	}

	Value::Value(const char *newCString) : type(STRING),
		data(new std::string(newCString)) {
	}

	Value::Value(int newInt) : type(INTEGER), data(new int(newInt)) {
	}

	Value::Value(double newDouble) : type(DOUBLE), data(new double(newDouble)) {
	}

	Value::Value(const Object &newObject) : type(OBJECT),
		data(new Object(newObject)) {
	}

	Value::Value(const Array &newArray) : type(ARRAY),
		data(new Array(newArray)) {
	}

	Value::Value(bool newBoolean) : type(BOOLEAN), data(new bool(newBoolean)) {
	}

	Value::Value(const Value &src) : type(src.type), data() {
		switch (type) {
		case STRING:
			data.stringValue = new std::string(*src.data.stringValue);
			break;

		case INTEGER:
			data.intValue = new int(*src.data.intValue);
			break;

		case DOUBLE:
			data.doubleValue = new double(*src.data.doubleValue);
			break;

		case OBJECT:
			data.objectValue = new Object(*src.data.objectValue);
			break;

		case ARRAY:
			data.arrayValue = new Array(*src.data.arrayValue);
			break;

		case BOOLEAN:
			data.boolValue = new bool(*src.data.boolValue);
			break;

		default:
			type = NULL_VALUE;
			break;
		}
	}

	Value::~Value() {
		clear();
	}

	Value &Value::operator=(const Value &src) {
		if (this != &src) {
			clear();
			type = src.type;

			switch (type) {
			case STRING:
				data.stringValue = new std::string(*src.data.stringValue);
				break;

			case INTEGER:
				data.intValue = new int(*src.data.intValue);
				break;

			case DOUBLE:
				data.doubleValue = new double(*src.data.doubleValue);
				break;

			case OBJECT:
				data.objectValue = new Object(*src.data.objectValue);
				break;

			case ARRAY:
				data.arrayValue = new Array(*src.data.arrayValue);
				break;

			case BOOLEAN:
				data.boolValue = new bool(*src.data.boolValue);
				break;

			default:
				type = NULL_VALUE;
				data.stringValue = NULL;
				break;
			}
		}

		return *this;
	}

	bool Value::operator==(const Value &rhs) const {
		bool result = true;

		if (this != &rhs) {
			if (type == rhs.type) {
				switch (type) {
				case STRING:
					result = (*data.stringValue == *data.stringValue);
					break;

				case INTEGER:
					result = (*data.intValue == *data.intValue);
					break;

				case DOUBLE:
					result = (*data.doubleValue == *data.doubleValue);
					break;

				case OBJECT:
					result = (*data.objectValue == *data.objectValue);
					break;

				case ARRAY:
					result = (*data.arrayValue == *data.arrayValue);
					break;

				case BOOLEAN:
					result = (*data.boolValue == *data.boolValue);
					break;

				default:
					break;
				}

			} else {
				result = false;
			}
		}

		return result;
	}

	bool Value::operator!=(const Value &rhs) const {
		return !(*this == rhs);
	}

	bool Value::operator<(const Value &rhs) const {
		bool result = false;

		if (this != &rhs) {
			if (type == rhs.type) {
				switch (type) {
				case STRING:
					result = (*data.stringValue < *data.stringValue);
					break;

				case INTEGER:
					result = (*data.intValue < *data.intValue);
					break;

				case DOUBLE:
					result = (*data.doubleValue < *data.doubleValue);
					break;

				case OBJECT:
					result = (*data.objectValue < *data.objectValue);
					break;

				case ARRAY:
					result = (*data.arrayValue < *data.arrayValue);
					break;

				case BOOLEAN:
					result = (*data.boolValue < *data.boolValue);
					break;

				default:
					break;
				}

			}
		}

		return result;
	}

	bool Value::operator<=(const Value &rhs) const {
		return *this < rhs || *this == rhs;
	}

	bool Value::operator>(const Value &rhs) const {
		bool result = false;

		if (this != &rhs) {
			if (type == rhs.type) {
				switch (type) {
				case STRING:
					result = (*data.stringValue > *data.stringValue);
					break;

				case INTEGER:
					result = (*data.intValue > *data.intValue);
					break;

				case DOUBLE:
					result = (*data.doubleValue > *data.doubleValue);
					break;

				case OBJECT:
					result = (*data.objectValue > *data.objectValue);
					break;

				case ARRAY:
					result = (*data.arrayValue > *data.arrayValue);
					break;

				case BOOLEAN:
					result = (*data.boolValue	> *data.boolValue);
					break;

				default:
					break;
				}

			}
		}

		return result;
	}

	bool Value::operator>=(const Value &rhs) const {
		return *this > rhs || *this == rhs;
	}

	Value &Value::operator[](const Object::key_type &key) {
		if (type != OBJECT) {
			clear();
			type = OBJECT;
			data.objectValue = new Object();
		}

		return (*data.objectValue)[key];
	}

	Value &Value::operator[](const char *key) {
		return operator[](std::string(key));
	}

	Value &Value::operator[](Array::size_type index) {
		if (type != ARRAY) {
			clear();
			type = ARRAY;
			data.arrayValue = new Array(index + 1);
		}
        
        assert(index <= (*data.arrayValue).size());
        
        if(index == (*data.arrayValue).size())
        {
            (*data.arrayValue).push_back(Value());
            return  (*data.arrayValue).back();
        }
        else
        {
            return (*data.arrayValue)[index];
        }
	}

	Value::Type Value::getType() const {
		return type;
	}

	bool Value::isString() const {
		return type == STRING;
	}

	bool Value::isInteger() const {
		return type == INTEGER;
	}

	bool Value::isDouble() const {
		return type == DOUBLE;
	}

	bool Value::isObject() const {
		return type == OBJECT;
	}

	bool Value::isArray() const {
		return type == ARRAY;
	}

	bool Value::isBoolean() const {
		return type == BOOLEAN;
	}

	bool Value::isNull() const {
		return type == NULL_VALUE;
	}

	const std::string &Value::getString() const {
		return (type == STRING) ? (*data.stringValue) : (EMPTY_STRING);
	}

	void Value::setString(std::string const &newString) {
		if (type == STRING) {
			*data.stringValue = newString;

		} else {
			clear();
			type = STRING;
			data.stringValue = new std::string(newString);
		}
	}

	int Value::getInt() const {
		return (type == INTEGER) ? (*data.intValue) : (EMPTY_INT);
	}

	void Value::setInt(int newInt) {
		if (type == INTEGER) {
			*data.intValue = newInt;

		} else {
			clear();
			type = INTEGER;
			data.intValue = new int(newInt);
		}
	}

	double Value::getDouble() const {
		return (type == DOUBLE) ? (*data.doubleValue) : (EMPTY_DOUBLE);
	}

	void Value::setDouble(double newDouble) {
		if (type == DOUBLE) {
			*data.doubleValue = newDouble;

		} else {
			clear();
			type = DOUBLE;
			data.doubleValue = new double(newDouble);
		}
	}

	const Object &Value::getObject() const {
		return (type == OBJECT) ? (*data.objectValue) : (EMPTY_OBJECT);
	}

	void Value::setObject(const Object &newObject) {
		if (type == OBJECT) {
			*data.objectValue = newObject;

		} else {
			clear();
			type = OBJECT;
			data.objectValue = new Object(newObject);
		}
	}

	const Array &Value::getArray() const {
		return (type == ARRAY) ? (*data.arrayValue) : (EMPTY_ARRAY);
	}

	void Value::setArray(const Array &newArray) {
		if (type == ARRAY) {
			*data.arrayValue = newArray;

		} else {
			clear();
			type = ARRAY;
			data.arrayValue = new Array(newArray);
		}
	}

	bool Value::getBoolean() const {
		return (type == BOOLEAN) ? (*data.boolValue) : (EMPTY_BOOL);
	}

	void Value::setBoolean(bool newBoolean) {
		if (type == BOOLEAN) {
			*data.boolValue = newBoolean;

		} else {
			clear();
			type = BOOLEAN;
			data.boolValue = new bool(newBoolean);
		}
	}

	void Value::setNull() {
		clear();
		type = NULL_VALUE;
		data.stringValue = NULL;
	}

	void Value::loadFromString(std::string const &json) {
		std::stringstream jsonStream(json);
		loadFromStream(jsonStream);
	}

	void Value::loadFromStream(std::istream &input) {
		char currentCharacter;

		// We check that the stream is in UTF-8.
		char encoding[2];
		input.get(encoding[0]);
		input.get(encoding[1]);

		if (encoding[0] != '\0' && encoding[1] != '\0') {
			// We put the characters back.
			input.putback(encoding[1]);
			input.putback(encoding[0]);

			// Boolean value used to stop reading characters after the value
			// is done loading.
			bool noErrors = true;

			while (noErrors && input.good()) {
				input.get(currentCharacter);

				if (input.good()) {
					if (currentCharacter == Structural::BEGIN_END_STRING) {
						// The value to be parsed is a string.
						setString("");
						readString(input, *data.stringValue);
						noErrors = false;

					} else if (currentCharacter == Structural::BEGIN_OBJECT) {
						// The value to be parsed is an object.
						setObject(Object());
						readObject(input, *data.objectValue);
						noErrors = false;

					} else if (currentCharacter == Structural::BEGIN_ARRAY) {
						// The value to be parsed is an array.
						setArray(Array());
						readArray(input, *data.arrayValue);
						noErrors = false;

					} else if (currentCharacter == Literals::NULL_STRING[0]) {
						// We try to read the literal 'null'.
						if (!input.eof()) {
							input.get(currentCharacter);

							if (currentCharacter == Literals::NULL_STRING[1]) {
								if (!input.eof()) {
									input.get(currentCharacter);

									if (currentCharacter == Literals::NULL_STRING[2]) {
										if (!input.eof()) {
											input.get(currentCharacter);

											if (currentCharacter == Literals::NULL_STRING[3]) {
												setNull();
												noErrors = false;

											} else {
												std::cout << "invalid characters found" << std::endl;
											}

										} else {
											std::cout << "json input ends incorrectly" << std::endl;
										}

									} else {
										std::cout << "invalid characters found" << std::endl;
									}

								} else {
									std::cout << "json input ends incorrectly" << std::endl;
								}

							} else {
								std::cout << "invalid characters found" << std::endl;
							}

						} else {
							std::cout << "json input ends incorrectly" << std::endl;
						}

					} else if (currentCharacter == Numbers::MINUS ||
					           (currentCharacter >= Numbers::DIGITS[0] && currentCharacter <= Numbers::DIGITS[9])) {
						// Numbers can't start with zeroes.
						input.putback(currentCharacter);
						readNumber(input, *this);
						noErrors = false;

					} else if (currentCharacter == Literals::TRUE_STRING[0]) {
						// We try to read the boolean literal 'true'.
						if (!input.eof()) {
							input.get(currentCharacter);

							if (currentCharacter == Literals::TRUE_STRING[1]) {
								if (!input.eof()) {
									input.get(currentCharacter);

									if (currentCharacter == Literals::TRUE_STRING[2]) {
										if (!input.eof()) {
											input.get(currentCharacter);

											if (currentCharacter == Literals::TRUE_STRING[3]) {
												setBoolean(true);
												noErrors = false;
											}
										}
									}
								}
							}
						}

					} else if (currentCharacter == Literals::FALSE_STRING[0]) {
						// We try to read the boolean literal 'false'.
						if (!input.eof()) {
							input.get(currentCharacter);

							if (currentCharacter == Literals::FALSE_STRING[1]) {
								if (!input.eof()) {
									input.get(currentCharacter);

									if (currentCharacter == Literals::FALSE_STRING[2]) {
										if (!input.eof()) {
											input.get(currentCharacter);

											if (currentCharacter == Literals::FALSE_STRING[3]) {
												if (!input.eof()) {
													input.get(currentCharacter);

													if (currentCharacter == Literals::FALSE_STRING[4]) {
														setBoolean(false);
														noErrors = false;
													}
												}
											}
										}
									}
								}
							}
						}

					} else if (!isWhiteSpace(currentCharacter)) {
						std::cout << "Invalid character found: '" << currentCharacter << "'" << std::endl;
					}
				}
			}

		} else {
			std::cout << "File is not in UTF-8, not parsing." << std::endl;
		}
	}

	void Value::loadFromFile(const std::string &filePath) {
		std::ifstream file;
		file.open(filePath.c_str());

		if (file.is_open()) {
			loadFromStream(file);
			file.close();

		} else {
			std::cout << "Failed to open file to load the json: " << filePath << std::endl;
		}
	}

	void Value::writeToStream(std::ostream &output, bool indent,
	                          bool escapeAll) const {
		this->output(output, indent, escapeAll);
	}

	void Value::writeToFile(const std::string &filePath, bool indent,
	                        bool escapeAll) const {
		std::ofstream file;
		file.open(filePath.c_str());

		if (file.is_open()) {
			writeToStream(file, indent, escapeAll);
			file.close();

		} else {
			std::cout << "Failed to open file to write the json into: " << filePath << std::endl;
		}
	}

	Value::ValueDataPointer::ValueDataPointer(): stringValue(NULL) {
	}

	Value::ValueDataPointer::ValueDataPointer(std::string *newStringValue) :
		stringValue(newStringValue) {
	}

	Value::ValueDataPointer::ValueDataPointer(int *newIntValue) :
		intValue(newIntValue) {
	}

	Value::ValueDataPointer::ValueDataPointer(double *newDoubleValue) :
		doubleValue(newDoubleValue) {
	}

	Value::ValueDataPointer::ValueDataPointer(Object *newObjectValue) :
		objectValue(newObjectValue) {
	}

	Value::ValueDataPointer::ValueDataPointer(Array *newArrayValue) :
		arrayValue(newArrayValue) {
	}

	Value::ValueDataPointer::ValueDataPointer(bool *newBoolValue) :
		boolValue(newBoolValue) {
	}

	bool Value::isHexDigit(char digit) {
		return (digit >= Numbers::DIGITS[0] && digit <= Numbers::DIGITS[9]) || (digit >= Numbers::DIGITS[10] && digit <= Numbers::DIGITS[15]) ||
		       (digit >= Numbers::DIGITS[16] && digit <= Numbers::DIGITS[21]);
	}

	bool Value::isWhiteSpace(char whiteSpace) {
		return whiteSpace == Whitespace::SPACE ||
		       whiteSpace == Whitespace::HORIZONTAL_TAB ||
		       whiteSpace == Whitespace::NEW_LINE ||
		       whiteSpace == Whitespace::CARRIAGE_RETURN;
	}

	void Value::readString(std::istream &input, std::string &result) {
		bool noErrors = true, noUnicodeError = true;
		char currentCharacter, tmpCharacter;
		std::stringstream constructing;
		std::string tmpStr(4, ' ');
		std::stringstream tmpSs;
		int32_t tmpInt;
		String32 tmpStr32;
		unsigned int tmpCounter;

		// As long as there aren't any errors and that we haven't reached the
		// end of the input stream.
		while (noErrors && !input.eof()) {
			input.get(currentCharacter);

			if (input.good()) {
				if (currentCharacter & 0x80) { // 0x80 --> 10000000
					// The character is part of an utf8 character.
					constructing << currentCharacter;

				} else if (currentCharacter == Strings::Json::Escape::BEGIN_ESCAPE) {
					if (!input.eof()) {
						input.get(tmpCharacter);

						switch (tmpCharacter) {
						case Strings::Json::Escape::QUOTATION_MARK:
							constructing << Strings::Std::QUOTATION_MARK;
							break;

						case Strings::Json::Escape::REVERSE_SOLIDUS:
							constructing << Strings::Std::REVERSE_SOLIDUS;
							break;

						case Strings::Json::Escape::SOLIDUS:
							constructing << Strings::Std::SOLIDUS;
							break;

						case Strings::Json::Escape::BACKSPACE:
							constructing << Strings::Std::BACKSPACE;
							break;

						case Strings::Json::Escape::FORM_FEED:
							constructing << Strings::Std::FORM_FEED;
							break;

						case Strings::Json::Escape::LINE_FEED:
							constructing << Strings::Std::LINE_FEED;
							break;

						case Strings::Json::Escape::CARRIAGE_RETURN:
							constructing << Strings::Std::CARRIAGE_RETURN;
							break;

						case Strings::Json::Escape::TAB:
							constructing << Strings::Std::TAB;
							break;

						case Strings::Json::Escape::BEGIN_UNICODE:
							// TODO: Check for utf16 surrogate pairs.
							tmpCounter = 0;
							tmpStr.clear();
							tmpStr = "    ";
							noUnicodeError = true;

							while (tmpCounter < 4 && !input.eof()) {
								input.get(tmpCharacter);

								if (isHexDigit(tmpCharacter)) {
									tmpStr[tmpCounter] = tmpCharacter;

								} else {
									noUnicodeError = false;
									std::cout << "Invalid \\u character, skipping it." << std::endl;
								}

								++tmpCounter;
							}

							if (noUnicodeError) {
								tmpSs.clear();
								tmpSs.str("");
								tmpSs << std::hex << tmpStr;
								tmpSs >> tmpInt;
								tmpStr32.clear();
								tmpStr32.push_back(tmpInt);
								tmpStr = Convert::encodeToUTF8(tmpStr32);
								constructing << tmpStr;
							}

							break;

						default:
							break;
						}
					}

				} else if (currentCharacter == '"') {
					result = constructing.str();
					noErrors = false;

				} else {
					constructing << currentCharacter;
				}
			}
		}
	}

	void Value::readObject(std::istream &input, Object &result) {
		bool noErrors = true;
		char currentCharacter;
		std::string tmpString;

		while (noErrors && !input.eof()) {
			input.get(currentCharacter);

			if (input.good()) {
				if (currentCharacter == Structural::BEGIN_END_STRING) {
					// We read the object's member's name.
					readString(input, tmpString);
					currentCharacter = input.peek();
					// We read white spaces until the next non white space.
					readToNonWhiteSpace(input, currentCharacter);

					if (!input.eof()) {

						// We make sure it's the right character.
						if (currentCharacter == Structural::NAME_SEPARATOR) {
							// We read until the value starts.
							readToNonWhiteSpace(input, currentCharacter);

							if (!input.eof()) {
								// We put the character back and we load the value
								// from the stream.
								input.putback(currentCharacter);
								result[tmpString].loadFromStream(input);

								while (!input.eof() && currentCharacter != Structural::VALUE_SEPARATOR &&
								       currentCharacter != Structural::END_OBJECT) {
									input.get(currentCharacter);
								}

								if (currentCharacter == Structural::END_OBJECT) {
									// We are done reading the object.
									noErrors = false;
								}
							}
						}
					}

				} else if (currentCharacter == Structural::END_OBJECT) {
					noErrors = false;

				} else if (!isWhiteSpace(currentCharacter)) {
					std::cout << "Expected '\"', got '" << currentCharacter << "', ignoring it." << std::endl;
				}
			}
		}
	}

	void Value::readArray(std::istream &input, Array &result) {
		bool notDone = true;
		char currentChar;

		while (notDone && !input.eof()) {
			input.get(currentChar);

			if (input.good()) {
				if (currentChar == Structural::END_ARRAY) {
					notDone = false;

				} else if (!isWhiteSpace(currentChar)) {
					input.putback(currentChar);
					result.push_back(Value());
					result.back().type = UNKNOWN;
					result.back().loadFromStream(input);

					if (result.back().type == UNKNOWN) {
						result.pop_back();
					}

					while (!input.eof() && currentChar != ',' &&
					       currentChar != Structural::END_ARRAY) {
						input.get(currentChar);
					}

					if (currentChar == Structural::END_ARRAY) {
						notDone = false;
					}
				}
			}
		}
	}

	void Value::readNumber(std::istream &input, JsonBox::Value &result) {
		bool notDone = true, inFraction = false, inExponent = false;
		char currentCharacter;
		std::stringstream constructing;

		if (!input.eof() && input.peek() == Numbers::DIGITS[0]) {
			// We make sure there isn't more than one zero.
			input.get(currentCharacter);

			if (input.peek() == '0') {
				notDone = false;

			} else {
				input.putback(currentCharacter);
			}
		}

		while (notDone && !input.eof()) {
			input.get(currentCharacter);

			if (currentCharacter == '-') {
				if (constructing.str().empty()) {
					constructing << currentCharacter;

				} else {
					std::cout << "Expected a digit, '.', 'e' or 'E', got '" << currentCharacter << "' instead, ignoring it." << std::endl;
				}

			} else if (currentCharacter >= '0' && currentCharacter <= '9') {
				constructing << currentCharacter;

			} else if (currentCharacter == '.') {
				if (!inFraction && !inExponent) {
					inFraction = true;
					constructing << currentCharacter;
				}

			} else if (currentCharacter == 'e' || currentCharacter == 'E') {
				if (!inExponent) {
					inExponent = true;
					constructing << currentCharacter;

					if (!input.eof() && (input.peek() == '-' || input.peek() == '+')) {
						input.get(currentCharacter);
						constructing << currentCharacter;
					}
				}

			} else {
				input.putback(currentCharacter);
				notDone = false;
			}
		}

		if (inFraction || inExponent) {
			double doubleResult;
			constructing >> doubleResult;
			result.setDouble(doubleResult);

		} else {
			int intResult;
			constructing >> intResult;
			result.setInt(intResult);
		}
	}

	void Value::readToNonWhiteSpace(std::istream &input, char &currentCharacter) {
		do {
			input.get(currentCharacter);
		} while (!input.eof() && isWhiteSpace(currentCharacter));
	}

	void Value::clear() {
		switch (type) {
		case STRING:
			delete data.stringValue;
			break;

		case INTEGER:
			delete data.intValue;
			break;

		case DOUBLE:
			delete data.doubleValue;
			break;

		case OBJECT:
			delete data.objectValue;
			break;

		case ARRAY:
			delete data.arrayValue;
			break;

		case BOOLEAN:
			delete data.boolValue;
			break;

		default:
			break;
		}
	}

	void Value::output(std::ostream &output, bool indent,
	                   bool escapeAll) const {
		if (indent) {
			if (escapeAll) {
				OutputFilter<SolidusEscaper> solidusEscaper(output.rdbuf());
				output.rdbuf(&solidusEscaper);
				output << *this;
				output.rdbuf(solidusEscaper.getDestination());

			} else {
				output << *this;
			}

		} else {
			OutputFilter<IndentCanceller> indentCanceller(output.rdbuf());
			output.rdbuf(&indentCanceller);

			if (escapeAll) {
				OutputFilter<SolidusEscaper> solidusEscaper(output.rdbuf());
				output.rdbuf(&solidusEscaper);
				output << *this;
				output.rdbuf(solidusEscaper.getDestination());

			} else {
				output << *this;
			}

			output.rdbuf(indentCanceller.getDestination());
		}
	}

	std::ostream &operator<<(std::ostream &output, const Value &v) {
		switch (v.type) {
		case Value::STRING:
			output << Structural::BEGIN_END_STRING << Value::escapeMinimumCharacters(v.getString()) << Structural::BEGIN_END_STRING;
			break;

		case Value::INTEGER:
			output << v.getInt();
			break;

		case Value::DOUBLE:
			output << v.getDouble();
			break;

		case Value::OBJECT:
			output << v.getObject();
			break;

		case Value::ARRAY:
			output << v.getArray();
			break;

		case Value::BOOLEAN:
			output << (v.getBoolean() ? Literals::TRUE_STRING : Literals::FALSE_STRING);
			break;

		case Value::NULL_VALUE:
			output << Literals::NULL_STRING;
			break;

		default:
			break;
		}

		return output;
	}
}
