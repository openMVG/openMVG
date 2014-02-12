#ifndef JB_OUTPUT_FILTER_H
#define JB_OUTPUT_FILTER_H

#include <cassert>

#include <limits>
#include <streambuf>

namespace JsonBox {
	/**
	 * Templated class used to filter output in an std::ostream. The custom
	 * mechanic of the filtering is easily implemented in the given Inserter. An
	 * output filter can be used to add indentation to an std::ostream, it can
	 * be used to log an std::ostream, etc.
	 * @tparam Inserter Type used as the inserter for the output filter.
	 */
	template <typename Inserter>
	class OutputFilter : public std::streambuf {
	public:
		/**
		 * Parameterized constructor.
		 * @param newDestination Pointer to the streambuf destination.
		 * @param newInserter Inserter to use to filter the output.
		 * @param newDeleteWhenFinished Used to know if the instance will have
		 * to delete its destination streambuf or not.
		 */
		OutputFilter(std::streambuf *newDestination, Inserter newInserter,
		             bool newDeleteWhenFinished = false) :
			destination(newDestination), inserter(newInserter),
			deleteWhenFinished(newDeleteWhenFinished) {
		}

		/**
		 * Parameterized constructor.
		 * @param newDestination Pointer to the streambuf destination.
		 * @param newInserter Inserter to use to filter the output.
		 */
		OutputFilter(std::streambuf *newDestination,
		             bool newDeleteWhenFinished = false) :
			destination(newDestination),
			deleteWhenFinished(newDeleteWhenFinished) {
		}

		/**
		 * Destructor. Takes care of deleting the destination streambuf if
		 * necessary.
		 */
		virtual ~OutputFilter() {
			// We delete the destination streambuf if necessary.
			if (deleteWhenFinished && destination) {
				delete destination;
			}
		}

		/**
		 * Actual function that calls the inserter to filter the output.
		 * @return Unspecified value not equal to traits::eof() on success,
		 * traits::eof() on failure.
		 */
		virtual int_type overflow(int_type ch) {
			int result = std::char_traits<char_type>::eof();

			// If the received character is invalid, we sync.
			if (ch == std::char_traits<char_type>::eof()) {
				result = sync();

			} else if (destination) {
				assert(ch >= 0 && ch <= static_cast<int_type>(std::numeric_limits<unsigned char>::max()));
				result = inserter(*destination, ch);
			}

			return result ;
		}

		/**
		 * Since it's an output filter, we don't need to do anything here.
		 */
		virtual int_type underflow() {
			return std::char_traits<char_type>::eof();
		}

		/**
		 * We don't need to do anything here. Calls the base class version.
		 * @return 0 on success, -1 otherwise. The base class version returns 0.
		 */
		virtual int_type sync() {
			return this->std::streambuf::sync();
		}

		/**
		 * We don't need to do anything here. Calls the base class version.
		 * @return Pointer to itself.
		 */
		virtual std::streambuf *setbuf(char *p, int len) {
			return this->std::streambuf::setbuf(p, len);
		}

		/**
		 * Gets the inserter's instance.
		 * @return Reference to the inserter.
		 */
		Inserter &getInserter() {
			return inserter;
		}

		/**
		 * Gets the destination streambuf.
		 * @return Pointer to the destination streambuf.
		 */
		std::streambuf *getDestination() const {
			return destination;
		}

	private:
		/// Pointer to the destination streambuf.
		std::streambuf *destination;

		/**
		 * Inserter to use to insert new characters in the destination
		 * streambuf.
		 * @see JsonBox::Filter<Inserter>::destination
		 */
		Inserter inserter;

		/**
		 * Bool used to know if the filter must delete its destination filter or
		 * not.
		 * @see JsonBox::Filter<Inserter>::destination
		 */
		bool deleteWhenFinished;
	};
}

#endif
