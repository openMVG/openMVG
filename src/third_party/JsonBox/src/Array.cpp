#include <JsonBox/Array.h>

#include <JsonBox/Value.h>
#include <JsonBox/Grammar.h>
#include <JsonBox/OutputFilter.h>
#include <JsonBox/Indenter.h>

namespace JsonBox {
	Array::Array(const allocator_type &alloc) : data(alloc) {
	}

	Array::Array(size_type count, const_reference value, const allocator_type &alloc) : data(count, value, alloc) {
	}

	Array::Array(const Array &other) : data(other.data) {
	}

	Array &Array::operator=(const Array &other) {
		data = other.data;
		return *this;
	}

	bool Array::operator==(const Array &rhs) const {
		return data == rhs.data;
	}

	bool Array::operator!=(const Array &rhs) const {
		return data != rhs.data;
	}

	bool Array::operator<(const Array &rhs) const {
		return data < rhs.data;
	}

	bool Array::operator<=(const Array &rhs) const {
		return data <= rhs.data;
	}

	bool Array::operator>(const Array &rhs) const {
		return data > rhs.data;
	}

	bool Array::operator>=(const Array &rhs) const {
		return data >= rhs.data;
	}

	void Array::assign(size_type count, const_reference value) {
		data.assign(count, value);
	}

	Array::allocator_type Array::get_allocator() const {
		return data.get_allocator();
	}

	Array::reference Array::at(size_type pos) {
		return data.at(pos);
	}

	Array::const_reference Array::at(size_type pos) const {
		return data.at(pos);
	}

	Array::reference Array::operator[](size_type pos) {
		return data[pos];
	}

	Array::const_reference Array::operator[](size_type pos) const {
		return data[pos];
	}

	Array::reference Array::front() {
		return data.front();
	}

	Array::const_reference Array::front() const {
		return data.front();
	}

	Array::reference Array::back() {
		return data.back();
	}

	Array::const_reference Array::back() const {
		return data.back();
	}

	Array::iterator Array::begin() {
		return data.begin();
	}

	Array::const_iterator Array::begin() const {
		return data.begin();
	}

	Array::iterator Array::end() {
		return data.end();
	}

	Array::const_iterator Array::end() const {
		return data.end();
	}

	Array::reverse_iterator Array::rbegin() {
		return data.rbegin();
	}

	Array::const_reverse_iterator Array::rbegin() const {
		return data.rbegin();
	}

	Array::reverse_iterator Array::rend() {
		return data.rend();
	}

	Array::const_reverse_iterator Array::rend() const {
		return data.rend();
	}

	bool Array::empty() const {
		return data.empty();
	}

	Array::size_type Array::size() const {
		return data.size();
	}

	Array::size_type Array::max_size() const {
		return data.max_size();
	}

	void Array::reserve(size_type size) {
		data.reserve(size);
	}

	Array::size_type Array::capacity() const {
		return data.capacity();
	}

	void Array::clear() {
		data.clear();
	}

	Array::iterator Array::insert(iterator pos, const_reference value) {
		return data.insert(pos, value);
	}

	void Array::insert(iterator pos, size_type count, const_reference value) {
		data.insert(pos, count, value);
	}

	Array::iterator Array::erase(iterator pos) {
		return data.erase(pos);
	}

	Array::iterator Array::erase(iterator first, iterator last) {
		return data.erase(first, last);
	}

	void Array::push_back(const_reference value) {
		data.push_back(value);
	}

	void Array::pop_back() {
		data.pop_back();
	}

	void Array::resize(size_type count, const_reference value) {
		data.resize(count, value);
	}

	void Array::swap(Array &other) {
		data.swap(other.data);
	}

	std::ostream &operator<<(std::ostream &output, const Array &a) {
		if (a.empty()) {
			output << Structural::BEGIN_ARRAY << Structural::END_ARRAY;

		} else {
			output << Structural::BEGIN_ARRAY << std::endl;
			OutputFilter<Indenter> indent(output.rdbuf());
			output.rdbuf(&indent);

			for (Array::const_iterator i = a.begin(); i != a.end(); ++i) {
				if (i != a.begin()) {
					output << Structural::VALUE_SEPARATOR << std::endl;
				}

				output << *i;
			}

			output.rdbuf(indent.getDestination());

			output << std::endl << Structural::END_ARRAY;
		}

		return output;
	}
}
