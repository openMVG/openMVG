#include <JsonBox/Object.h>

#include <JsonBox/Value.h>
#include <JsonBox/Grammar.h>
#include <JsonBox/OutputFilter.h>
#include <JsonBox/Indenter.h>

namespace JsonBox {
	Object::Object(const key_compare &comp, const allocator_type &alloc) : data(comp, alloc) {
	}

	Object::Object(const Object &other) : data(other.data) {
	}

	Object &Object::operator=(const Object &other) {
		data = other.data;
		return *this;
	}
	
	bool Object::operator==(const Object &rhs) const {
		return data == rhs.data;
	}
	
	bool Object::operator!=(const Object &rhs) const {
		return data != rhs.data;
	}
	
	bool Object::operator<(const Object &rhs) const {
		return data < rhs.data;
	}
	
	bool Object::operator<=(const Object &rhs) const {
		return data <= rhs.data;
	}
	
	bool Object::operator>(const Object &rhs) const {
		return data > rhs.data;
	}
	
	bool Object::operator>=(const Object &rhs) const {
		return data >= rhs.data;
	}

	Object::allocator_type Object::get_allocator() const {
		return data.get_allocator();
	}

	Object::mapped_type &Object::operator[](const key_type &key) {
		return data[key];
	}

	Object::iterator Object::begin() {
		return data.begin();
	}

	Object::const_iterator Object::begin() const {
		return data.begin();
	}

	Object::iterator Object::end() {
		return data.end();
	}

	Object::const_iterator Object::end() const {
		return data.end();
	}

	Object::reverse_iterator Object::rbegin() {
		return data.rbegin();
	}

	Object::const_reverse_iterator Object::rbegin() const {
		return data.rbegin();
	}

	Object::reverse_iterator Object::rend() {
		return data.rend();
	}

	Object::const_reverse_iterator Object::rend() const {
		return data.rend();
	}

	bool Object::empty() const {
		return data.empty();
	}

	Object::size_type Object::size() const {
		return data.size();
	}

	Object::size_type Object::max_size() const {
		return data.max_size();
	}

	void Object::clear() {
		data.clear();
	}

	std::pair<Object::iterator, bool> Object::insert(const_reference value) {
		return data.insert(value);
	}

	Object::iterator Object::insert(iterator hint, const_reference value) {
		return data.insert(hint, value);
	}

	void Object::erase(iterator position) {
		data.erase(position);
	}

	void Object::erase(iterator first, iterator last) {
		data.erase(first, last);
	}

	Object::size_type Object::erase(const key_type &key) {
		return data.erase(key);
	}

	void Object::swap(Object &other) {
		data.swap(other.data);
	}

	Object::size_type Object::count(const key_type &key) const {
		return data.count(key);
	}

	Object::iterator Object::find(const key_type &key) {
		return data.find(key);
	}

	Object::const_iterator Object::find(const key_type &key) const {
		return data.find(key);
	}

	std::pair<Object::iterator, Object::iterator> Object::equal_range(const key_type &key) {
		return data.equal_range(key);
	}

	std::pair<Object::const_iterator, Object::const_iterator> Object::equal_range(const key_type &key) const {
		return data.equal_range(key);
	}

	Object::iterator Object::lower_bound(const key_type &key) {
		return data.lower_bound(key);
	}

	Object::const_iterator Object::lower_bound(const key_type &key) const {
		return data.lower_bound(key);
	}

	Object::iterator Object::upper_bound(const key_type &key) {
		return data.upper_bound(key);
	}

	Object::const_iterator Object::upper_bound(const key_type &key) const {
		return data.upper_bound(key);
	}

	Object::key_compare Object::key_comp() const {
		return data.key_comp();
	}

	std::ostream &operator<<(std::ostream &output, const Object &o) {
		// If the object is empty, we simply write "{}".
		if (o.empty()) {
			output << Structural::BEGIN_OBJECT << Structural::END_OBJECT;

		} else {
			output << Structural::BEGIN_OBJECT << std::endl;
			OutputFilter<Indenter> indent(output.rdbuf());
			output.rdbuf(&indent);

			// For each item in the object.
			for (Object::const_iterator i = o.begin(); i != o.end(); ++i) {
				if (i != o.begin()) {
					output << Structural::VALUE_SEPARATOR << std::endl;
				}

				// We print the name of the attribute and its value.
				output << Structural::BEGIN_END_STRING << Value::escapeMinimumCharacters(i->first) << Structural::BEGIN_END_STRING << Whitespace::SPACE << Structural::NAME_SEPARATOR << Whitespace::SPACE << i->second;
			}

			output.rdbuf(indent.getDestination());

			output << std::endl << Structural::END_OBJECT;
		}

		return output;
	}
}
