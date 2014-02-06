#ifndef JB_ARRAY_H
#define JB_ARRAY_H

#include <iostream>
#include <vector>

#include "Value.h"

namespace JsonBox {
	/**
	 * Represents an array of values in JSON. It's a vector with added methods.
	 * So it can be used the same way as a standard STL vector, but can be more
	 * easily output in a stream.
	 * @see JsonBox::Value
	 */
	class Array {
	public:
		typedef std::vector<Value> container;
		typedef container::value_type value_type;
		typedef container::allocator_type allocator_type;
		typedef container::size_type size_type;
		typedef container::difference_type difference_type;
		typedef container::reference reference;
		typedef container::const_reference const_reference;
		typedef container::pointer pointer;
		typedef container::const_pointer const_pointer;
		typedef container::iterator iterator;
		typedef container::const_iterator const_iterator;
		typedef container::reverse_iterator reverse_iterator;
		typedef container::const_reverse_iterator const_reverse_iterator;

		Array(const allocator_type &alloc = allocator_type());

		explicit Array(size_type count, const_reference value = value_type(), const allocator_type &alloc = allocator_type());

		template <typename InputIterator>
		Array(InputIterator first, InputIterator last, const allocator_type &alloc = allocator_type()) : data(first, last) {
		}

		Array(const Array &other);

		Array &operator=(const Array &other);

		bool operator==(const Array &rhs) const;

		bool operator!=(const Array &rhs) const;

		bool operator<(const Array &rhs) const;

		bool operator<=(const Array &rhs) const;

		bool operator>(const Array &rhs) const;

		bool operator>=(const Array &rhs) const;

		void assign(size_type count, const_reference value);

		template <typename InputIterator>
		void assign(InputIterator first, InputIterator last) {
			data.assign(first, last);
		}

		allocator_type get_allocator() const;

		reference at(size_type pos);

		const_reference at(size_type pos) const;

		reference operator[](size_type pos);

		const_reference operator[](size_type pos) const;

		reference front();

		const_reference front() const;

		reference back();

		const_reference back() const;

		iterator begin();

		const_iterator begin() const;

		iterator end();

		const_iterator end() const;

		reverse_iterator rbegin();

		const_reverse_iterator rbegin() const;

		reverse_iterator rend();

		const_reverse_iterator rend() const;

		bool empty() const;

		size_type size() const;

		size_type max_size() const;

		void reserve(size_type size);

		size_type capacity() const;

		void clear();

		iterator insert(iterator pos, const_reference value);

		void insert(iterator pos, size_type count, const_reference value);

		template <typename InputIterator>
		void insert(iterator pos, InputIterator first, InputIterator last) {
			data.insert(pos, first, last);
		}

		iterator erase(iterator pos);

		iterator erase(iterator first, iterator last);

		void push_back(const_reference value);

		void pop_back();

		void resize(size_type count, const_reference value = value_type());

		void swap(Array &other);
	private:
		container data;
	};

	/**
	 * Output operator overload for the JSON array. Outputs in standard JSON
	 * format. Indents the output and esapes the minimum characters.
	 * @param output Output stream in which to write the array's JSON.
	 * @param a Array to output into the stream.
	 * @return Output stream filled with the JSON code.
	 */
	std::ostream &operator<<(std::ostream &output, const Array &a);
}

#endif
