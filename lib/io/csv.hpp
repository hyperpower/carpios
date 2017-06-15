#ifndef _CSV_HPP_
#define _CSV_HPP_

#include "io_define.hpp"

#include <stdexcept>
#include <string>
#include <vector>
#include <list>
#include <sstream>

namespace carpio {

namespace csv {

class Error: public std::runtime_error {

public:
	Error(const std::string &msg) :
		std::runtime_error(std::string("CSVparser : ").append(msg)) {
	}
};

class Row {
public:
	typedef std::string Str;
	typedef std::fstream fst;
public:
	Row(const std::vector<Str> &);
	~Row(void);

public:
	uInt size(void) const;
	void push(const Str &);
	bool set(const Str &, const Str &);

private:
	const std::vector<Str> _header;
	std::vector<Str> _values;

public:

	template<typename T>
	const T get_value(uInt pos) const {
		if (pos < _values.size()) {
			T res;
			std::stringstream ss;
			ss << _values[pos];
			ss >> res;
			return res;
		}
		throw Error("can't return this value (doesn't exist)");
	}
	const Str operator[](uInt) const;
	const Str operator[](const Str &valueName) const;
	friend std::ostream& operator<<(std::ostream& os, const Row &row);
	friend std::ofstream& operator<<(std::ofstream& os, const Row &row);
};

enum DataType {
	eFILE = 0, ePURE = 1
};

class Parser {
public:
	typedef std::string Str;
	typedef std::fstream fst;
public:
	/*
	 * Str input can be the file name or pure data
	 * if input is a file
	 *     data is filename
	 *     type = eFILE
	 * if input is data
	 *     data is data
	 *     type = ePURE
	 */
	Parser(const Str &, const DataType &type = eFILE, char sep = ',');
	~Parser(void);

public:
	Row &get_row(uInt row) const;
	uInt row_count(void) const;
	uInt column_count(void) const;
	std::vector<Str> get_header(void) const;
	const Str get_header_element(uInt pos) const;
	const Str &get_file_name(void) const;

public:
	bool delete_row(uInt row);
	bool add_row(uInt pos, const std::vector<Str> &);
	void sync(void) const;

protected:
	void parse_header(void);
	void parse_content(void);

private:
	Str _file;
	const DataType _type;
	const char _sep;
	std::vector<Str> _originalFile;
	std::vector<Str> _header;
	std::vector<Row *> _content;

public:
	Row &operator[](uInt row) const;
};


}

}

#endif
