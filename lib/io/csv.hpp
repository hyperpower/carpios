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
	Row(const std::vector<Str> &header) :
			_header(header) {

	}
	~Row(void) {
	}
	uInt size(void) const {
		return _values.size();
	}
	void push(const Str &value) {
		_values.push_back(value);
	}
	bool set(const std::string &key, const std::string &value) {
		std::vector<std::string>::const_iterator it;
		int pos = 0;

		for (it = _header.begin(); it != _header.end(); it++) {
			if (key == *it) {
				_values[pos] = value;
				return true;
			}
			pos++;
		}
		return false;
	}
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
	const Str operator[](uInt valuePosition) const {
		if (valuePosition < _values.size()) {
			return _values[valuePosition];
		}
		throw Error("can't return this value (doesn't exist)");
	}
	const Str operator[](const Str &valueName) const {
		std::vector<std::string>::const_iterator it;
		int pos = 0;

		for (it = _header.begin(); it != _header.end(); it++) {
			if (valueName == *it)
				return _values[pos];
			pos++;
		}
		throw Error("can't return this value (doesn't exist)");
	}
	friend inline std::ostream& operator<<(std::ostream& os, const Row &row);
	friend inline std::ofstream& operator<<(std::ofstream& os, const Row &row);
};

inline std::ostream &operator<<(std::ostream &os, const Row &row) {
	for (uInt i = 0; i != row._values.size(); i++) {
		os << row._values[i] << " | ";
	}
	return os;
}

inline std::ofstream &operator<<(std::ofstream &os, const Row &row) {
	for (uInt i = 0; i != row._values.size(); i++) {
		os << row._values[i];
		if (i < row._values.size() - 1)
			os << ",";
	}
	return os;
}

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
	Parser(const Str &data, const DataType &type = eFILE, char sep = ',') :
			_type(type), _sep(sep) {
		std::string line;
		if (type == eFILE) {
			_file = data;
			std::ifstream ifile(_file.c_str());
			if (ifile.is_open()) {
				while (ifile.good()) {
					getline(ifile, line);
					if (line != "")
						_originalFile.push_back(line);
				}
				ifile.close();

				if (_originalFile.size() == 0)
					throw Error(std::string("No Data in ").append(_file));

				parse_header();
				parse_content();
			} else
				throw Error(std::string("Failed to open ").append(_file));
		} else {
			std::istringstream stream(data);
			while (std::getline(stream, line))
				if (line != "")
					_originalFile.push_back(line);
			if (_originalFile.size() == 0)
				throw Error(std::string("No Data in pure content"));

			parse_header();
			parse_content();
		}
	}
	~Parser(void) {
		std::vector<Row *>::iterator it;

		for (it = _content.begin(); it != _content.end(); it++)
			delete *it;
	}

public:
	Row &get_row(uInt rowPosition) const {
		if (rowPosition < _content.size()) {
			return *(_content[rowPosition]);
		}
		throw Error("can't return this row (doesn't exist)");
	}
	uInt size_row(void) const {
		return _content.size();
	}
	uInt size_column(void) const {
		return _header.size();
	}
	std::vector<Str> get_header(void) const {
		return _header;
	}
	const Str get_header_element(uInt pos) const {
		if (pos >= _header.size())
			throw Error("can't return this header (doesn't exist)");
		return _header[pos];
	}
	const Str &get_file_name(void) const {
		return _file;
	}

public:
	bool delete_row(uInt row) {
		if (row < _content.size()) {
			delete *(_content.begin() + row);
			_content.erase(_content.begin() + row);
			return true;
		}
		return false;
	}
	bool add_row(uInt pos, const std::vector<Str> &r) {
		Row *row = new Row(_header);

		for (auto it = r.begin(); it != r.end(); it++)
			row->push(*it);

		if (pos <= _content.size()) {
			_content.insert(_content.begin() + pos, row);
			return true;
		}
		return false;
	}
	void sync(void) const {
		if (_type == DataType::eFILE) {
			std::ofstream f;
			f.open(_file, std::ios::out | std::ios::trunc);

			// header
			uInt i = 0;
			for (auto it = _header.begin(); it != _header.end(); it++) {
				f << *it;
				if (i < _header.size() - 1)
					f << ",";
				else
					f << std::endl;
				i++;
			}

			for (auto it = _content.begin(); it != _content.end(); it++)
				f << **it << std::endl;
			f.close();
		}
	}

protected:
	void parse_header(void) {
		std::stringstream ss(_originalFile[0]);
		std::string item;

		while (std::getline(ss, item, _sep))
			_header.push_back(item);
	}

	void parse_content(void) {
		std::vector<std::string>::iterator it;

		it = _originalFile.begin();
		it++; // skip header

		for (; it != _originalFile.end(); it++) {
			bool quoted = false;
			int tokenStart = 0;
			uInt i = 0;

			Row *row = new Row(_header);

			for (; i != it->length(); i++) {
				if (it->at(i) == '"')
					quoted = ((quoted) ? (false) : (true));
				else if (it->at(i) == ',' && !quoted) {
					row->push(it->substr(tokenStart, i - tokenStart));
					tokenStart = i + 1;
				}
			}

			//end
			row->push(it->substr(tokenStart, it->length() - tokenStart));

			// if value(s) missing
			if (row->size() != _header.size())
				throw Error("corrupted data !");
			_content.push_back(row);
		}
	}

private:
	Str _file;
	const DataType _type;
	const char _sep;
	std::vector<Str> _originalFile;
	std::vector<Str> _header;
	std::vector<Row *> _content;
public:
	Row &operator[](uInt rowPosition) const{
		return Parser::get_row(rowPosition);
	}
};

}

}

#endif
