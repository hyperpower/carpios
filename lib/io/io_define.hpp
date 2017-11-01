#ifndef IO_DEFINE_H_
#define IO_DEFINE_H_

#include <type_define.hpp>
#include <iostream>
#include <fstream>
#include <unistd.h>

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <list>
#include <map>

namespace carpio {
template<class V1>
std::string ToString(V1 a) {
	std::ostringstream sst;
	sst << a;
	return sst.str();
}
/*
 *  out put a and be to string,
 *  class V1 and v2 must overload operator<<
 *  sep is the separator
 *  For example:
 *  a = 1.3
 *  b = 1.4
 *  sep = " "
 *  return  1.3 1.4
 *
 */
template<class V1, class V2>
std::string ToString(V1 a, V2 b, const std::string sep) {
	std::ostringstream sst;
	sst << a << sep << b;
	return sst.str();
}

template<class V1, class V2, class V3>
std::string ToString(V1 a, V2 b, V3 c, const std::string sep) {
	std::ostringstream sst;
	sst << a << sep << b << sep << c;
	return sst.str();
}

template<class V1, class V2, class V3, class V4, class V5, class V6, class V7>
std::string ToString(V1 a, V2 b, V3 c, V4 d, V5 e, V6 f, V7 g,
		const std::string sep) {
	std::ostringstream sst;
	sst << a << sep; //1
	sst << b << sep; //2
	sst << c << sep; //3
	sst << d << sep; //4
	sst << e << sep; //5
	sst << f << sep; //6
	sst << g;        //7
	return sst.str();
}

inline bool file_access_check( //
		const std::string &filename, //
		int mode //
		) {
	if (mode < 0 || mode > 7) {
		std::cerr << " >! Input mode is wrong  =" << mode
				<< ", it should be from 0 to 7" << std::endl;
		return false;
	}

	//  int _access(const char *path, int mode);
	//  returns 0 if the file has the given mode,
	//  it returns -1 if the named file does not exist or is not accessible in
	//  the given mode
	// mode = 0 (F_OK) (default): checks file for existence only
	// mode = 1 (X_OK): execution permission
	// mode = 2 (W_OK): write permission
	// mode = 4 (R_OK): read permission
	// mode = 6       : read and write permission
	// mode = 7       : read, write and execution permission
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
	if (_access(filename.c_str(), mode) == 0)
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
	if (access(filename.c_str(), mode) == 0)
			#endif
			{
		return true;
	} else {
		return false;
	}
}

class TextFile {
public:
	typedef std::string str;
	typedef std::list<str> lines;
	typedef std::fstream fst;
	typedef std::map<std::string, std::string> dict;
protected:
	/// a text file is consisted from three parts
	///  1. filename
	///  2. a dictionary: illustrate how to read this text file
	///  3. contents

	str _filename;
	lines _content;
	dict _config;
public:
	TextFile() :
			_filename(""), _content() {
	}
	TextFile(const str& filename) :
			_filename(filename), _content() {
	}
	TextFile(const str& filename, const lines& content) {
		_filename = filename;
		_content = content;
	}

	void add_line(const str& line) {
		return _content.push_back(line);
	}

	void _open_read(fst& ins) {
		ins.open(this->_filename.c_str(), std::ifstream::in);
		if (!ins.is_open()) {
			std::cerr << "!> Open file error! " << this->_filename.c_str()
					<< " \n";
			exit(-1);
		}
	}

	void _open_write(fst& outs) {
		outs.open(this->_filename.c_str(), std::fstream::out);
		if (!outs.is_open()) {
			std::cerr << "!> Open file error! " << this->_filename.c_str()
					<< " \n";
			exit(-1);
		}
	}

	void read() {
		fst ins;
		this->_open_read(ins);
		ins.seekg(0, std::ios::beg);
		while (!ins.eof()) {
			str sline;
			getline(ins, sline, '\n');
			this->_content.push_back(sline);
		}
	}

	lines& content(){
		return this->_content;
	}

	void parse_config() {
		for (auto& line : _content) {
			std::vector<std::string> tokens;
			Tokenize(line, tokens);
			int count = 0;
			for (auto& str : tokens) {
				if (count == 0 && str == "##") {
					// this line is a dict line
					this->_config[tokens[1]] = tokens[3];
					break;
				}
			}
		}
	}

	void show_config() {
		for (auto& str : _config) {
			std::cout << str.first << " : " << str.second << "\n";
		}
	}

	std::string get_config(const std::string& key) const {
		typename dict::const_iterator iter = this->_config.find(key);
		if (iter == this->_config.end()) {
			ASSERT_MSG(false, "Not found key");
			return "";
		} else {
			return iter->second;
		}
	}

	void write() {
		fst outs;
		this->_open_write(outs);
		outs.seekg(0, std::ios::beg);
		for (lines::iterator iter = _content.begin(); iter != _content.end();
				++iter) {
			outs << (*iter) << "\n";
		}
	}



	template<class ContainerT>
	static void Tokenize(const std::string& str, ContainerT& tokens,
			const std::string& delimiters = " ", bool trimEmpty = true)
			{
		std::string::size_type pos, lastPos = 0, length = str.length();

		using value_type = typename ContainerT::value_type;
		using size_type = typename ContainerT::size_type;

		while (lastPos < length + 1)
		{
			pos = str.find_first_of(delimiters, lastPos);
			if (pos == std::string::npos) {
				pos = length;
			}

			if (pos != lastPos || !trimEmpty)
				tokens.push_back(value_type(str.data() + lastPos,
						(size_type) pos - lastPos));

			lastPos = pos + 1;
		}
	}

};

}

#endif /* IO_H_ */
