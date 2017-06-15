#ifndef IO_DEFINE_H_
#define IO_DEFINE_H_

#include <type_define.hpp>
#include <iostream>
#include <fstream>
#include <unistd.h>

#include <string>
#include <fstream>
#include <sstream>
#include <list>

namespace carpio {
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
protected:
	str _filename;
	lines _content;

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

	void add_line(const str& line){
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

	void write() {
		fst outs;
		this->_open_write(outs);
		outs.seekg(0, std::ios::beg);
		for (lines::iterator iter = _content.begin(); iter != _content.end();
				++iter) {
			outs << (*iter) << "\n";
		}
	}

};

}

#endif /* IO_H_ */
