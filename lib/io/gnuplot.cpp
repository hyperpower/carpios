#include "gnuplot.h"

namespace carpio {

Gnuplot::Gnuplot() {
	_init();
}

/*
 * Opens up a gnuplot session, ready to receive commands
 */
void Gnuplot::_init() {
	//------------------------------------------------------------------------------
	//
	// initialize data
	//

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
	std::string this->m_sGNUPlotFileName = "pgnuplot.exe";
	std::string this->m_sGNUPlotPath = "C:/program files/gnuplot/bin/";
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
	this->m_sGNUPlotFileName = "gnuplot";
	this->m_sGNUPlotPath = "/usr/local/bin";
#endif

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
	std::string this->terminal_std = "windows";
#elif ( defined(unix) || defined(__unix) || defined(__unix__) ) && !defined(__APPLE__)
	this->terminal_std = "wxt";
#elif defined(__APPLE__)
	this->terminal_std = "qt";
#endif

	// char * getenv ( const char * name );  get value of environment variable
	// Retrieves a C string containing the value of the environment variable
	// whose name is specified as argument.  If the requested variable is not
	// part of the environment list, the function returns a NULL pointer.
#if ( defined(unix) || defined(__unix) || defined(__unix__) ) && !defined(__APPLE__)
	if (getenv("DISPLAY") == nullptr) {
		this->_valid = false;
		std::cerr << "!> Can't find DISPLAY variable! " << " \n";
	}
#endif

	// if gnuplot not available
	if (!this->_get_program_path()) {
		this->_valid = false;
		std::cerr << "!> Can't find gnuplot! " << " \n";
	}

	//
	// open pipe
	//
	std::string tmp = this->m_sGNUPlotPath + "/" + this->m_sGNUPlotFileName
			+ " -persist";

	// FILE *popen(const char *command, const char *mode);
	// The popen() function shall execute the command specified by the string
	// command, create a pipe between the calling program and the executed
	// command, and return a pointer to a stream that can be used to either read
	// from or write to the pipe.
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
	gnucmd = _popen(tmp.c_str(),"w");
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
	gnucmd = popen(tmp.c_str(), "w");
#endif

	// popen() shall return a pointer to an open stream that can be used to read
	// or write to the pipe.  Otherwise, it shall return a null pointer and may
	// set errno to indicate the error.
	if (!gnucmd) {
		this->_valid = false;
		std::cerr << "!> Couldn't open connection to gnuplot! " << " \n";
	}

	this->_valid = true;

	//set terminal type
	cmd("set output");
	cmd("set terminal " + this->terminal_std + " enhanced font 'Helvetica,12'");

	return;
}

Gnuplot& Gnuplot::set_terminal_pdf(const std::string& filename, double x,
		double y, const std::string& font,
		int fontsize ) {

	this->terminal_std = "pdf";
	std::stringstream sst;
	sst<< "set terminal " << this->terminal_std <<
			" enhanced font '" << font
			<< "," << fontsize << "'"
			<< "size "<< x <<", " <<y;
	cmd(sst.str());
	cmd("set output '"+ filename + "'");
	return *this;
}

Gnuplot& Gnuplot::set_terminal_png(const std::string& filename, double x,
		double y, const std::string& font,
		int fontsize ) {

	this->terminal_std = "pngcairo";
	std::stringstream sst;
	sst<< "set terminal " << this->terminal_std <<
			" enhanced font '" << font
			<< "," << fontsize << "'"
			<< "size "<< x <<", " <<y;
	cmd(sst.str());
	cmd("set output '"+ filename + "'");
	return *this;
}

Gnuplot& Gnuplot::set_terminal_jpeg(const std::string& filename, double x,
		double y, const std::string& font,
		int fontsize ) {

	this->terminal_std = "jpeg";
	std::stringstream sst;
	sst<< "set terminal " << this->terminal_std <<
			" enhanced font '" << font
			<< "," << fontsize << "'"
			<< "size "<< x <<", " <<y;
	cmd(sst.str());
	cmd("set output '"+ filename + "'");
	return *this;
}

/*
 *  Find out if a command lives in m_sGNUPlotPath or in PATH
 */
bool Gnuplot::_get_program_path() {
	//
	// first look in m_sGNUPlotPath for Gnuplot
	//
	std::string tmp = this->m_sGNUPlotPath + "/" + this->m_sGNUPlotFileName;

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
	if ( this->_file_exists(tmp,0) ) // check existence
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
	if (this->_file_exists(tmp, 1)) // check existence and execution permission
#endif
			{
		return true;
	}

	//
	// second look in PATH for Gnuplot
	//
	char *path;
	// Retrieves a C string containing the value of environment variable PATH
	path = getenv("PATH");

	if (path == NULL) {
		std::cerr << "!> PATH is not set! " << " \n";
		return false;
	} else {
		std::list<std::string> ls;
		//split path (one long string) into list ls of strings
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
		stringtok(ls,idx,";");
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
		stringtok(ls, path, ":");
#endif
		// scan list for Gnuplot program files
		for (std::list<std::string>::const_iterator i = ls.begin();
				i != ls.end(); ++i) {
			tmp = (*i) + "/" + this->m_sGNUPlotFileName;
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
			if ( this->_file_exists(tmp,0) ) // check existence
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
			if (this->_file_exists(tmp, 1)) // check existence and execution permission
#endif
					{
				this->m_sGNUPlotPath = *i; // set m_sGNUPlotPath
				return true;
			}
		}

		std::cerr << "Can't find gnuplot neither in PATH nor in \""
				<< this->m_sGNUPlotPath << " \n";

		this->m_sGNUPlotPath = " ";
		return false;
	}
}

bool Gnuplot::_file_exists(const std::string &filename, int mode) {
	if (mode < 0 || mode > 7) {
		throw std::runtime_error(
				"In function \"Gnuplot::file_exists\": mode\
                has to be an integer between 0 and 7");
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
//
// Sends a command to an active gnuplot session
//
Gnuplot& Gnuplot::cmd(const std::string &cmdstr) {
	if (!(this->_valid)) {
		return *this;
	}

	// int fputs ( const char * str, FILE * stream );
	// writes the string str to the stream.
	// The function begins copying from the address specified (str) until it
	// reaches the terminating null character ('\0'). This final
	// null-character is not copied to the stream.
	fputs((cmdstr + "\n").c_str(), gnucmd);

	// int fflush ( FILE * stream );
	// If the given stream was open for writing and the last i/o operation was
	// an output operation, any unwritten data in the output buffer is written
	// to the file.  If the argument is a null pointer, all open files are
	// flushed.  The stream remains open after this call.
	fflush(gnucmd);

	return *this;
}

Gnuplot& Gnuplot::set_equal_ratio() {
	std::ostringstream cmdstr;
	cmdstr << "set size ratio -1";
	cmd(cmdstr.str());

	return *this;

}
/// turns on/off log scaling for the specified xaxis (logscale is not set by default)
Gnuplot& Gnuplot::set_xlogscale(const double base) {
	std::ostringstream cmdstr;

	cmdstr << "set logscale x " << base;
	cmd(cmdstr.str());

	return *this;
}
/// turns on/off log scaling for the specified yaxis (logscale is not set by default)
Gnuplot& Gnuplot::set_ylogscale(const double base) {
	std::ostringstream cmdstr;

	cmdstr << "set logscale y " << base;
	cmd(cmdstr.str());

	return *this;
}
/// turns on/off log scaling for the specified zaxis (logscale is not set by default)
Gnuplot& Gnuplot::set_zlogscale(const double base) {
	std::ostringstream cmdstr;

	cmdstr << "set logscale z " << base;
	cmd(cmdstr.str());

	return *this;
}

Gnuplot& Gnuplot::set_label(int tag, const std::string & label, const double& x,
		const double& y, const std::string& append) {
	std::ostringstream cmdstr;
	cmdstr << "set label " << tag << " \"" << label << "\" at first " << x
			<< ",first " << y << " " << append;
	//std::cout<<cmdstr.str();
	cmd(cmdstr.str());
	return *this;
}

/// set x axis label
Gnuplot& Gnuplot::set_ylabel(const std::string &label) {
	std::ostringstream cmdstr;

	cmdstr << "set xlabel \"" << label << "\"";
	cmd(cmdstr.str());

	return *this;
}
/// set y axis label
Gnuplot& Gnuplot::set_xlabel(const std::string &label) {
	std::ostringstream cmdstr;

	cmdstr << "set ylabel \"" << label << "\"";
	cmd(cmdstr.str());

	return *this;
}
/// set z axis label
Gnuplot& Gnuplot::set_zlabel(const std::string &label) {
	std::ostringstream cmdstr;

	cmdstr << "set zlabel \"" << label << "\"";
	cmd(cmdstr.str());

	return *this;
}

//------------------------------------------------------------------------------
//
// set range
//
// set the xrange
Gnuplot& Gnuplot::set_xrange(const double iFrom, const double iTo) {
	std::ostringstream cmdstr;
	cmdstr << "set xrange[" << iFrom << ":" << iTo << "]";
	cmd(cmdstr.str());

	return *this;
}

Gnuplot& Gnuplot::set_xrange_reverse(const double iFrom, const double iTo) {
	std::ostringstream cmdstr;
	cmdstr << "set xrange[" << iFrom << ":" << iTo << "] reverse";
	cmd(cmdstr.str());

	return *this;
}

Gnuplot& Gnuplot::set_yrange(const double iFrom, const double iTo) {
	std::ostringstream cmdstr;
	cmdstr << "set yrange[" << iFrom << ":" << iTo << "]";
	cmd(cmdstr.str());
	return *this;
}

Gnuplot& Gnuplot::set_yrange_reverse(const double iFrom, const double iTo) {
	std::ostringstream cmdstr;
	cmdstr << "set yrange[" << iFrom << ":" << iTo << "] reverse";
	cmd(cmdstr.str());

	return *this;
}

Gnuplot& Gnuplot::set_zrange(const double iFrom, const double iTo) {
	std::ostringstream cmdstr;
	cmdstr << "set zrange[" << iFrom << ":" << iTo << "]";
	cmd(cmdstr.str());

	return *this;
}

Gnuplot& Gnuplot::set_zrange_reverse(const double iFrom, const double iTo) {
	std::ostringstream cmdstr;
	cmdstr << "set zrange[" << iFrom << ":" << iTo << "] reverse";
	cmd(cmdstr.str());

	return *this;
}

Gnuplot& Gnuplot::set(const std::string& str) {
	std::ostringstream cmdstr;
	cmdstr << "set " << str;
	cmd(cmdstr.str());

	return *this;
}

Gnuplot& Gnuplot::set_cbrange(const double iFrom, const double iTo) {
	std::ostringstream cmdstr;
	cmdstr << "set cbrange[" << iFrom << ":" << iTo << "]";
	cmd(cmdstr.str());

	return *this;
}

Gnuplot& Gnuplot::set_palette_blue_red() {
	std::ostringstream cmdstr;
	cmdstr
			<< "set palette defined(\
	0       0.2314  0.2980  0.7529,\
	0.03125 0.2667  0.3529  0.8000,\
	0.0625  0.3020  0.4078  0.8431,\
	0.09375 0.3412  0.4588  0.8824,\
	0.125   0.3843  0.5098  0.9176,\
	0.15625 0.4235  0.5569  0.9451,\
	0.1875  0.4667  0.6039  0.9686,\
	0.21875 0.5098  0.6471  0.9843,\
	0.25    0.5529  0.6902  0.9961,\
	0.28125 0.5961  0.7255  1.0000,\
	0.3125  0.6392  0.7608  1.0000,\
	0.34375 0.6824  0.7882  0.9922,\
	0.375   0.7216  0.8157  0.9765,\
	0.40625 0.7608  0.8353  0.9569,\
	0.4375  0.8000  0.8510  0.9333,\
	0.46875 0.8353  0.8588  0.9020,\
	0.5     0.8667  0.8667  0.8667,\
	0.53125 0.8980  0.8471  0.8196,\
	0.5625  0.9255  0.8275  0.7725,\
	0.59375 0.9451  0.8000  0.7255,\
	0.625   0.9608  0.7686  0.6784,\
	0.65625 0.9686  0.7333  0.6275,\
	0.6875  0.9686  0.6941  0.5804,\
	0.71875 0.9686  0.6510  0.5294,\
	0.75    0.9569  0.6039  0.4824,\
	0.78125 0.9451  0.5529  0.4353,\
	0.8125  0.9255  0.4980  0.3882,\
	0.84375 0.8980  0.4392  0.3451,\
	0.875   0.8706  0.3765  0.3020,\
	0.90625 0.8353  0.3137  0.2588,\
	0.9375  0.7961  0.2431  0.2196,\
	0.96875 0.7529  0.1569  0.1843,\
	1       0.7059  0.0157  0.1490\
	)";
	cmd(cmdstr.str());

	return *this;
}

//
// Destructor: needed to delete temporary files
//
Gnuplot::~Gnuplot() {
	//remove_tmpfiles();

	// A stream opened by popen() should be closed by pclose()
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
	if (_pclose(gnucmd) == -1) {
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
	if (pclose(gnucmd) == -1) {
#endif
		std::cerr << " >! Problem closing communication to gnuplot \n";
	}
}

int GnuplotShow(const std::list<Gnuplot_actor>& lga) {
	//
	Gnuplot gp;
	gp.set_equal_ratio();
	std::ostringstream ss;
	ss << "plot ";
	for (std::list<Gnuplot_actor>::const_iterator iter = lga.begin();
			iter != lga.end(); ++iter) {
		if (iter->empty_style()) {
			ss << "\"-\" " << iter->command() << "with lines lw 2";
		} else {
			ss << "\"-\" " << iter->command() << iter->style();
		}

		if (lga.size() >= 2 && (iter != (--lga.end()))) {
			ss << ",\\\n";
		}
	}
	gp.cmd(ss.str() + "\n");
	ss.str("");
	for (std::list<Gnuplot_actor>::const_iterator iter = lga.begin();
			iter != lga.end(); ++iter) {
		gp.output_inline_data((*iter));
	}
	return _SUCCESS;
}

int GnuplotShow(Gnuplot& gp, const std::list<Gnuplot_actor>& lga) {
	//
	//gp.set_equal_ratio();
	std::ostringstream ss;
	ss << "plot ";
	for (std::list<Gnuplot_actor>::const_iterator iter = lga.begin();
			iter != lga.end(); ++iter) {
		if (iter->empty_style()) {
			ss << "\"-\" " << iter->command() << "with lines lw 1";
		} else {
			ss << "\"-\" " << iter->command() << iter->style();
		}

		if (lga.size() >= 2 && (iter != (--lga.end()))) {
			ss << ",\\\n";
		}
	}
	gp.cmd(ss.str() + "\n");
	ss.str("");
	for (std::list<Gnuplot_actor>::const_iterator iter = lga.begin();
			iter != lga.end(); ++iter) {
		gp.output_inline_data((*iter));
	}
	return _SUCCESS;
}
}
