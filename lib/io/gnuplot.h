#ifndef _GNUPLOT_H_
#define _GNUPLOT_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>              // for std::ostringstream
#include <stdexcept>
#include <cstdio>
#include <stdio.h>
#include <cstdlib>              // for getenv()
#include <list>                 // for std::list
#include <stdlib.h>
#include <type_define.hpp>
#include "io_define.hpp"
#include <memory>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
//defined for 32 and 64-bit environments
#include <io.h>                // for _access(), _mktemp()
#define GP_MAX_TMP_FILES  27   // 27 temporary files it's Microsoft restriction
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
//all UNIX-like OSs (Linux, *BSD, MacOSX, Solaris, ...)
#include <unistd.h>            // for access(), mkstemp()
#define GP_MAX_TMP_FILES  64
#else
#error unsupported or unknown operating system
#endif

namespace carpio {

//
class Gnuplot_actor {
protected:
	std::string _scmd;  //style comand
	std::string _pcmd;  //plot command
	std::list<std::string> _data;
public:
	friend class Gnuplot;
	/*
	 *  Constructor
	 */
	Gnuplot_actor() :
			_scmd(""), _pcmd(""), _data() {
	}
	Gnuplot_actor(const std::string& pcmd, const std::list<std::string>& data) :
			_scmd(""), _pcmd(pcmd), _data(data) {
	}
	Gnuplot_actor(const std::string& pcmd, const std::string& scmd,
			const std::list<std::string>& data) :
			_scmd(scmd), _pcmd(pcmd), _data(data) {
	}
	Gnuplot_actor(const Gnuplot_actor& ga) :
			_scmd(ga._scmd), _pcmd(ga._pcmd), _data(ga._data) {
	}
	/*
	 *  is empty
	 */
	bool empty() const {
		if (_pcmd == "") {
			return true;
		} else {
			return false;
		}
	}
	bool empty_style() const {
		if (_scmd == "") {
			return true;
		} else {
			return false;
		}
	}
	std::string& command() {
		return _pcmd;
	}
	const std::string& command() const {
		return _pcmd;
	}
	std::string& style() {
		return _scmd;
	}
	const std::string& style() const {
		return _scmd;
	}
	std::list<std::string>& data() {
		return _data;
	}
	const std::list<std::string>& data() const {
		return _data;
	}
	void clear() {
		_pcmd = "";
		_scmd = "";
		_data.clear();
	}
	void show_command() const {
		std::cout << "Actor command :" << _pcmd << "\n";
	}
	void show_data() const {
		for (std::list<std::string>::const_iterator ci = _data.begin();
				ci != _data.end(); ++ci) {
			std::cout << (*ci) << std::endl;
		}
	}

};




class Gnuplot {
public:
	typedef std::shared_ptr<Gnuplot_actor> spActor;
	typedef std::list<std::shared_ptr<Gnuplot_actor> > list_spActor;
protected:
	/*
	 * \brief pointer to the stream that can be used to write to the pipe
	 */
	FILE *gnucmd;
	///\brief name of executed GNUPlot file
	std::string m_sGNUPlotFileName;
	///\brief gnuplot path
	std::string m_sGNUPlotPath;
	///\brief standart terminal, used by showonscreen
	std::string terminal_std;

	list_spActor _actors; // if plot() has external parameter, inner actors will be disabled

	bool _valid;
	/*
	 * Opens up a gnuplot session, ready to receive commands
	 */
	void _init();
	/*
	 * Find out if a command lives in m_sGNUPlotPath or in PATH
	 */
	bool _get_program_path();

	/*
	 * Check if file exists
	 */
	bool _file_exists(const std::string &filename, int mode);
public:
	Gnuplot();
	~Gnuplot();

	/*
	 * Sends a command to an active gnuplot session
	 */
	Gnuplot& cmd(const std::string &cmdstr);

	/// turns grid on/off
	inline Gnuplot& set_grid() {
		cmd("set grid");
		return *this;
	}

	/// grid is not set by default
	inline Gnuplot& unset_grid() {
		cmd("unset grid");
		return *this;
	}
	inline Gnuplot& set_view(int rot_x, int rot_z, double scale,
			double scale_z) {
		std::ostringstream sst;
		sst << "set view " << rot_x << ", " << rot_z << ", " << scale << ", "
				<< scale_z;
		cmd(sst.str());
		return *this;
	}
	// -----------------------------------------------
	/// set the mulitplot mode
	///
	/// \param ---
	///
	/// \return <-- reference to the gnuplot object
	// -----------------------------------------------
	inline Gnuplot& set_multiplot() {
		cmd("set multiplot");
		return *this;
	}

	// -----------------------------------------------
	/// unsets the mulitplot mode
	///
	/// \param ---
	///
	/// \return <-- reference to the gnuplot object
	// -----------------------------------------------
	inline Gnuplot& unset_multiplot() {
		cmd("unset multiplot");
		return *this;
	}

	// -----------------------------------------------------------------------
	/// \brief sets and clears the title of a gnuplot session
	///
	/// \param title --> the title of the plot [optional, default == ""]
	///
	/// \return <-- reference to the gnuplot object
	// -----------------------------------------------------------------------
	inline Gnuplot& set_title(const std::string &title = "") {
		std::string cmdstr;
		cmdstr = "set title \"";
		cmdstr += title;
		cmdstr += "\"";
		cmd(cmdstr);
		return *this;
	}

	//----------------------------------------------------------------------------------
	///\brief Clears the title of a gnuplot session
	/// The title is not set by default.
	///
	/// \param ---
	///
	/// \return <-- reference to the gnuplot object
	// ---------------------------------------------------------------------------------
	inline Gnuplot& unset_title() {
		this->set_title();
		return *this;
	}

	/// set x axis label
	Gnuplot& set_xlabel(const std::string &label = "x");
	/// set y axis label
	Gnuplot& set_ylabel(const std::string &label = "y");
	/// set z axis label
	Gnuplot& set_zlabel(const std::string &label = "z");

	Gnuplot& set_equal_ratio();

	Gnuplot& set_label(int, const std::string &, const double&, const double&,
			const std::string &append = "");

	//------------------------------------------------------------------------------
	//
	Gnuplot& set_terminal_pdf(const std::string& filename, double x = 400,
			double y = 300, const std::string& font = "Helvetica",
			int fontsize = 12);

	Gnuplot& set_terminal_png(const std::string& filename, double x = 400,
			double y = 300, const std::string& font = "Helvetica",
			int fontsize = 12);

	Gnuplot& set_terminal_jpeg(const std::string& filename, double x = 400,
				double y = 300, const std::string& font = "Helvetica",
				int fontsize = 12);

	// set
	Gnuplot& set(const std::string& str);

	Gnuplot& set_palette_blue_red();
	// set range
	// set the xrange
	Gnuplot& set_xrange(const double iFrom, const double iTo);
	Gnuplot& set_xrange_reverse(const double iFrom, const double iTo);
	Gnuplot& set_yrange(const double iFrom, const double iTo);
	Gnuplot& set_yrange_reverse(const double iFrom, const double iTo);
	Gnuplot& set_zrange(const double iFrom, const double iTo);
	Gnuplot& set_zrange_reverse(const double iFrom, const double iTo);
	Gnuplot& set_cbrange(const double iFrom, const double iTo);

	/// turns on/off log scaling for the specified xaxis (logscale is set by default)
	Gnuplot& set_xlogscale(const double base = 10);
	/// turns on/off log scaling for the specified yaxis (logscale is set by default)
	Gnuplot& set_ylogscale(const double base = 10);
	/// turns on/off log scaling for the specified zaxis (logscale is set by default)
	Gnuplot& set_zlogscale(const double base = 10);
	/*
	 *  plot
	 */
	template<typename CONTAINER>
	Gnuplot& plot_1(     //
			const CONTAINER& x,  //
			const std::string &str = "") {
		std::ostringstream ss;
		ss << "plot \"-\" using 1 " << str << "\n";
		cmd(ss.str());
		ss.str("");
		for (typename CONTAINER::const_iterator it = x.begin(); it != x.end();
				++it) {
			ss << (*it) << "\n";
			cmd(ss.str());
			ss.str("");
		}
		cmd("e\n");
		return *this;
	}
	template<typename X, typename Y>
	Gnuplot& plot_2(   //  type has [] and size()
			const X& x, //
			const Y& y, //
			const std::string &str = "") {  //
		if (x.size() != y.size()) {
			std::cerr << " >Warning! The containers' size are not equal. \n";
			std::cerr << " >Warning! x =" << x.size() << " y =" << y.size
					<< " \n";
		}
		// inline data
		std::ostringstream sst;
		//
		sst << "plot \"-\" using 1:2 " << str;
		cmd(sst.str());
		sst.str("");
		typename X::const_iterator iterx = x.begin();
		typename Y::const_iterator itery = y.begin();
		if (x.size() >= y.size()) {
			for (; itery != y.end();) {
				sst << (*iterx) << " " << (*itery);
				sst << "\n";
				cmd(sst.str());
				sst.str("");
				iterx++;
				itery++;
			}
		} else {
			for (; iterx != x.end();) {
				sst << (*iterx) << " " << (*itery);
				sst << "\n";
				cmd(sst.str());
				sst.str("");
				iterx++;
				itery++;
			}
		}
		cmd("e\n");
		return *this;
	}

	Gnuplot& plot(const Gnuplot_actor& actor) {
		if (actor.empty()) {
			std::cerr << " >Warning! The Gnuplot actor is empty! \n";
			return *this;
		}
		// inline data
		std::ostringstream sst;
		//
		sst << "plot \"-\" " << actor._pcmd;
		cmd(sst.str());
		sst.str("");
		cmd("\n");
		for (std::list<std::string>::const_iterator iter = actor._data.begin();
				iter != actor._data.end(); ++iter) {
			sst << (*iter);
			cmd(sst.str());
			sst.str("");
			cmd("\n");
		}
		cmd("e\n");
		return *this;
	}

	Gnuplot& plot(const std::list<std::shared_ptr<Gnuplot_actor> >& lga) {
		if (lga.empty()) {
			std::cerr << " >Warning! The Gnuplot actor is empty! \n";
			return *this;
		}
		std::ostringstream ss;
		ss << "plot ";
		for (std::list<std::shared_ptr<Gnuplot_actor> >::const_iterator iter =
				lga.begin(); iter != lga.end(); ++iter) {
			const std::shared_ptr<Gnuplot_actor> spa = (*iter);
			if (spa->empty_style()) {
				ss << "\"-\" " << spa->command() << "with lines lw 1";
			} else {
				ss << "\"-\" " << spa->command() << spa->style();
			}

			if (lga.size() >= 2 && (iter != (--lga.end()))) {
				ss << ",\\\n";
			}
		}
		cmd(ss.str() + "\n");
		ss.str("");
		for (std::list<std::shared_ptr<Gnuplot_actor> >::const_iterator iter =
				lga.begin(); iter != lga.end(); ++iter) {
			const std::shared_ptr<Gnuplot_actor> spa = (*iter);
			output_inline_data((*spa));
		}
		return *this;
	}

	Gnuplot& plot(const std::list<Gnuplot_actor>& lga) {
		if (lga.empty()) {
			std::cerr << " >Warning! The Gnuplot actor is empty! \n";
			return *this;
		}
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
		cmd(ss.str() + "\n");
		ss.str("");
		for (std::list<Gnuplot_actor>::const_iterator iter = lga.begin();
				iter != lga.end(); ++iter) {
			output_inline_data((*iter));
		}
		return *this;
	}

	Gnuplot& splot(const std::list<Gnuplot_actor>& lga) {
		if (lga.empty()) {
			std::cerr << " >Warning! The Gnuplot actor is empty! \n";
			return *this;
		}
		std::ostringstream ss;
		ss << "splot ";
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
		cmd(ss.str() + "\n");
		ss.str("");
		for (std::list<Gnuplot_actor>::const_iterator iter = lga.begin();
				iter != lga.end(); ++iter) {
			output_inline_data((*iter));
		}
		return *this;
	}
	Gnuplot& splot(const Gnuplot_actor& actor) {
		if (actor.empty()) {
			std::cerr << " >Warning! The Gnuplot actor is empty! \n";
			return *this;
		}
		// inline data
		std::ostringstream sst;
		//
		sst << "plot \"-\" " << actor._pcmd;
		cmd(sst.str());
		sst.str("");
		cmd("\n");
		for (std::list<std::string>::const_iterator iter = actor._data.begin();
				iter != actor._data.end(); ++iter) {
			sst << (*iter);
			cmd(sst.str());
			sst.str("");
			cmd("\n");
		}
		cmd("e\n");
		return *this;
	}

	Gnuplot& output_inline_data(const Gnuplot_actor& actor) {
		std::ostringstream sst;
		for (std::list<std::string>::const_iterator iter = actor._data.begin();
				iter != actor._data.end(); ++iter) {
			sst << (*iter);
			cmd(sst.str());
			sst.str("");
		}
		cmd("e\n");
		return *this;
	}

	// inner actors
	Gnuplot& add(spActor actor){
		this->_actors.push_back(actor);
		return *this;
	}
	Gnuplot& clear(){
		this->_actors.clear();
		return *this;
	}


	Gnuplot& plot(){
		this->plot(this->_actors);
		this->clear();
		return *this;
	}

	Gnuplot& save_cmd(const std::string& filename){
		this->cmd("save \'" + filename + "\'");
		return *this;
	}



}
;

//------------------------------------------------------------------------------
//
// A string tokenizer taken from http://www.sunsite.ualberta.ca/Documentation/
// /Gnu/libstdc++-2.90.8/html/21_strings/stringtok_std_h.txt
//
template<typename Container>
void stringtok(Container &container, std::string const &in,
		const char * const delimiters = " \t\n") {
	const std::string::size_type len = in.length();
	std::string::size_type i = 0;

	while (i < len) {
		// eat leading whitespace
		i = in.find_first_not_of(delimiters, i);

		if (i == std::string::npos)
			return;   // nothing left but white space

		// find the end of the token
		std::string::size_type j = in.find_first_of(delimiters, i);

		// push token
		if (j == std::string::npos) {
			container.push_back(in.substr(i));
			return;
		} else
			container.push_back(in.substr(i, j - i));
		// set up for next loop
		i = j + 1;
	}
	return;
}
int GnuplotShow(const std::list<Gnuplot_actor>& lga);
int GnuplotShow(Gnuplot&, const std::list<Gnuplot_actor>& lga);

namespace GnuplotActor {

typedef Float Cvt;
typedef Float Vt;

typedef std::shared_ptr<Gnuplot_actor> spActor;
typedef std::list<spActor> list_spActor;



}



}

#endif
