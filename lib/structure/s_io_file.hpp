#ifndef _S_IO_FILE_HPP
#define _S_IO_FILE_HPP

#include "s_grid.hpp"
#include "s_define.hpp"
#include "s_vector.hpp"
#include "s_data.hpp"
#include <fstream>

#include "utility/format.h"
#include "io/io_define.hpp"

namespace structure {

template<St DIM>
int Output(const std::string& filename, const Scalar_<DIM>& data) {
	// Open a file
	FILE* fp = std::fopen(filename.c_str(), "w");
	if (!fp) {
		std::perror("File opening failed");
		return 0;
	}
	typedef typename Scalar_<DIM>::Grid Grid;
	Grid& grid = *(data.get_grid());
	// format first line
	typename Grid::Index n = grid.n();
	fmt::print(fp, "# SIZE:{0:d} DIM:{1:d} NX:{2:d} NY:{3:d} NZ:{4:d}\n",
			data.num_cells(), DIM, n(_X_), n(_Y_), n(_Z_));
	//fmt::print(fp, "# DIM:{0:d}", DIM);
	//
	for (typename Grid::Ijk ijk = grid.begin_ijk(); !ijk.is_end(); ++ijk) {
		for (St d = 0; d < 3; ++d) { // always 3 d
			Vt cor = grid.c_(d, ijk.current());
			fmt::print(fp, "{0:10.5f},", cor);
		}
		Vt val = data(ijk.current());
		fmt::print(fp, "{0:10.5f}\n", val);
	}
	std::fclose(fp);
}

template<St DIM>
int Output_Scalar(const std::string& filename, const Scalar_<DIM>& data) {
	// Open a file
	typedef typename Scalar_<DIM>::Grid Grid;
	Grid& grid = *(data.get_grid());
	// format first line
	typename Grid::Index n = grid.n();
	carpio::TextFile txtf(filename);
	txtf.add_line(fmt::format("## Size: {0:d}", data.num_cells()));
	txtf.add_line(fmt::format("## Dim : {0:d}", DIM));
	txtf.add_line(fmt::format("## NX :  {0:d}", n(_X_)));
	if (DIM >= 2) {
		txtf.add_line(fmt::format("## NY :  {0:d}", n(_Y_)));
	}
	if (DIM == 3) {
		txtf.add_line(fmt::format("## NZ :  {0:d}", n(_Z_)));
	}
	//fmt::print(fp, "# DIM:{0:d}", DIM);
	//
	for (typename Grid::Ijk ijk = grid.begin_ijk(); !ijk.is_end(); ++ijk) {
		fmt::MemoryWriter mw;
		for (St d = 0; d < 3; ++d) { // always 3 d
			Vt cor = grid.c_(d, ijk.current());
			mw.write("{0:10.6e}, ", cor);
		}
		Vt val = data(ijk.current());
		mw.write("{0:10.6e}", val);
		txtf.add_line(mw.str());
	}
	txtf.write();
}

inline int Output_PointData(const std::string& filename, Vt val, Vt x, Vt y,
		Vt z) {
	carpio::TextFile txtf(filename);
	txtf.add_line(fmt::format("## Size: {0:d}", 1));
	//
	fmt::MemoryWriter mw;
	mw.write("{0:10.5f},", x);
	mw.write("{0:10.5f},", y);
	mw.write("{0:10.5f},", z);
	mw.write("{0:10.5f}", val);
	txtf.add_line(mw.str());
	txtf.write();
	return 0;
}




}

#endif
