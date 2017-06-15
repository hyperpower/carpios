#ifndef __TEST_CSV_HPP_
#define __TEST_CSV_HPP_

#include "io/gnuplot.h"

#include "gtest/gtest.h"
#include "io/csv.hpp"
#include <math.h>

namespace carpio {

TEST(CSV, test1) {
	std::cout << "CSV Test \n";
	//try {
	csv::Parser file = csv::Parser("./test/input_files/wiki.csv");

	std::cout << file[0][0] << std::endl; // display : 1997

	std::cout << file[0] << std::endl;    // display : 1997 | Ford | E350

	std::cout << file[1]["Model"] << std::endl; // display : Cougar

	std::cout << file.row_count() << std::endl; // display : 2
	std::cout << file.column_count() << std::endl; // display : 3

	std::cout << file.get_header_element(3) << std::endl; // display : Model
	//} catch (csv::Error &e) {
	//	std::cerr << e.what() << std::endl;
	//}
}

TEST(CSV, test2) {
	std::cout << "CSV Test 2\n";
	try {
		std::string data =
				"Year,Make,Model,Description,Price\n1997,Ford,E350,\"ac, abs, moon\",3000.00\n1999,Chevy,\"Venture \"\"Extended Edition\"\"\",\"\",4900.00\n1999,Chevy,\"Venture \"\"Extended Edition, Very Large\"\"\",,5000.00\n1996,Jeep,Grand Cherokee,\"MUST SELL! air, moon roof, loaded\",4799.00";
		csv::Parser file = csv::Parser(data, csv::ePURE);
		std::cout << file[0][0] << std::endl;       // display : 1997
		std::cout << file[0] << std::endl;          // display : 1997 | Ford | E350
		std::cout << file[1]["Model"] << std::endl; // display : Cougar
		std::cout << file.row_count() << std::endl; // display : 2
		std::cout << file.column_count() << std::endl; // display : 3

		std::cout << file.get_header_element(3) << std::endl; // display : Model
	} catch (csv::Error &e) {
		std::cerr << e.what() << std::endl;
	}
}

}

#endif
