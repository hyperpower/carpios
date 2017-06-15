
#ifndef __TEST_MATRIX_SMALL_HPP_
#define __TEST_MATRIX_SMALL_HPP_

#include "gtest/gtest.h"
#include "../../lib/algebra/matrix_small.hpp"

#include <math.h>
#include <string>
#include <memory>

namespace carpio {

TEST(MatrixS, construct) {
	typedef MatrixS<double, 3, 3> Mat;
	std::cout<<"Build a matrix small \n";
	Mat mat;
	std::cout<<" i len = "<< mat.size_i() << "\n";
	std::cout<<" j len = "<< mat.size_j() << "\n";
	Mat::size_type il = mat.size_i();
	Mat::size_type jl = mat.size_j();
	EXPECT_TRUE(il==3);
	EXPECT_TRUE(jl==3);
}




}

#endif
