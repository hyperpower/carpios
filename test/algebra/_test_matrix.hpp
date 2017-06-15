
#ifndef __TEST_MATRIX_HPP_
#define __TEST_MATRIX_HPP_

#include "gtest/gtest.h"
#include "../../lib/algebra/matrix.hpp"

#include <math.h>
#include <string>
#include <memory>

namespace carpio {

TEST(Matrix, construct) {
	typedef MatrixV<double> Mat;
	std::cout<<"Build a matrix \n";
	Mat mat(10, 10);
	std::cout<<" i len = "<< mat.size_i() << "\n";
	std::cout<<" j len = "<< mat.size_j() << "\n";
	Mat::size_type il = mat.size_i();
	Mat::size_type jl = mat.size_j();
	EXPECT_TRUE(il==10);
	EXPECT_TRUE(jl==10);
}




}

#endif
