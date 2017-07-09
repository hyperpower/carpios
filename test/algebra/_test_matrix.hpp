#ifndef __TEST_MATRIX_HPP_
#define __TEST_MATRIX_HPP_

#include "gtest/gtest.h"
#include "../../lib/algebra/matrix.hpp"
#include "../../lib/utility/clock.h"
#ifdef OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <string>
#include <memory>

namespace carpio {

typedef MatrixV<double> Mat;

TEST(Matrix, construct) {

	std::cout << "Build a matrix \n";
	Mat mat(10, 10);
	std::cout << " i len = " << mat.size_i() << "\n";
	std::cout << " j len = " << mat.size_j() << "\n";
	Mat::size_type il = mat.size_i();
	Mat::size_type jl = mat.size_j();
	EXPECT_TRUE(il == 10);
	EXPECT_TRUE(jl == 10);
}

Mat One(St x, St y) {
	Mat res(x, y);
	res.assign(1);
	return std::move(res);
}
Mat Two(St x, St y) {
	Mat res(x, y);
	res.assign(2);
	return std::move(res);
}

double test_mat_add(int n) {
	int nx = n;
	int ny = n;
	Mat a = One(nx, ny);
	Mat b = Two(nx, ny);
	tick_t t = Clock::Tick();
	Mat c = a + b;
	return Clock::MillisecondsSinceD(t);
}

TEST(Matrix, add) {
	int numt = 1;
#ifdef OPENMP
	fmt::print("Open MP add matrix test\n");
	for (numt = 1; numt < 5; numt++) {
		omp_set_num_threads(numt);
#endif
		for (int i = 0; i < 3; i++) {
			int n = 10;
			n *= std::pow(10, i);
			double t = test_mat_add(n);
			fmt::print("Threads : {:5d} n : {:5d} time :{:8.5f}ms\n", numt, n, t);
		}
#ifdef OPENMP
	}
#endif
	//c.show();

}

}

#endif
