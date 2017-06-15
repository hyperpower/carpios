#ifndef _TEST_POLYNOMIAL_H_
#define _TEST_POLYNOMIAL_H_

#include "algebra/algebra.hpp"
#include "algebra/vector_list.hpp"
#include "utility/random.h"
//#include "../calculation/expression.hpp"
#include "gtest/gtest.h"
#include "io/gnuplot_actor.h"

#include <string>
#include "utility/clock.h"

namespace carpio {
TEST(test_polynomial,simple1) {
	typedef Polynomial_<Float, std::string, int> Poly;
	Poly poly;
	Poly::Term t1(1.0, "b", 1);
	Poly::Term t2(4.0, "a", 1);
	Poly::Term t5(4.0, "c", 1);
	Poly::Term t3(5.0, "a", 1);
	Poly::Term t4(6.0, "a", 1);
	std::cout<< "===== test concise =====\n";
	poly.insert(t1);
	poly.insert(t2);
	poly.insert(t3);
	poly.insert(t4);
	poly.insert(t5);
	poly.show();
	std::cout<< " concise ==>\n";
	poly.concise();
	EXPECT_EQ(3, poly.size());
	poly.show();
	std::cout<<" ===== =====\n";
}

TEST(test_polynomial,const1) {
	typedef Polynomial_<Float, std::string, int> Poly;
	Poly poly;
	Poly::Term t1(3.0, "a", 1);
	Poly::Term t2(4.0, "a", 1);
	Poly::Term t3(5.0, "a", 0);
	Poly::Term t4(6.0, "a", 0);
	poly.insert(t1);
	poly.insert(t2);
	poly.insert(t3);
	poly.insert(t4);
	//EXPECT_EQ(3, poly.size());
	poly.show();
}

TEST(DISABLED_test_polynomial,plus_ploy_DISABLED) {
	typedef Polynomial_<Float, std::string, int> Poly;
	Poly poly;
	Poly::Term t1(1.0, "a", 1);
	Poly::Term t2(2.0, "b", 1);
	Poly::Term t3(3.0, "d", 0);
	poly.insert(t1);
	poly.insert(t2);
	poly.insert(t3);
	Poly poly2;
	Poly::Term t12(1.0, "a", 1);
	Poly::Term t22(2.0, "b", 1);
	Poly::Term t32(2.0, "b", 0);
	poly2.insert(t12);
	poly2.insert(t22);
	poly2.insert(t32);

	//EXPECT_EQ(3, poly.size());
	poly.show();
	poly2.show();
	poly.plus(poly2);
	poly.show();
	std::cout << "concise ------ \n" << std::endl;
	poly.concise();
	poly.show();
}

TEST(DISABLED_test_polynomial,plus_ploy2) {
	typedef Polynomial_<Float, int, int> Poly;
	Poly poly;
	Poly::Term t1(1.0, 10, 1);
	Poly::Term t2(2.0, 11, 0);
	poly.insert(t1);
	poly.insert(t2);
	Poly poly2;
	Poly::Term t12(1.0, 13, 1);
	Poly::Term t22(2.0, 14, 0);
	poly2.insert(t12);
	poly2.insert(t22);

	//EXPECT_EQ(3, poly.size());
	poly.show();
	poly2.show();
	poly.plus(poly2);
	poly.show();
	std::cout << "concise ------ \n" << std::endl;
	poly.concise();
	poly.show();
}

TEST(test_polynomial,vector_list_base) {
	Clock t;
	VectorList_<int> vl(10);
	std::cout << " size :" << vl.size() << std::endl;
	std::cout << " cap  :" << vl.capacity() << std::endl;
	//vl[1] = 5;
	vl.push_back(3);
	for (VectorList_<int>::iterator iter = vl.begin(); iter != vl.end();
			++iter) {
		std::cout << *iter << std::endl;
	}

}

Polynomial_<Float, int, int> get_a_ploy() {
	typedef Polynomial_<Float, int, int> Poly;
	Poly res;
	Poly::Term t1(1.0, 10, 1);
	Poly::Term t2(2.0, 11, 0);
	res.insert(t1);
	res.insert(t2);

	return res;
}

TEST(test_polynomial,move_contructor) {
	typedef Polynomial_<Float, int, int> Poly;
	Poly p = get_a_ploy();
	p.show();
}


void poly_add(long num_repeat) {
	for (long n = 0; n < num_repeat; n++) {
		Random::randomSeed();
		// generate the first ploynomial
		// class COE, class TERM, class EXP,
		Polynomial_<Float, std::string, int> poly;
		int num_term = Random::nextInt(1, 5);
		for (int i = 0; i < num_term; i++) {
			poly.insert(Random::nextFloat(-100, 100),
					Random::nextString(1, "abcdefghijkz"),
					Random::nextInt(0, 1));
		}
		Polynomial_<Float, std::string, int> poly2;
		num_term = Random::nextInt(1, 5);
		for (int i = 0; i < num_term; i++) {
			poly2.insert(Random::nextFloat(-100, 100),
					Random::nextString(1, "abcdefghijkz"),
					Random::nextInt(0, 1));
		}
		poly.plus(poly2);
	}
}

TEST(DISABLED_test_polynomial,random_number) {
	//Polynomial2_<std::string, Float> poly("z");
	//Expression2_<Float, Float, 3> exp;
	std::vector<long> a = { 10, 100, 1000, 5000};
	Clock c;
	c.start();
	for (auto& num : a) {
		poly_add(num);
		c.break_point("Poly RBtree", num);
	}
	c.show();
	Clock c2;
	c2.start();
	for (auto& num : a) {
		poly_add(num);
		c2.break_point("Poly Hash", num);
	}
	c2.show();
	// draw
	// show ================================
	GnuplotActor::list_spActor lga;
	lga.push_back(GnuplotActor::Clock_Wall(c));
	lga.push_back(GnuplotActor::Clock_Wall(c2));
	Gnuplot gp;
	gp.set_equal_ratio();
	//gp.set_xrange(2.0,3.0);
	//gp.set_yrange(1.5,2.5);
	gp.set_xlogscale(10);
	gp.set_ylogscale(10);
	gp.plot(lga);
	//delete shape

}

}
#endif
