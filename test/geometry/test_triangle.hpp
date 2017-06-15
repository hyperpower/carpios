#ifndef __TEST_TRIANGLE_HPP_
#define __TEST_TRIANGLE_HPP_

#include "gtest/gtest.h"
#include "geometry/_triangle.hpp"
#include "geometry/_tri_tri_intersect.hpp"
#include "utility/random.h"
#include "utility/clock.h"

#include <functional>
#include <io/plotly.h>
#include <io/plotly_actor.h>

#include <math.h>
#include <string>
#include <memory>

namespace carpio {

TEST(Triangle, construct) {
	typedef Triangle_<double, 2> Tri2;
	typedef Point_<double, 2> Poi2;
	Poi2 pa(1, 1);
	Poi2 pb(2, 1);
	Poi2 pc(4, 3);

	Tri2 tri(pa, pb, pc);
	tri.show();
	std::cout << "Max x = " << tri.max(_X_) << std::endl;
	std::cout << "Max y = " << tri.max(_Y_) << std::endl;
	std::cout << "Min x = " << tri.min(_X_) << std::endl;
	std::cout << "Min y = " << tri.min(_Y_) << std::endl;

}

TEST(Triangle, intersect) {
	typedef Triangle_<double, 3> Tri;
	typedef Point_<double, 3> Poi;
	Poi pa(0, 0, 0);
	Poi pb(1, 0, 0);
	Poi pc(0.5, 1, 0);
	Tri tri1(pa, pb, pc);

	Poi pa2(0.5, 0.5, -1.0);
	Poi pb2(0.51, 0.0, 1.0);
	Poi pc2(0.52, 1.0, 1.0);
	Tri tri2(pa2, pb2, pc2);

	std::cout << "Is intersect 97 : " << TriTriIsect_97(tri1, tri2)
			<< std::endl;
	std::cout << "Is intersect 02 : " << TriTriIsect_02(tri1, tri2)
			<< std::endl;

	Plotly p;
	p.add(PlotlyActor::Triangle(tri1));
	p.add(PlotlyActor::Triangle(tri2));
	//p.plot();

}




TEST(Triangle, intersectm) {
	typedef Triangle_<double, 3> Tri;
	typedef Point_<double, 3> Poi;

	std::function<Tri(double, double)> fun = [](double min, double max) {
		bool valid = false;
		Tri tri;
		while(!valid) {
			Random::randomSeed();
			for(Tri::size_type i=0;i<3;i++) {
				tri.p(i)[0] = Random::nextDouble(min, max);
				tri.p(i)[1] = Random::nextDouble(min, max);
				tri.p(i)[2] = Random::nextDouble(min, max);
			}
			valid = tri.is_valid();
		}
		return tri;
	};
	double min = 1;
	double max = 10;

	std::list<double> ltime97;
	std::list<double> ltime02;
	std::list<double> lnum;

	std::vector<int> vnum = { 1, 10, 50, 100, 500, 1000, 2000, 3000, 4000, 5000,
			6000, 7000, 10000 };
	for (int n : vnum) {
		tick_t time97 = 0;
		tick_t time02 = 0;
		for (int i = 0; i < n; i++) {
			Tri tri1 = fun(min, max);
			//std::cout << "Is valid     1  : " << tri1.is_valid();

			Tri tri2 = fun(min, max);
			//std::cout << "  |   2  : " << tri2.is_valid() << std::endl;

			tick_t s97 = Clock::Tick();
			int res97 = TriTriIsect_97(tri1, tri2);
			tick_t t97 = Clock::Tick() - s97;
			time97 += t97;

			tick_t s02 = Clock::Tick();
			int res02 = TriTriIsect_02(tri1, tri2);
			tick_t t02 = Clock::Tick() - s02;
			time02 += t02;

			//std::cout << "Is intersect 97 : " << res97;
			//std::cout << "  |   02 : " << res02 << std::endl;
			ASSERT_EQ(res97, res02);
		}
		ltime97.push_back(Clock::TicksToMillisecondsD(time97));
		ltime02.push_back(Clock::TicksToMillisecondsD(time02));
		lnum.push_back(n);
	}

	Plotly p;
	PlotlyActor::spPA actor = PlotlyActor::XY(lnum, ltime97);
	actor->set_name("method 97");
	p.add(actor);
	actor = PlotlyActor::XY(lnum, ltime02);
	actor->set_name("method 02");
	p.add(actor);
	//p.plot();
}

TEST(Triangle, catch_diff) {
	typedef Triangle_<double, 3> Tri;
	typedef Point_<double, 3> Poi;

	std::function<Tri(double, double)> fun = [](double min, double max) {
		bool valid = false;
		Tri tri;
		while(!valid) {
			Random::randomSeed();
			for(Tri::size_type i=0;i<3;i++) {
				tri.p(i)[0] = Random::nextDouble(min, max);
				tri.p(i)[1] = Random::nextDouble(min, max);
				tri.p(i)[2] = Random::nextDouble(min, max);
			}
			valid = tri.is_valid();
		}
		return tri;
	};
	for (int i = 0; i < 1000; i++) {

		double min = 1;
		double max = 10;
		Tri tri1 = fun(min, max);
		//std::cout << "Is valid     1  : " << tri1.is_valid();

		Tri tri2 = fun(min, max);
		//std::cout << "  |   2  : " << tri2.is_valid() << std::endl;

		int res97 = TriTriIsect_97(tri1, tri2);

		int res02 = TriTriIsect_02(tri1, tri2);

		//std::cout << "Is intersect 97 : " << res97;
		//std::cout << "  |   02 : " << res02 << std::endl;
		if (res97 != res02) {
			tri1.show();
			tri2.show();
			Plotly p;
			p.add(PlotlyActor::Triangle(tri1));
			p.add(PlotlyActor::Triangle(tri2));
			p.plot();
			ASSERT_EQ(res97, res02);
			break;
		}
	}
	std::cout << "End catch difference \n";
}

TEST(Triangle, intersect_cal) {
	typedef Triangle_<double, 3> Tri;
	typedef Point_<double, 3> Poi;
	Poi pa(0, 0, 0);
	Poi pb(1, 0, 0);
	Poi pc(0.5, 1, 0);
	Tri tri1(pa, pb, pc);

	Poi pa2(1.5, 0.5, -1.0);
	Poi pb2(1.51, 0.0, 1.0);
	Poi pc2(1.52, 1.0, 1.0);
	Tri tri2(pa2, pb2, pc2);

	std::cout << "Is intersect 02 : " << TriTriIsect_02(tri1, tri2)
			<< std::endl;
	bool cp;
	Poi ps, pe;
	int res = TriTriIsect_Segment(tri1, tri2, cp, ps, pe);
	ps.show();
	pe.show();
	std::cout << "coplane  : " << cp << std::endl;

	Plotly p;
	p.add(PlotlyActor::Triangle(tri1));
	p.add(PlotlyActor::Triangle(tri2));
	if(res == 1){
		p.add(PlotlyActor::Segment(ps, pe));
	}
	//p.plot();

}

TEST(Triangle, intersect_cal_special_case) {
	// this test case show that the intersection point is one point
	// the function will reture two equal points
	typedef Triangle_<double, 3> Tri;
	typedef Point_<double, 3> Poi;
	Poi pa(0, 0, 0);
	Poi pb(1, 0, 0);
	Poi pc(0.5, 1, 0);
	Tri tri1(pa, pb, pc);

	Poi pa2(0.5, 0.5,  0.0);
	Poi pb2(0.51, 0.0, 1.0);
	Poi pc2(0.52, 1.0, 1.0);
	Tri tri2(pa2, pb2, pc2);

	std::cout << "Is intersect 02 : " << TriTriIsect_02(tri1, tri2)
			<< std::endl;
	bool cp;
	Poi ps, pe;
	int res = TriTriIsect_Segment(tri1, tri2, cp, ps, pe);
	ps.show();
	pe.show();
	std::cout << "coplane  : " << cp << std::endl;

	Plotly p;
	p.add(PlotlyActor::Triangle(tri1));
	p.add(PlotlyActor::Triangle(tri2));
	if(res == 1){
		p.add(PlotlyActor::Segment(ps, pe));
	}
	//p.plot();

}

TEST(Triangle, intersect_cal_special_case2) {
	typedef Triangle_<double, 3> Tri;
	typedef Point_<double, 3> Poi;
	Poi pa(0, 0, 0);
	Poi pb(1, 0, 0);
	Poi pc(0.5, 1, 0);
	Tri tri1(pa, pb, pc);

	Poi pa2(0.5, 0.5,  -1.0);
	Poi pb2(0.51, 0.0, 0.0);
	Poi pc2(0.52, 1.0, 0.0);
	Tri tri2(pa2, pb2, pc2);

	std::cout << "Is intersect 02 : " << TriTriIsect_02(tri1, tri2)
			<< std::endl;
	bool cp;
	Poi ps, pe;
	int res = TriTriIsect_Segment(tri1, tri2, cp, ps, pe);
	ps.show();
	pe.show();
	std::cout << "coplane  : " << cp << std::endl;

	Plotly p;
	p.add(PlotlyActor::Triangle(tri1));
	p.add(PlotlyActor::Triangle(tri2));
	if(res == 1){
		p.add(PlotlyActor::Segment(ps, pe));
	}
	p.plot();
}

}
#endif
