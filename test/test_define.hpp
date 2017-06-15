#ifndef _TEST_DEFINE_H_
#define _TEST_DEFINE_H_

//#include "../io/io_gnuplot_domain.h"

#include "domain/domain.hpp"
#include "../src/calculation/poisson.hpp"
#include "../src/calculation/interpolate.h"
#include "gtest/gtest.h"
#include <random>
#include <math.h>
using namespace std;
namespace carpio {

//inline Float generate_random_number(Float min, Float max) {
//	//srand((unsigned)time(0));
//	Float floor = min, ceiling = max, range = (ceiling - floor);
//	return floor + ((range * rand()) / (RAND_MAX + 1.0));
//}

typedef Domain_<Float, Float, 2> Domain;

inline Float error_1(Domain& d, St ires, St iexact) {
	Float norm1 = 0;
	Float svol = 0;
	for (typename Domain::Grid::iterator_leaf iterf = d.grid().begin_leaf();
			iterf != d.grid().end_leaf(); ++iterf) {
		Float res = iterf->cdva(ires);
		Float exa = iterf->cdva(iexact);
		Float vol = iterf->volume();
		Float err = res - exa;
		norm1 += (Abs(err) * vol);
		svol += vol;
	}
	return norm1 / svol;
}

inline Float error_2(Domain& d, St ires, St iexact) {
	Float norm2 = 0;
	Float svol = 0;
	for (typename Domain::Grid::iterator_leaf iterf = d.grid().begin_leaf();
			iterf != d.grid().end_leaf(); ++iterf) {
		Float res = iterf->cdva(ires);
		Float exa = iterf->cdva(iexact);
		Float vol = iterf->volume();
		Float err = res - exa;
		norm2 += (err * err * vol);
		svol += vol;
	}
	return sqrt(norm2) / svol;
}

inline Float error_i(Domain& d, St ires, St iexact) {
	Float normi = 0;
	for (typename Domain::Grid::iterator_leaf iterf = d.grid().begin_leaf();
			iterf != d.grid().end_leaf(); ++iterf) {
		Float res = iterf->cdva(ires);
		Float exa = iterf->cdva(iexact);
		Float err = res - exa;
		if (iterf == d.grid().begin_leaf()) {
			normi = Abs(err);
		} else {
			if (normi < Abs(err)) {
				normi = Abs(err);
			}
		}
	}
	return normi;
}

inline Float cal_order(Float ec, Float ef) {
	return log(Abs(ec) / Abs(ef)) / log(2);
}



//inline void show_splot_surface(Domain& d, St idx) {
//	std::list<Gnuplot_actor> slga;
//	Gnuplot_actor sga;
//	//GnuplotActor_LeafNodesSurface(sga, d.grid(), idx);
//	//GnuplotActor_NodesSurface(sga, pn, 0);
//	slga.push_back(sga);
//	//GnuplotActor_GhostNodesSurface(sga, d.ghost(), idx);
//	slga.push_back(sga);
//	//sga.show_data();
//	Gnuplot sgp;
//	//sgp.set_equal_ratio();
//	sgp.set_view(45, 10, 1, 1);
//	sgp.set_palette_blue_red();
//	sgp.set("ticslevel 0");
//	//sgp.set_xrange(1.4, 2.0);
//	//sgp.set_yrange(1.4, 2.0);
//	sgp.splot(slga);
//}
//
//inline void show_plot_contour(Domain& d, St idx, std::string title = "") {
//	std::list<Gnuplot_actor> lga;
//	Gnuplot_actor ga;
//	//GnuplotActor_LeafNodesContours(ga, d.grid(), idx);
//	//GnuplotActor_NodesSurface(sga, pn, 0);
//	lga.push_back(ga);
//	//GnuplotActor_GhostNodesContours(ga, d.ghost(), idx);
//	lga.push_back(ga);
//	//sga.show_data();
//	Gnuplot gp;
//	gp.set_title(title);
//	//sgp.set_equal_ratio();
//	//sgp.set_view(45, 10, 1, 1);
//	gp.set_palette_blue_red();
//	//sgp.set("ticslevel 0");
//	//sgp.set_xrange(1.4, 2.0);
//	//sgp.set_yrange(1.4, 2.0);
//	gp.plot(lga);
//}

//inline void show_veo_field(Domain& d, St idxu, St idxv) {
//	std::list<Gnuplot_actor> lga;
//	Gnuplot_actor ga;
//	//GnuplotActor_Velocity_field(ga, d.grid(), idxu, idxv, 1);
//	lga.push_back(ga);
//	//GnuplotActor_LeafNodes(ga, d.grid());
//	lga.push_back(ga);
//
//	//sga.show_data();
//	Gnuplot gp;
//	gp.set_equal_ratio();
//	//sgp.set_view(45, 10, 1, 1);
//	gp.set_palette_blue_red();
//	//sgp.set("ticslevel 0");
//	//sgp.set_xrange(1.4, 2.0);
//	//sgp.set_yrange(1.4, 2.0);
//	gp.plot(lga);
//}

//inline void show_value_on_line(Domain_<Float, Float, 2>& domain, Axes aix,
//		Float v, St i) {
//	typedef Domain_<Float, Float, 2> Domain;
//	typedef Domain_<Float, Float, 2>::pNode pNode;
//	std::list<Gnuplot_actor> lga;
//	Gnuplot_actor ga;
//
//	//GnuplotActor_GhostNodesContours(ga, domain.ghost(), 4);
//	//lga.push_back(ga);
//	//GnuplotActor_GhostNodesDataIndex(ga, domain.ghost());
//	//lga.push_back(ga);
//	//GnuplotActor_LeafNodes(ga, domain.grid(), 4);
//	//lga.push_back(ga);
//	std::list<pNode> lpn = domain.grid().get_leaf(aix, v);
//
//	//GnuplotActor_NodesValues(ga, lpn, i, VerticalAxes2D(aix));
//	//GnuplotActor_NodesContours(ga, lpn, 4);
//	lga.push_back(ga);
//	//GnuplotActor_LeafNodes(ga, domain.grid());
//	//lga.push_back(ga);
//	//GnuplotActor_Shape2D(ga, shape, 0);
//	//lga.push_back(ga);
//	Gnuplot gp;
//	//gp.set_equal_ratio();
//	gp.plot(lga);
//}
//
//inline void output_value_on_line(const std::string& filename,	//
//		Domain_<Float, Float, 2>& domain,  //
//		Axes aix, Float v, St idx) {
//	//typedef Domain_<Float, Float, 2> Domain;
//	typedef Domain_<Float, Float, 2>::pNode pNode;
//	typedef Interpolate_<Float, Float, 2> Inter;
//	typedef Point_<Float, 2> Point;
//	typedef Expression_<Float, Float,2> Exp;
//
//	// build an array
//	St n = 300;
//	ArrayListV<Float> arr = Linspace(  //
//			domain.grid().min(VerticalAxes2D(aix)),  //
//			domain.grid().max(VerticalAxes2D(aix)), n);
//
//	// Interpolate
//	ArrayListV<Float> arrval(n);
//	std::fstream fss;
//	fss.open(filename, std::fstream::out);
//	for (St i = 0; i < n; ++i) {
//		Point p;
//		p[aix] = v;
//		p[VerticalAxes2D(aix)] = arr[i];
//		pNode pn = domain.grid().get_pnode(p.x(), p.y());
//
//		Exp exp = Inter::OnPlane(pn, _X_, _Y_, p.x(), p.y(), 1);
//		arrval[i] = exp.substitute(idx);
//		fss << arr[i] << " " << arrval[i] << "\n";
//	}
//	fss.close();
//}

//inline Float set_v_1(Float x, Float y, Float z) {
//	return 1.0;
//}
//inline Float set_v_0(Float x, Float y, Float z) {
//	return 0.0;
//}

}

#endif
