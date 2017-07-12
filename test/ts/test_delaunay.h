#ifndef TEST_DELAUNAY_H_
#define TEST_DELAUNAY_H_

#include <iostream>
#include "ts/ts_AABBox.h"
#include "ts/ts_BBTree.h"
#include "ts/ts_intersect.h"
#include "ts/ts_triangle.h"
#include "ts/ts_delaunay.h"
#include "ts/ts_surface.h"
#include "ts/ts_creation.h"
#include "io/plotly.h"
#include "io/plotly_actor.h"

#include <cstdlib>
#include <iostream>
#include <ctime>
#include <assert.h>
#include <gtest/gtest.h>

using namespace std;

namespace TS {

TEST(ts_delaunay, read) {
	typedef Delaunay<double, 2> D2;
	typedef D2::list_spcPoi LV;
	typedef D2::spcPoi spcPoi;
	typedef D2::Poi Poi;

	int nx = 5;
	int ny = 5;
	double cosa = cos(0.3);
	double sina = sin(0.4);

	LV listv;
	for (int i = 0; i < nx; i++) {
		double x = (double) i / (double) (nx - 1);
		for (int j = 0; j < ny; j++) {
			double y = (double) j / (double) (nx - 1);
			spcPoi spp(new Poi(cosa * x - sina * y, sina * x + cosa * y, 0.0));
			listv.push_back(spp);
		}
	}
	D2 d;
	auto supert = d.triangle_enclosing(listv, 1.0);
	D2::Sur sur;
	sur.add_face(supert);


	carpio::Plotly p;
	p.add(carpio::PlotlyActor::ScatterPoints(listv));
	p.add(carpio::PlotlyActor::SurfaceWireFrame(sur));
	//p.add(carpio::PlotlyActor::BBTree(bbtree, 3));
	//p.plot();
}

}

#endif
