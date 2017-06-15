#ifndef TEST_SURFACE_CONSTRUCTOR_2D_HPP_
#define TEST_SURFACE_CONSTRUCTOR_2D_HPP_

#include <cstdlib>
#include <iostream>
#include <ctime>
#include <assert.h>
#include <gtest/gtest.h>
#include "io/plotly.h"
#include "io/plotly_actor.h"

#include "ts/ts.h"

using namespace carpio;

namespace TS {

TEST(ts_2d, construction) {
	typedef double vt;
	typedef std::shared_ptr<Surface<vt, 2> > pSur;
	typedef Surface<vt, 2> Sur;
	std::cout << "test construct =====" << std::endl;
	Float r = 1;
	uInt n = 10;
	pSur sur = ConstructCircle2(n, r);
	st count_edge = sur->count_edge();

	//std::cout << " num face = " << count_face << std::endl;
	std::cout << " num edge   = " << count_edge << std::endl;
	std::cout << " num vertex = " << sur->count_vertex() << std::endl;

	//Sur.show_vertex();
	ASSERT_EQ(sur->size_face(), n);
	sur->transfer(1.0, 2.2,3.3);
	//ASSERT_EQ(Sur->size_edge(), 20);

	Plotly p;
	p.add(carpio::PlotlyActor::SurfaceWireFrame(*sur));
	p.plot();
	//std::list<vtkSmartPointer<vtkProp> > actors;
	//actors.push_back(vtk_new_actor(Sur));
	//actors.push_back(vtk_new_actor_normal(Sur));
	//actors.push_back(vtk_new_actor_axes(0, 0, 0));
	// vtk_show_actor(actors);
}

TEST(ts_2d, show) {
	std::cout << "test show =====" << std::endl;

}

}

#endif
