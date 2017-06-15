#ifndef TEST_SURFACE_CONSTRUCTOR_3D_HPP_
#define TEST_SURFACE_CONSTRUCTOR_3D_HPP_

#include <cstdlib>
#include <iostream>
#include <ctime>
#include <assert.h>
#include <gtest/gtest.h>
#include "io/plotly.h"
#include "io/plotly_actor.h"

namespace carpio {

const std::string FILES = "./test/input_files/";

TEST(ts_3d, construction) {
	typedef double vt;
	std::cout << "test construct =====" << std::endl;
	Float r = 1.0;
	uInt n = 10;
	std::shared_ptr<TS::Surface<Float, 3> > sur = TS::ConstructCircle3(n, r);

	sur->show_vertex();
	ASSERT_EQ(sur->size_face(), 10);
	ASSERT_EQ(sur->size_edge(), 20);

	//Plotly p;
	//p.add(PlotlyActor::Surface(*Sur));
	//p.plot();
}

TEST(ts_3d, cone) {
	typedef double vt;
	std::cout << "test construct =====" << std::endl;
	vt r = 1;
	uInt n = 100;
	vt zbottom = 0;
	vt zpex = 1;
	auto psur = TS::ConstructCone(n, r, zbottom, zpex);

	psur->show_vertex();
	ASSERT_EQ(psur->size_face(), n);
	//ASSERT_EQ(Sur.size_edge(), 20);

	Plotly p;
	auto actor = PlotlyActor::SurfaceWireFrame(*psur);
	actor->set_opacity(0.8);
	p.add(actor);
	//p.plot();
}

TEST(ts_3d, read_head) {
	typedef double vt;
	typedef TS::Surface<vt, 3> Sur;
	typedef std::shared_ptr<TS::Surface<vt, 3> > pSur;
	std::cout << "test construct =====" << std::endl;
	pSur psur(new Sur());
	psur->load_gts_file(FILES + "icosa.gts");
	//Sur.output_vtk("out.vtk");
	std::cout << "num faces " << psur->size_face() << "\n";
	std::cout << "++++++++++--------------------------------\n";

	Plotly p;
	auto actor = PlotlyActor::Surface(*psur);
	p.add(actor);
	actor->set_opacity(0.7);
	auto actor2 = PlotlyActor::SurfaceWireFrame(*psur);
	p.add(actor2);
	//p.plot();
}

TEST(ts_3d, walk) {
	typedef double vt;
	typedef TS::Surface<vt, 3> Sur;
	typedef TS::Face<vt, 3> Fac;
	typedef std::shared_ptr<TS::Surface<vt, 3> > pSur;
	std::cout << "test construct =====" << std::endl;
	pSur psur(new Sur());
	psur->load_gts_file(FILES + "icosa.gts");
	//Sur.output_vtk("out.vtk");
	std::cout << "num faces " << psur->size_face() << "\n";
	std::cout << "++++++++++--------------------------------\n";

	int c = 0;
	typename Sur::Fun_Fac fun = [&c](Sur::Fac& face) {
		std::cout<<" face : "<< c <<"\n";
		c++;
	};
	psur->walk_each_face(fun);

	Plotly p;
	auto actor = PlotlyActor::Surface(*psur);
	p.add(actor);
	actor->set_opacity(0.7);
	auto actor2 = PlotlyActor::SurfaceWireFrame(*psur);
	p.add(actor2);
	//p.add(PlotlyActor::TriangleNormal(*psur));
	p.add(PlotlyActor::TriangleNormal(*(*(psur->begin_face()))));
//	p.plot();
}

TEST(ts_3d, compatiable) {
	typedef double vt;
	typedef TS::Surface<vt, 3> Sur;
	typedef TS::Face<vt, 3> Fac;
	typedef std::shared_ptr<TS::Surface<vt, 3> > pSur;
	std::cout << "test construct =====" << std::endl;
	vt r = 1;
	uInt n = 10;
	vt zbottom = 0;
	vt zpex = 0.1;
	auto psur = TS::ConstructCone(n, r, zbottom, zpex);
	//Sur.output_vtk("out.vtk");
	std::cout << "num faces " << psur->size_face() << "\n";
	std::cout << "++++++++++--------------------------------\n";

	bool res = psur->is_compatible();
	std::cout << " compatible : " << res << "\n";

	psur->revert();

	Plotly p;
	auto actor = PlotlyActor::Surface(*psur);
	p.add(actor);
	actor->set_opacity(0.7);
	auto actor2 = PlotlyActor::SurfaceWireFrame(*psur);
	p.add(actor2);
	p.add(PlotlyActor::TriangleNormal(*psur));
	//p.add(PlotlyActor::TriangleNormal(*(*(psur->begin_face()))));
	p.plot();
}

TEST(ts_3d, compatiable_read) {
	typedef double vt;
	typedef TS::Surface<vt, 3> Sur;
	typedef TS::Face<vt, 3> Fac;
	typedef std::shared_ptr<TS::Surface<vt, 3> > pSur;
	std::cout << "test construct =====" << std::endl;
	pSur psur(new Sur());
	psur->load_gts_file(FILES + "icosa.gts");
	//Sur.output_vtk("out.vtk");
	std::cout << "num faces " << psur->size_face() << "\n";
	std::cout << "++++++++++--------------------------------\n";

	bool res = psur->is_compatible();
	std::cout << " compatible : " << res << "\n";

	//psur->revert();

	Plotly p;
	auto actor = PlotlyActor::Surface(*psur);
	p.add(actor);
	actor->set_opacity(0.7);
	//auto actor2 = PlotlyActor::SurfaceWireFrame(*psur);
	//p.add(actor2);
	p.add(PlotlyActor::TriangleNormal(*psur, 0.2));
	//p.add(PlotlyActor::TriangleNormal(*(*(psur->begin_face()))));
	//p.plot();
}

}

#endif
