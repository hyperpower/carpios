/************************
 //  \file   test_AABBox.h
 //  \brief
 // 
 //  \author czhou
 //  \date   17 juin 2015 
 ***********************/
#ifndef TEST_AABBOX_H_
#define TEST_AABBOX_H_
#include <iostream>
#include "ts/ts_AABBox.h"
#include "ts/ts_BBTree.h"
#include "ts/ts_intersect.h"
#include "ts/ts_triangle.h"
//#include "ts/ts_io.h"
#include "ts/ts_surface.h"
#include "ts/ts_surface_constructor.h"
#include "io/plotly.h"
#include "io/plotly_actor.h"

#include <cstdlib>
#include <iostream>
#include <ctime>
#include <assert.h>
#include <gtest/gtest.h>

using namespace std;
namespace TS {

const std::string FILES = "./test/input_files/";

TEST(ts_AABBox, read) {
	typedef double vt;
	typedef TS::Surface<vt, 3> Sur;
	//typedef std::shared_ptr<TS::Surface<vt, 3> > pSur;
	//pSur psur(new Sur());
	//Sur sur(FILES + "head.gts");
	auto psur = ConstructFromGtsFile(FILES + "head.gts", 0.0);
	auto bsur = ConstructFromGtsFile(FILES + "icosa.gts", 0.0);
	bsur->transfer(3.0, 1.0, 1.0);
	//Sur.output_vtk("out.vtk");
	//cout << "num faces " << Sur.faces.size() << "\n";
	//cout << "++++++++++--------------------------------\n";
	Set<AABBox<vt, 3> > set_box;
	int i = 0;
	for (auto iter = psur->begin_face(); iter != psur->end_face(); ++iter) {
		auto ptr = (*iter);   //pFac
		AABBox<vt, 3> tbox(ptr.get());  //pTri
		set_box.insert(tbox);
		i++;
	}
	Set<AABBox<vt, 3> > set_box2;
	i = 0;
	for (auto iter = bsur->begin_face(); iter != bsur->end_face(); ++iter) {
		auto ptr = (*iter);   //pFac
		AABBox<vt, 3> tbox(ptr.get());  //pTri
		set_box.insert(tbox);
		i++;
	}
	/*************** show a box
	 auto pf = *(psur->begin_face());   //pFac
	 AABBox<vt, 3> tb(pf.get());        //pTri
	 carpio::Plotly p;
	 p.add(carpio::PlotlyActor::Surface(*psur));
	 p.add(carpio::PlotlyActor::AABBox(tb));t
	 p.plot();
	 */

	std::cout << "num boxes " << set_box.size() << "\n";

	BBTree<AABBox<vt, 3> > bbtree(set_box);

	carpio::Plotly p;
	p.add(carpio::PlotlyActor::Surface(*psur));
	p.add(carpio::PlotlyActor::Surface(*bsur));
	p.add(carpio::PlotlyActor::BBTree(bbtree, 3));
//	p.plot();

	//bbtree.output_vtk_height("h.vtk", 1);
	//cout << "tree size " << bbtree.size() << "\n";
	//vtk_show(bbtree,bbtree.height()-1);
	//set_box.begin()->output_vtk("begin.vtk");
	//std::list<vtkSmartPointer<vtkProp> > actors;
	//actors.push_back(vtk_new_actor(bbtree.begin()->box));
	//actors.push_back(vtk_new_actor_axes(0, 0, 0));
	//vtk_show_actor(actors);
	std::cout << "is Orientable  " << psur->is_orientable() << std::endl;
	//Surface<Float, 3> Sur_t("triangle.gts");
	//Sur_t.output_vtk("out_tri.vtk");
	//auto iter = Sur_t.faces.begin();
	//AABBox<Float, 3> bs((*iter));
	//bs.output_vtk("bs.vtk");

	//cout << "inter " << do_intersect_box_box(&bbtree, &bs) << endl;
	//List<AABBox<Float, 3>*> lres;
	//do_intersect_box_obj(&bbtree, &bs, lres);
	//std::list<vtkSmartPointer<vtkProp> > actors;
	//actors.push_back(vtk_new_actor(lres));
	//actors.push_back(vtk_new_actor(Sur));
	//actors.push_back(vtk_new_actor_normal(Sur));

	//actors.push_back(vtk_new_actor_axes(0, 0, 0));
	//vtk_show_actor(actors);

	//output_vtk("inter.vtk", lres);
	//cout << "end test ===============\n";

}

TEST(AABox, boxtri) {
	// tri and its box -------------------------
	auto psur = ConstructFromGtsFile(FILES + "icosa.gts", 0.0);
	auto iter = psur->begin_face();
	//std::advance(iter, 4);
	AABBox<Float, 3> bt((*iter).get());

	// box set ---------------------------------
	AABBox<Float, 3> bs(EMPTY, nullptr, 0, 0, 0, 1, 2, 3);

	carpio::Plotly p;
	p.add(carpio::PlotlyActor::Surface(*psur));
	p.add(carpio::PlotlyActor::Face(*(iter->get())));
	p.add(carpio::PlotlyActor::AABBox(bs));
	//p.plot();
	std::cout << "interset tri :" << bs.are_overlapping((*iter).get())
			<< std::endl;
	std::cout << "interset box :" << bs.are_overlapping(bt) << std::endl;
	cout << "end test ===============\n";
}

TEST(AABox, box_tree) {
	// tri and its box -------------------------
	typedef Float vt;
	auto psur = ConstructFromGtsFile(FILES + "icosa.gts", 0.0);
	auto iter = psur->begin_face();
	//std::advance(iter, 4);
	Set<AABBox<vt, 3> > set_box;
	int i = 0;
	for (auto iter = psur->begin_face(); iter != psur->end_face(); ++iter) {
		auto ptr = (*iter);   //pFac
		AABBox<vt, 3> tbox(ptr.get());  //pTri
		set_box.insert(tbox);
		i++;
	}
	BBTree<AABBox<vt, 3> > bbtree(set_box);

	// box set ---------------------------------
	AABBox<Float, 3> bs(EMPTY, nullptr, 0, 0, 0, 1, 2, 3);

	carpio::Plotly p;
	p.add(carpio::PlotlyActor::Surface(*psur));
	p.add(carpio::PlotlyActor::BBTree(bbtree, 0));
	p.add(carpio::PlotlyActor::AABBox(bs));
	//p.plot();
	cout << "end test ===============\n";
}

void test_triangle() {
	Float v0[3];
	Float v1[3];
	Float v2[3];

	Float u0[3];
	Float u1[3];
	Float u2[3];
	std::srand(std::time(0)); // use current time as seed for random generator
	for (int c = 0; c < 1; c++) {
		cout << c << " !!!--------------------------- \n";
		for (int i = 0; i < 3; i++) {
			v0[i] = std::rand() / float(RAND_MAX);
			v1[i] = std::rand() / float(RAND_MAX);
			v2[i] = std::rand() / float(RAND_MAX);
			u0[i] = std::rand() / float(RAND_MAX);
			u1[i] = std::rand() / float(RAND_MAX);
			u2[i] = std::rand() / float(RAND_MAX);
		}

		int res1, res2;
		//res1 = TriTriIsect_Guigue(v0, v1, v2, u0, u1, u2);
		//res2 = NoDivTriTriIsect(v0, v1, v2, u0, u1, u2);
		//return code change
		if (res1 == -1 && res2 == 0) {
			cout << "pass \n";
		} else if (res1 == 1 && res2 == 1) {
			cout << "pass \n";
		} else {
			//output_vtk("t1.vtk", v0, v1, v2);
			//output_vtk("t2.vtk", u0, u1, u2);
			cout << res1 << "  " << res2 << endl;
			exit(0);
		}
	}
}

void test_triangle_cal() {
	Float v0[3];
	Float v1[3];
	Float v2[3];

	Float u0[3];
	Float u1[3];
	Float u2[3];
	std::srand(std::time(0)); // use current time as seed for random generator
	for (int c = 0; c < 1; c++) {
		cout << c << " !!!--------------------------- \n";
		for (int i = 0; i < 3; i++) {
			v0[i] = std::rand() / float(RAND_MAX);
			v1[i] = std::rand() / float(RAND_MAX);
			v2[i] = std::rand() / float(RAND_MAX);
			u0[i] = std::rand() / float(RAND_MAX);
			u1[i] = std::rand() / float(RAND_MAX);
			u2[i] = std::rand() / float(RAND_MAX);
		}

		int res1, res2;
		Float* x = NULL;
		Float* y = NULL;
		Float* z = NULL;
		short len = 0;
		//res1 = TriTriIsect_Guigue_calculation(v0, v1, v2, u0, u1, u2, x, y, z,
		//		len);
		cout << "Len " << len << endl;
		Float resv0[] = { x[0], y[0], z[0] };
		Float resv1[] = { x[1], y[1], z[1] };
		//output_vtk("seg.vtk", resv0, resv1);

		//res2 = NoDivTriTriIsect(v0, v1, v2, u0, u1, u2);
		//return code change
		if (res1 == -1 && res2 == 0) {
			cout << "pass \n";
		} else if (res1 == 1 && res2 == 1) {
			cout << "pass \n";
		} else {
			cout << res1 << "  " << res2 << endl;
			exit(0);
		}
		//output_vtk("t1.vtk", v0, v1, v2);
		//output_vtk("t2.vtk", u0, u1, u2);
	}
}

void test_triangle_sigle() {
	Float v0[3];
	Float v1[3];
	Float v2[3];

	Float u0[3];
	Float u1[3];
	Float u2[3];
	std::srand(std::time(0)); // use current time as seed for random generator

	cout << " !!!--------------------------- \n";
	v0[0] = 0;
	v0[1] = 0;
	v0[2] = 0;
	v1[0] = 1;
	v1[1] = 0;
	v1[2] = 0;
	v2[0] = 0;
	v2[1] = 1;
	v2[2] = 0;

	u0[0] = 0.25;
	u0[1] = 0.25;
	u0[2] = 1;
	u1[0] = 0.5;
	u1[1] = 0;
	u1[2] = -1;
	u2[0] = 0;
	u2[1] = 0.5;
	u2[2] = -1;

	int res1, res2;
	Float* x = NULL;
	Float* y = NULL;
	Float* z = NULL;
	short len = 0;
	//res1 = TriTriIsect_Guigue_calculation(v0, v1, v2, u0, u1, u2, x, y, z, len);
	cout << "Len " << len << endl;
	Float resv0[] = { x[0], y[0], z[0] };
	Float resv1[] = { x[1], y[1], z[1] };
	//output_vtk("seg.vtk", resv0, resv1);
	//res2 = NoDivTriTriIsect(v0, v1, v2, u0, u1, u2);
	//return code change
	if (res1 == -1 && res2 == 0) {
		cout << "pass \n";
	} else if (res1 == 1 && res2 == 1) {
		cout << "pass \n";
	} else {
		//	output_vtk("t1.vtk", v0, v1, v2);
		//	output_vtk("t2.vtk", u0, u1, u2);
		cout << res1 << "  " << res2 << endl;
		exit(0);
	}
	//output_vtk("t1.vtk", v0, v1, v2);
	//output_vtk("t2.vtk", u0, u1, u2);

}

}

#endif /* TEST_AABBOX_H_ */
