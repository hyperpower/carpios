#ifndef TEST_OVERLAP_H_
#define TEST_OVERLAP_H_
#include <iostream>
#include "ts/ts.h"

#include "io/plotly.h"
#include "io/plotly_actor.h"
#include <gtest/gtest.h>

namespace TS {

const std::string FILES = "./test/input_files/";

TEST(AABox, box_tree) {
	// tri and its box -------------------------
	typedef Float vt;
	typedef Creation<vt, 3> Cr;
	auto psur = Cr::FromGtsFile(FILES + "icosa.gts");
	auto iter = psur->begin_face();
	//std::advance(iter, 4);
	Set<AABBox<vt, 3> > set_box;
	int i = 0;
	for (auto iter = psur->begin_face(); iter != psur->end_face(); ++iter) {
		auto ptr = (*iter);   //pFac
		AABBox<vt, 3> tbox(ptr.get());//pTri
		set_box.insert(tbox);
		i++;
	}
	BBTree<AABBox<vt, 3> > bbtree(set_box);

	// box set ---------------------------------
	AABBox<Float, 3> bs(EMPTY, nullptr, 3, 3, 3, 4, 5, 6);

	bool res = Intersection<vt, 3>::DetectBox(bbtree, bs);
	std::cout<< "overlap :" << res << std::endl;

	carpio::Plotly p;
	p.add(carpio::PlotlyActor::Surface(*psur));
	p.add(carpio::PlotlyActor::BBTreeLevel(bbtree, 1));
	p.add(carpio::PlotlyActor::AABBox(bs));
	//p.plot();
	std::cout << "end test ===============\n";
}

TEST(AABox, box_tree_tree) {
	// tri and its box -------------------------
	typedef Float vt;
	typedef Creation<vt, 3> Cr;
	auto psur = Cr::FromGtsFile(FILES + "icosa.gts");
	psur->transfer(2,0,0);
	Set<AABBox<vt, 3> > set_box1;
	for (auto iter = psur->begin_face(); iter != psur->end_face(); ++iter) {
		auto ptr = (*iter);   //pFac
		AABBox<vt, 3> tbox(ptr.get());//pTri
		set_box1.insert(tbox);
	}
	BBTree<AABBox<vt, 3> > bbtree1(set_box1);

	// box set ---------------------------------
	typedef Float vt;
	typedef Creation<vt, 3> Cr;
	auto bsur = Cr::FromGtsFile(FILES + "head.gts");
	Set<AABBox<vt, 3> > set_box2;
	for (auto iter = bsur->begin_face(); iter != bsur->end_face(); ++iter) {
		auto ptr = (*iter);   //pFac
		AABBox<vt, 3> tbox(ptr.get());//pTri
		set_box2.insert(tbox);
	}
	BBTree<AABBox<vt, 3> > bbtree2(set_box2);
	// test box overlapping
	bool res = Intersection<vt, 3>::DetectBox(bbtree1, bbtree2);
	std::cout<< "overlap box :" << res << std::endl;
	// test the object in the box overlapping
	// here, the object is triangle
	bool res2 = Intersection<vt, 3>::DetectObj(bbtree1, bbtree2);
	std::cout<< "overlap obj :" << res2 << std::endl;

	carpio::Plotly p;
	p.add(carpio::PlotlyActor::Surface(*psur));
	p.add(carpio::PlotlyActor::BBTree(bbtree1, 0));
	p.add(carpio::PlotlyActor::BBTree(bbtree2, 0));
	//p.plot();
	std::cout << "end test ===============\n";
}

TEST(AABox, tri_tri) {
	// tri and its box -------------------------
	typedef Float vt;
	typedef Creation<vt, 3> Cr;
	auto psur = Cr::FromGtsFile(FILES + "tri1.gts");
	auto pf1 = *(psur->begin_face());
	//psur->transfer(2,0,0);
	AABBox<vt, 3> box1(pf1.get());//pTri

	auto bsur = Cr::FromGtsFile(FILES + "tri2.gts");
	auto pf2 = *(bsur->begin_face());
	AABBox<vt, 3> box2(pf2.get());//pTri

	bool res = Intersection<vt, 3>::DetectObj(box1, box2);
	Point<vt, 3> pa, pb;
	bool res2= Intersection<vt, 3>::CalInersectBox(box1, box2, pa, pb);
	std::cout<< "Point a" << std::endl;
	pa.show();
	std::cout<< "Point b" << std::endl;
	pb.show();
	std::cout<< "touch one p     :" << res << std::endl;
	std::cout<< "touch one p res2:" << res2 << std::endl;

	// ----------------------------------------
	bsur->transfer(0.0, 0.0, 1.0);
	pf2 = *(bsur->begin_face());
	box2.reset(pf2.get());
	res = Intersection<vt, 3>::DetectObj(box1, box2);
	res2= Intersection<vt, 3>::CalInersectBox(box1, box2, pa, pb);
	std::cout<< "Point a" << std::endl;
	pa.show();
	std::cout<< "Point b" << std::endl;
	pb.show();
	std::cout<< "not touch      :" << res << std::endl;
	std::cout<< "not touch res2 :" << res2 << std::endl;
	// ----------------------------------------
	bsur->transfer(0.0, 0.0, -1.3);
	pf2 = *(bsur->begin_face());
	box2.reset(pf2.get());
	res = Intersection<vt, 3>::DetectObj(box1, box2);
	res2= Intersection<vt, 3>::CalInersectBox(box1, box2, pa, pb);
	std::cout<< "Point a" << std::endl;
	pa.show();
	std::cout<< "Point b" << std::endl;
	pb.show();
	std::cout<< "3              :" << res << std::endl;
	std::cout<< "3              :" << res2 << std::endl;
	// ----------------------------------------
	bsur->transfer(0.0, 0.0, -0.7);
	pf2 = *(bsur->begin_face());
	box2.reset(pf2.get());
	res = Intersection<vt, 3>::DetectObj(box1, box2);
	res2= Intersection<vt, 3>::CalInersectBox(box1, box2, pa, pb);
	std::cout<< "Point a" << std::endl;
	pa.show();
	std::cout<< "Point b" << std::endl;
	pb.show();
	std::cout<< "edge touch     :" << res  << std::endl;
	std::cout<< "edge touch     :" << res2 << std::endl;

	carpio::Plotly p;
	p.add(carpio::PlotlyActor::Surface(*psur));
	p.add(carpio::PlotlyActor::Surface(*bsur));
	//p.add(carpio::PlotlyActor::BBTree(bbtree1, 0));
	//p.add(carpio::PlotlyActor::BBTree(bbtree2, 0));
	//p.plot();
	std::cout << "end test ===============\n";
}

TEST(AABox, box_tree_tree_vector) {
	// tri and its box -------------------------
	typedef Float vt;
	typedef Creation<vt, 3> Cr;
	auto psur = Cr::FromGtsFile(FILES + "icosa.gts");
	psur->transfer(2,0,0);
	Set<AABBox<vt, 3> > set_box1;
	for (auto iter = psur->begin_face(); iter != psur->end_face(); ++iter) {
		auto ptr = (*iter);           //pFac
		AABBox<vt, 3> tbox(ptr.get());//pTri
		set_box1.insert(tbox);
	}
	BBTree<AABBox<vt, 3> > bbtree1(set_box1);

	// box set ---------------------------------
	typedef Float vt;
	auto bsur = Cr::FromGtsFile(FILES + "head.gts");
	Set<AABBox<vt, 3> > set_box2;
	for (auto iter = bsur->begin_face(); iter != bsur->end_face(); ++iter) {
		auto ptr = (*iter);   //pFac
		AABBox<vt, 3> tbox(ptr.get());//pTri
		set_box2.insert(tbox);
	}
	BBTree<AABBox<vt, 3> > bbtree2(set_box2);
	std::vector<AABBox<vt, 3> > vec1, vec2;
	// test the object in the box overlapping
	// here, the object is triangle
	// get all of the intersection boxes in two list;
	bool res2 = Intersection<vt, 3>::DetectObj(bbtree1, bbtree2, vec1, vec2);
	std::cout<< "overlap obj :" << res2 << std::endl;

	carpio::Plotly p;
	//p.add(carpio::PlotlyActor::Surface(*psur));
	//p.add(carpio::PlotlyActor::Surface(*bsur));
	//p.add(carpio::PlotlyActor::BBTree(bbtree1, 0));
	//p.add(carpio::PlotlyActor::BBTree(bbtree2, 0));
	//p.add(carpio::PlotlyActor::AABBox(vec1));
	//p.add(carpio::PlotlyActor::AABBox(vec2));
	p.add(carpio::PlotlyActor::Triangle(vec1));
	p.add(carpio::PlotlyActor::Triangle(vec2));
	p.add(carpio::PlotlyActor::TriangleNormal(vec1));
	p.add(carpio::PlotlyActor::TriangleNormal(vec2));
	//p.plot();
	std::cout << "end test ===============\n";
}

TEST(AABox, box_vector_vector) {
	// tri and its box -------------------------
	typedef Float vt;
	typedef Creation<vt, 3> Cr;
	auto psur = Cr::FromGtsFile(FILES + "icosa.gts");
	psur->transfer(2,0,0);
	Set<AABBox<vt, 3> > set_box1;
	for (auto iter = psur->begin_face(); iter != psur->end_face(); ++iter) {
		auto ptr = (*iter);           //pFac
		AABBox<vt, 3> tbox(ptr.get());//pTri
		set_box1.insert(tbox);
	}
	BBTree<AABBox<vt, 3> > bbtree1(set_box1);

	// box set ---------------------------------
	typedef Float vt;
	auto bsur = Cr::FromGtsFile(FILES + "head.gts");
	Set<AABBox<vt, 3> > set_box2;
	for (auto iter = bsur->begin_face(); iter != bsur->end_face(); ++iter) {
		auto ptr = (*iter);   //pFac
		AABBox<vt, 3> tbox(ptr.get());//pTri
		set_box2.insert(tbox);
	}
	BBTree<AABBox<vt, 3> > bbtree2(set_box2);
	std::vector<AABBox<vt, 3> > vec1, vec2;
	// test the object in the box overlapping
	// here, the object is triangle
	// get all of the intersection boxes in two list;
	bool res2 = Intersection<vt, 3>::DetectObj(bbtree1, bbtree2, vec1, vec2);


	std::cout<< "overlap obj :" << res2 << std::endl;

	carpio::Plotly p;
	//p.add(carpio::PlotlyActor::Surface(*psur));
	//p.add(carpio::PlotlyActor::Surface(*bsur));
	//p.add(carpio::PlotlyActor::BBTree(bbtree1, 0));
	//p.add(carpio::PlotlyActor::BBTree(bbtree2, 0));
	//p.add(carpio::PlotlyActor::AABBox(vec1));
	//p.add(carpio::PlotlyActor::AABBox(vec2));
	p.add(carpio::PlotlyActor::Triangle(vec1));
	p.add(carpio::PlotlyActor::Triangle(vec2));
	p.add(carpio::PlotlyActor::TriangleNormal(vec1));
	p.add(carpio::PlotlyActor::TriangleNormal(vec2));
	p.plot();
	std::cout << "end test ===============\n";
}



}

#endif
