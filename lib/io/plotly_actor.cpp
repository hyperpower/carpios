#include "plotly_actor.h"

namespace carpio {

namespace PlotlyActor {

spPA_scatter XY(const ArrListd& x, const ArrListd& y) {
	spPA_scatter res = spPA_scatter(new Plotly_actor_scatter(x, y));
	return res;
}
spPA_scatter XY(const std::list<double>& x, const std::list<double>& y) {
	spPA_scatter res = spPA_scatter(new Plotly_actor_scatter(x, y));
	return res;
}

spPA_scatter3d Segment(const Poid3& s, const Poid3& e) {
	ArrListd ax(2);
	ArrListd ay(2);
	ArrListd az(2);
	ax[0] = s[0];
	ax[1] = e[0];
	ay[0] = s[1];
	ay[1] = e[1];
	az[0] = s[2];
	az[1] = e[2];
	spPA_scatter3d res = spPA_scatter3d(new Plotly_actor_scatter3d(ax, ay, az));
	return res;
}

spPA_mesh3d Triangle(const Trid3& r) {
	ArrListd ax(3);
	ArrListd ay(3);
	ArrListd az(3);
	for (Trid3::size_type i = 0; i < 3; i++) {
		ax[i] = r.p(i, _X_);
		ay[i] = r.p(i, _Y_);
		az[i] = r.p(i, _Z_);
	}
	spPA_mesh3d res = spPA_mesh3d(new Plotly_actor_mesh3d(ax, ay, az));
	return res;
}

spPA_mesh3d Triangle(const std::vector<AABBoxd3>& vbox) {
	Listd lx, li;
	Listd ly, lj;
	Listd lz, lk;
	int count = 0;
	for (auto iter = vbox.begin(); iter != vbox.end(); ++iter) {
		AABBoxd3 box = (*iter);
		ASSERT(box.type == TS::TRIANGLE || box.type == TS::FACE);
		typename AABBoxd3::pTri pf = CAST(typename AABBoxd3::pTri,
				box.utp_obj());
		for (Surd3::size_type i = 0; i < 3; i++) {
			lx.push_back(pf->get_vertex(i)->x());
			ly.push_back(pf->get_vertex(i)->y());
			lz.push_back(pf->get_vertex(i)->z());
		}
		li.push_back(count);
		lj.push_back(count + 1);
		lk.push_back(count + 2);
		count += 3;
	}

	spPA_mesh3d res = spPA_mesh3d(new Plotly_actor_mesh3d(lx, ly, lz));
	res->set_ijk(li, lj, lk);
	return res;
}

spPA_scatter3d TriangleNormal(const std::vector<AABBoxd3>& vbox) {
	Listd lx;
	Listd ly;
	Listd lz;
	int count = 0;
	for (auto iter = vbox.begin(); iter != vbox.end(); ++iter) {
		AABBoxd3 box = (*iter);
		ASSERT(box.type == TS::TRIANGLE || box.type == TS::FACE);
		typename AABBoxd3::pTri pt = CAST(typename AABBoxd3::pTri,
				box.utp_obj());
		typename AABBoxd3::Poi pc = pt->centroid();
		typename AABBoxd3::Poi pn = pt->normal();
		lx.push_back(pc[0]);
		ly.push_back(pc[1]);
		lz.push_back(pc[2]);
		lx.push_back(pn[0] + pc[0]);
		ly.push_back(pn[1] + pc[1]);
		lz.push_back(pn[2] + pc[2]);

		count += 1;
	}

	spPA_scatter3d res = spPA_scatter3d(
			new Plotly_actor_scatter3d(lx, ly, lz, 2));
	res->set_mode("lines");
	return res;
}

spPA_scatter3d TriangleNormal(const Surd3& surface, double scale) {
	Listd lx;
	Listd ly;
	Listd lz;
	int count = 0;
	for (auto iter = surface.begin_face(); iter != surface.end_face(); ++iter) {
		typename Surd3::spFac spf = (*iter);
		typename AABBoxd3::Poi pc = spf->centroid();
		typename AABBoxd3::Poi pn = spf->normal();
		lx.push_back(pc[0]);
		ly.push_back(pc[1]);
		lz.push_back(pc[2]);
		lx.push_back(pc[0] + pn[0] * scale);
		ly.push_back(pc[1] + pn[1] * scale);
		lz.push_back(pc[2] + pn[2] * scale);

		count += 1;
	}

	spPA_scatter3d res = spPA_scatter3d(
			new Plotly_actor_scatter3d(lx, ly, lz, 2));
	res->set_mode("lines");
	return res;
}

spPA_scatter3d TriangleNormal(const Faced3& face) {
	Listd lx;
	Listd ly;
	Listd lz;
	int count = 0;
	typename AABBoxd3::Poi pc = face.centroid();
	typename AABBoxd3::Poi pn = face.normal();
	lx.push_back(pc[0]);
	ly.push_back(pc[1]);
	lz.push_back(pc[2]);
	lx.push_back(pn[0] + pc[0]);
	ly.push_back(pn[1] + pc[1]);
	lz.push_back(pn[2] + pc[2]);

	count += 1;

	spPA_scatter3d res = spPA_scatter3d(
			new Plotly_actor_scatter3d(lx, ly, lz, 2));
	res->set_mode("lines");
	return res;
}

spPA_mesh3d Face(const Faced3& r) {
	ArrListd ax(3);
	ArrListd ay(3);
	ArrListd az(3);
	for (Trid3::size_type i = 0; i < 3; i++) {
		ax[i] = r.get_vertex(i, TS::_X);
		ay[i] = r.get_vertex(i, TS::_Y);
		az[i] = r.get_vertex(i, TS::_Z);
	}
	spPA_mesh3d res = spPA_mesh3d(new Plotly_actor_mesh3d(ax, ay, az));
	return res;
}

spPA_mesh3d Surface(const Surd3& sur) {
	Listd lx, li;
	Listd ly, lj;
	Listd lz, lk;
	int count = 0;
	for (auto iter = sur.faces.begin(); iter != sur.faces.end(); ++iter) {
		typename Surd3::spFac pf = (*iter);
		for (Surd3::size_type i = 0; i < 3; i++) {
			lx.push_back(pf->get_vertex(i)->x());
			ly.push_back(pf->get_vertex(i)->y());
			lz.push_back(pf->get_vertex(i)->z());
		}
		li.push_back(count);
		lj.push_back(count + 1);
		lk.push_back(count + 2);
		count += 3;
	}

	spPA_mesh3d res = spPA_mesh3d(new Plotly_actor_mesh3d(lx, ly, lz));
	res->set_ijk(li, lj, lk);
	return res;
}

spPA_scatter3d ScatterPoints(const list_spcPoid2& list_p) {
	Listd lx;
	Listd ly;
	Listd lz;
	for (auto iter = list_p.begin(); iter != list_p.end(); ++iter) {
		auto pf = (*iter);
		for (Surd3::size_type i = 0; i < 3; i++) {
			lx.push_back(pf->x());
			ly.push_back(pf->y());
			lz.push_back(0.0);
		}
	}
	spPA_scatter3d res = spPA_scatter3d(
			new Plotly_actor_scatter3d(lx, ly, lz, 1));
	res->set_mode("points");
	return res;
}

spPA_scatter3d SurfaceWireFrame(const Trid2& tri) {
	Listd lx;
	Listd ly;
	Listd lz;
	//for (auto iter = sur.faces.begin(); iter != sur.faces.end(); ++iter) {
	//	typename Surd3::spFac pf = (*iter);
	for (Trid3::size_type i = 0; i < 3; i++) {
		lx.push_back(tri.get_vertex(i)->x());
		ly.push_back(tri.get_vertex(i)->y());
		lz.push_back(0.0);
	}
	lx.push_back(tri.get_vertex(0)->x());
	ly.push_back(tri.get_vertex(0)->y());
	lz.push_back(0.0);

	//}
	spPA_scatter3d res = spPA_scatter3d(
			new Plotly_actor_scatter3d(lx, ly, lz, 4));
	res->set_mode("lines");
	return res;
}

spPA_scatter3d SurfaceWireFrame(const Surd3& sur) {
	Listd lx;
	Listd ly;
	Listd lz;
	for (auto iter = sur.faces.begin(); iter != sur.faces.end(); ++iter) {
		typename Surd3::spFac pf = (*iter);
		for (Surd3::size_type i = 0; i < 3; i++) {
			lx.push_back(pf->get_vertex(i)->x());
			ly.push_back(pf->get_vertex(i)->y());
			lz.push_back(pf->get_vertex(i)->z());
		}
	}
	spPA_scatter3d res = spPA_scatter3d(
			new Plotly_actor_scatter3d(lx, ly, lz, 3));
	res->set_mode("lines");
	return res;

}

spPA_scatter SurfaceWireFrame(const Surd2& sur) {
	Listd lx;
	Listd ly;
	for (auto iter = sur.faces.begin(); iter != sur.faces.end(); ++iter) {
		typename Surd2::spFac pf = (*iter);
		for (Surd3::size_type i = 0; i < 3; i++) {
			lx.push_back(pf->get_vertex(i)->x());
			ly.push_back(pf->get_vertex(i)->y());
		}
		lx.push_back(pf->get_vertex(0)->x());
		ly.push_back(pf->get_vertex(0)->y());
	}
	spPA_scatter res = spPA_scatter(new Plotly_actor_scatter(lx, ly, 4));
	res->set_mode("lines");
	return res;
}

const TS::Location _ORDER[24][3] = { { TS::_M, TS::_M, TS::_M }, //
		{ TS::_P, TS::_M, TS::_M }, //
		{ TS::_P, TS::_M, TS::_M }, //
		{ TS::_P, TS::_P, TS::_M }, //
		{ TS::_P, TS::_P, TS::_M }, //
		{ TS::_M, TS::_P, TS::_M }, //
		{ TS::_M, TS::_P, TS::_M }, //
		{ TS::_M, TS::_M, TS::_M }, //

		{ TS::_M, TS::_M, TS::_P }, //
		{ TS::_P, TS::_M, TS::_P }, //
		{ TS::_P, TS::_M, TS::_P }, //
		{ TS::_P, TS::_P, TS::_P }, //
		{ TS::_P, TS::_P, TS::_P }, //
		{ TS::_M, TS::_P, TS::_P }, //
		{ TS::_M, TS::_P, TS::_P }, //
		{ TS::_M, TS::_M, TS::_P }, //

		{ TS::_M, TS::_M, TS::_M }, //
		{ TS::_M, TS::_M, TS::_P }, //
		{ TS::_P, TS::_M, TS::_M }, //
		{ TS::_P, TS::_M, TS::_P }, //
		{ TS::_P, TS::_P, TS::_M }, //
		{ TS::_P, TS::_P, TS::_P }, //
		{ TS::_M, TS::_P, TS::_M }, //
		{ TS::_M, TS::_P, TS::_P }, //
		};//

spPA_scatter3d AABBox(const AABBoxd3& box) {
	Listd lx;
	Listd ly;
	Listd lz;
	for (int i = 0; i < 24; i++) {
		lx.push_back(box.get(TS::_X, _ORDER[i][0]));
		ly.push_back(box.get(TS::_Y, _ORDER[i][1]));
		lz.push_back(box.get(TS::_Z, _ORDER[i][2]));
	}
	spPA_scatter3d res = spPA_scatter3d(
			new Plotly_actor_scatter3d(lx, ly, lz, 2));
	res->set_mode("lines");
	return res;

}

spPA_scatter3d AABBox(const std::vector<AABBoxd3>& vbox) {
	Listd lx;
	Listd ly;
	Listd lz;
	for (auto iter = vbox.begin(); iter != vbox.end(); ++iter) {
		const AABBoxd3& box = (*iter);
		for (int i = 0; i < 24; i++) {
			lx.push_back(box.get(TS::_X, _ORDER[i][0]));
			ly.push_back(box.get(TS::_Y, _ORDER[i][1]));
			lz.push_back(box.get(TS::_Z, _ORDER[i][2]));
		}
	}
	spPA_scatter3d res = spPA_scatter3d(
			new Plotly_actor_scatter3d(lx, ly, lz, 2));
	res->set_mode("lines");
	return res;
}

spPA_scatter3d BBTree(const BBTreed3& tree, BBTreed3::size_type height) {
	Listd lx;
	Listd ly;
	Listd lz;

	BBTreed3::Fun fun = [&lx, &ly, &lz, &height](BBTreed3::pNode pn) {
		if (pn->height() == height) {
			for (int i = 0; i < 24; i++) {
				lx.push_back(pn->box().get(TS::_X, _ORDER[i][0]));
				ly.push_back(pn->box().get(TS::_Y, _ORDER[i][1]));;
				lz.push_back(pn->box().get(TS::_Z, _ORDER[i][2]));
			}
		}
	};
	tree.InOrder(fun);

	spPA_scatter3d res = spPA_scatter3d(
			new Plotly_actor_scatter3d(lx, ly, lz, 2));
	res->set_mode("lines");
	return res;
}

spPA_scatter3d BBTreeLevel(const BBTreed3& tree, BBTreed3::size_type level) {
	Listd lx;
	Listd ly;
	Listd lz;

	BBTreed3::Fun fun = [&lx, &ly, &lz, &level](BBTreed3::pNode pn) {
		if (pn->level() == level) {
			for (int i = 0; i < 24; i++) {
				lx.push_back(pn->box().get(TS::_X, _ORDER[i][0]));
				ly.push_back(pn->box().get(TS::_Y, _ORDER[i][1]));
				lz.push_back(pn->box().get(TS::_Z, _ORDER[i][2]));
			}
		}
	};
	tree.InOrder(fun);

	spPA_scatter3d res = spPA_scatter3d(
			new Plotly_actor_scatter3d(lx, ly, lz, 2));
	res->set_mode("lines");
	return res;
}

}

}
