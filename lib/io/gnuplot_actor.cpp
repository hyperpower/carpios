#include "gnuplot_actor.h"

namespace carpio {

namespace GnuplotActor {

/*
 * generate a spActor
 *
 * x : number of tests
 * y : wall time
 */
spActor Clock_Wall(const carpio::Clock& c) {
	spActor actor = spActor(new Gnuplot_actor());
	actor->command() = "using 1:2 title \"" + c.name() + "\" ";
	actor->style() = "with linespoints lw 2";
	auto iter_t = c.begin_wall();
	auto iter_n = c.begin_num();
	int i = 0;
	for (; iter_t != c.end_wall();) {
		//                      "number of tests  time"
		tick_t dt;
		if (i == 0) {
			dt = Clock::TicksInBetween((*iter_t), c.start_time_cpu());
		} else {
			auto prev_t = std::prev(iter_t, 1);
			dt = Clock::TicksInBetween(*iter_t, *prev_t);
		}
		actor->data().push_back(
				ToString(*iter_n, Clock::TicksToSecondsD(dt), " "));
		iter_t++;
		iter_n++;
		i++;
	}
	return actor;
}

int DataPushBack(std::list<std::string>& ldata, const Cell_2D& c) {
	ldata.push_back(ToString(c.get(_M_, _X_), c.get(_M_, _Y_), " "));
	ldata.push_back(ToString(c.get(_P_, _X_), c.get(_M_, _Y_), " "));
	ldata.push_back("");
	ldata.push_back(ToString(c.get(_P_, _X_), c.get(_M_, _Y_), " "));
	ldata.push_back(ToString(c.get(_P_, _X_), c.get(_P_, _Y_), " "));
	ldata.push_back("");
	ldata.push_back(ToString(c.get(_P_, _X_), c.get(_P_, _Y_), " "));
	ldata.push_back(ToString(c.get(_M_, _X_), c.get(_P_, _Y_), " "));
	ldata.push_back("");
	ldata.push_back(ToString(c.get(_M_, _X_), c.get(_P_, _Y_), " "));
	ldata.push_back(ToString(c.get(_M_, _X_), c.get(_M_, _Y_), " "));
	ldata.push_back("");
	return _SUCCESS;
}

int DataPushBack_Contour(std::list<std::string>& ldata, const Node_2D* pn,
		St idx) {
	ldata.push_back(
			ToString(pn->cp(_X_), pn->cp(_Y_), pn->p(_M_, _X_), pn->p(_P_, _X_),
					pn->p(_M_, _Y_), pn->p(_P_, _Y_), pn->cda(idx), " ")); //point" "));
	return _SUCCESS;
}

spActor Cell(const Cell_2D& c) {
	spActor actor = spActor(new Gnuplot_actor());
	actor->clear();
	actor->command() = "using 1:2 title \"\" ";
	DataPushBack(actor->data(), c);
	return actor;
}

spActor Node(const Node_2D& node) {
	spActor actor = spActor(new Gnuplot_actor());
	actor->clear();
	actor->command() = "using 1:2 title \"\" ";
	DataPushBack(actor->data(), *(node.cell));
	return actor;
}

spActor LeafNodes(const Grid_2D& g) {
	spActor actor = spActor(new Gnuplot_actor());
	actor->command() = "using 1:2 title \"\" ";
	for (Grid_2D::const_iterator_leaf iter = g.begin_leaf();
			iter != g.end_leaf(); ++iter) {
		DataPushBack(actor->data(), *(iter->cell));
	}
	return actor;
}

spActor RootNodes(const Grid_2D& g) {
	spActor actor = spActor(new Gnuplot_actor());
	actor->command() = "using 1:2 title \"\" ";
	for (Grid_2D::const_iterator iter = g.begin(); iter != g.end(); ++iter) {
		Grid_2D::pNode proot = (*iter);
		if (proot != nullptr) {
			DataPushBack(actor->data(), *(proot->cell));
		}
	}
	return actor;
}

spActor LeafNodesContour(const Grid_2D& g, St idx) {
	spActor actor = spActor(new Gnuplot_actor());
	actor->command() = "using 1:2:3:4:5:6:7 title \"\" ";
	actor->style() = "with boxxy fs solid palette";
	for (Grid_2D::const_iterator_leaf iter = g.begin_leaf();
			iter != g.end_leaf(); ++iter) {
		const Node_2D* pn = iter.get_pointer();
		DataPushBack_Contour(actor->data(), pn, idx);
	}
	return actor;
}

spActor GhostNodes(const Ghost_2D& g) {
	spActor actor = spActor(new Gnuplot_actor());
	actor->command() = "using 1:2 title \"\" ";
	for (typename Ghost_2D::const_iterator iter = g.begin(); iter != g.end();
			++iter) {
		DataPushBack(actor->data(), *(iter->second.pghost->cell));
	}
	return actor;
}

spActor GhostNodesContour_BoundaryIndex(
		const Ghost_2D& g) {
	spActor actor = spActor(new Gnuplot_actor());
	actor->command() = "using 1:2:3:4:5:6:7 title \"\" ";
	actor->style() = "with boxxy fs solid palette";
	typedef typename Ghost_2D::GhostNode Node;
	std::function<void(const Node&)> fun = [&](const Node& n) {
		//
			typename Ghost_2D::GhostVal gval=n.second;
			typename Ghost_2D::pNode pn = gval.pghost;
			typename Ghost_2D::GhostID gid = n.first;
			//assume segments for each shape less than 10000
			Float v = gval.shape_idx * 10 + gval.seg_idx;
			//std::cout<< v << " " <<gval.seg_idx << " " << gval.shape_idx<<"\n";
			if (pn != nullptr) {
				actor->data().push_back(
						ToString(pn->cp(_X_), pn->cp(_Y_), pn->p(_M_, _X_), pn->p(_P_, _X_),
								pn->p(_M_, _Y_), pn->p(_P_, _Y_), v, " "));
			}
		};
	g.for_each_node(fun);

	return actor;
}

spActor GhostNodesContours(const Ghost_2D& g, St vi) {
	spActor actor = spActor(new Gnuplot_actor());
	actor->command() = "using 1:2:3:4:5:6:7 title \"\" ";
	actor->style() = "with boxxy fs solid palette";
	typedef typename Ghost_2D::GhostNode Node;
	std::function<void(const Node&)> fun = [&](const Node& n) {
		//
			typename Ghost_2D::GhostVal gval=n.second;
			typename Ghost_2D::pNode pn = gval.pghost;
			typename Ghost_2D::GhostID gid = n.first;
			//assume segments for each shape less than 10000
			Float v = pn->cda(vi);
			//std::cout<< v << " " <<gval.seg_idx << " " << gval.shape_idx<<"\n";
			if (pn != nullptr) {
				actor->data().push_back(
						ToString(pn->cp(_X_), pn->cp(_Y_), pn->p(_M_, _X_), pn->p(_P_, _X_),
								pn->p(_M_, _Y_), pn->p(_P_, _Y_), v, " "));
			}
		};
	g.for_each_node(fun);

	return actor;
}

spActor Shape(const Shape_2D& g, St base_idx) {
	spActor actor = spActor(new Gnuplot_actor());
	actor->command() = "using 1:2:3 title \"\" ";
	actor->style() = "with lines lc variable";
	if (g.empty()) {
		actor->data().push_back("");
		return actor;
	}
	typedef typename Shape2D::S2D::Point Poi;
	for (St i = 0; i < g.size_vertexs(); ++i) {
		const Poi& p = g.v(i);
		actor->data().push_back(ToString(p.x(), p.y(), i + base_idx, " "));
	}
	const Poi& pstart = g.v(0);
	actor->data().push_back(ToString(pstart.x(), pstart.y(), base_idx, " "));
	actor->data().push_back("");
	return actor;
}

}

}
