#ifndef _PLOTLY_ACTOR_H_
#define _PLOTLY_ACTOR_H_

#include <io/plotly.h>
#include "algebra/array_list.hpp"
#include "geometry/geometry.hpp"

#include "ts/ts.h"

#include "plotly.h"
#include <Python.h>
#include <map>
#include <memory>

namespace carpio {


namespace PlotlyActor{

typedef std::shared_ptr<Plotly_actor> spPA;
typedef std::shared_ptr<Plotly_actor_scatter> spPA_scatter;
typedef std::shared_ptr<Plotly_actor_scatter3d> spPA_scatter3d;
typedef std::shared_ptr<Plotly_actor_mesh3d> spPA_mesh3d;

// working with Arraylist
typedef ArrayListV<double> ArrListd;
typedef std::list<double> Listd;

spPA_scatter XY(const ArrListd& x, const ArrListd& y);
spPA_scatter XY(const std::list<double>& x, const std::list<double>& y);

//
typedef Point_<double, 3> Poid3;
spPA_scatter3d Segment(const Poid3& s, const Poid3& e);

// working with geometry/triangle
typedef Triangle_<double, 3> Trid3;
spPA_mesh3d Triangle(const Trid3&);


// working with ts
typedef TS::Surface<double, 3> Surd3;
typedef TS::Surface<double, 2> Surd2;
typedef TS::Triangle<double, 2> Trid2;
typedef TS::Point<double, 2> Poid2;
typedef std::list<std::shared_ptr<const Poid2> > list_spcPoid2;
typedef TS::Face<double, 3> Faced3;
spPA_mesh3d Surface(const Surd3&);
spPA_scatter SurfaceWireFrame(const Surd2&);
spPA_mesh3d Face(const Faced3&);
spPA_scatter3d SurfaceWireFrame(const Surd3&);
spPA_scatter3d SurfaceWireFrame(const Trid2&);
spPA_scatter ScatterPoints(const list_spcPoid2&);
// working with AABBox
typedef TS::AABBox<double, 3> AABBoxd3;
typedef TS::BBTree<AABBoxd3> BBTreed3;
spPA_scatter3d AABBox(const AABBoxd3& box);
spPA_scatter3d AABBox(const std::vector<AABBoxd3>& vbox);
spPA_mesh3d Triangle(const std::vector<AABBoxd3>&);
spPA_scatter3d TriangleNormal(const std::vector<AABBoxd3>&);
spPA_scatter3d TriangleNormal(const Surd3&, double scale = 1.0);
spPA_scatter3d TriangleNormal(const Faced3&);
spPA_scatter3d BBTree(const BBTreed3& tree, BBTreed3::size_type hight);
spPA_scatter3d BBTreeLevel(const BBTreed3& tree, BBTreed3::size_type level);

}

}

#endif
