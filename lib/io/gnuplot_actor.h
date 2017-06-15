#ifndef _GNUPLOT_ACTOR_H_
#define _GNUPLOT_ACTOR_H_

#include "domain/cell.hpp"
#include "domain/grid.hpp"
#include "domain/boundary.hpp"
#include "domain/shape.hpp"
#include "gnuplot.h"
#include <memory>
#include "utility/clock.h"

namespace carpio {

namespace GnuplotActor {

typedef Float Cvt;
typedef Float Vt;

typedef std::shared_ptr<Gnuplot_actor> spActor;
typedef std::list<spActor> list_spActor;
/*
 * generate a spActor
 *
 * x : number of tests
 * y : wall time
 */
spActor Clock_Wall(const carpio::Clock& c);

// work with cell -------------------------------
typedef Cell_<Cvt, 1> Cell_1D;
typedef Cell_<Cvt, 2> Cell_2D;
typedef Cell_<Cvt, 3> Cell_3D;

spActor Cell(const Cell_2D& c);

// work with node
typedef Node_<Cvt, Vt, 1> Node_1D;
typedef Node_<Cvt, Vt, 2> Node_2D;
typedef Node_<Cvt, Vt, 3> Node_3D;

spActor Node(const Node_2D& node);

// ----------------------------------------------

// work with gird -------------------------------
typedef Grid_<Cvt, Vt, 1> Grid_1D;
typedef Grid_<Cvt, Vt, 2> Grid_2D;
typedef Grid_<Cvt, Vt, 3> Grid_3D;

typedef Ghost_<Cvt, Vt, 2> Ghost_2D;
typedef Ghost_<Cvt, Vt, 3> Ghost_3D;

spActor LeafNodes(const Grid_2D&);
spActor RootNodes(const Grid_2D&);
spActor LeafNodesContour(const Grid_2D&, St idx);
spActor GhostNodes(const Ghost_2D& g);
spActor GhostNodesContours(const Ghost_2D& g, St vi);
spActor GhostNodesContour_BoundaryIndex(const Ghost_2D& g);


// work with shape ------------------------------
typedef Shape_<Cvt, 2> Shape_2D;

spActor Shape(const Shape_2D& shape, St base_idx);

}

}

#endif
