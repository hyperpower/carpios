#ifndef _TEST_HPP_
#define _TEST_HPP_

#include <iostream>
using namespace std;
#include <gtest/gtest.h>
#include <glog/logging.h>

//#include "../test/utility/_test_any.hpp"
//#include "../test/algebra/_test_matrix.hpp"
//#include "../test/algebra/_test_matrix_small.hpp"
//#include "../test/algebra/_test_polynomial.hpp"
//#include "../test/algebra/_test_solver.hpp"

//#include "../test/domain/_test_grid_2D.hpp"
//#include "../test/domain/_test_domain.hpp"

//#include "../test/io/_test_csv.hpp"

//#include "../test/calculation/_test_equation.hpp"
//#include "../test/calculation/_test_poisson.hpp"
//#include "../test/calculation/_test_ns_explicit.hpp"
//#include "../test/calculation/_test_interpolate.hpp"

#include "../test/ts/test_surface_constructor_3d.hpp"
//#include "../test/ts/test_AABBox.h"
//#include "../test/ts/test_overlap.h"
//#include "../test/ts/test_delaunay.h"

//#include "../test/structure/_test_s_poisson.hpp"
//#include "../test/structure/_test_s_advection.hpp"
//#include "../test/structure/_test_s_operation.hpp"
//#include "../test/structure/_test_s_ns.hpp"
//#include "../test/structure/_test_s_grid.hpp"


//#include "../test/geometry/test_triangle.hpp"
//#include "../test/geometry/test_polygon.hpp"
//#include "../test/geometry/test_contour.hpp"
//#include "../test/io/_test_plotly.hpp"

//#include "../test/parallel/test_mpi.hpp"

namespace carpio {

int RunTests(int argc, char **argv) {
	//	::testing::GTEST_FLAG(filter) = "Poisson.adpgrid_source";
	//::testing::GTEST_FLAG(filter) = "Grid2D.no_leaf";
	::testing::InitGoogleTest(&argc, argv);
	//FLAGS_log_dir = "./";
	FLAGS_logtostderr = true;
	FLAGS_stderrthreshold = 0;
	google::InitGoogleLogging(argv[0]);
	int res = RUN_ALL_TESTS();
	return res;
}

}

#endif
