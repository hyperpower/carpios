#include "structure/s_grid.hpp"
#include "structure/s_poisson.hpp"
#include "structure/s_operation.hpp"
#include "structure/s_io_file.hpp"
#include <math.h>
#include <iostream>
#include <memory>

using namespace structure;


const short DIM = 1;
typedef Poi_<DIM> Poi;
typedef Grid_<DIM> Grid;
typedef std::shared_ptr<Grid_<DIM> > spGrid;
typedef Scalar_<DIM> Scalar;
typedef std::shared_ptr<Scalar_<DIM> > spScalar;
typedef Index_<DIM> Index;
typedef Operation_<DIM> Op;

int main() {
    Poi pmin(0.0);
    Poi pmax(3.0);
    Index mn(20);
    spGrid spg(new Grid(pmin, pmax, mn, 3));
    spScalar sps(new Scalar(spg));

    // function of initialize variable on the scalar field
    typename Op::Function fun = [](Vt t, Vt x, Vt y, Vt z) {
        return sin( cos(std::pow(x,2)+std::pow(sin(x/2.0),2))) - cos(std::pow(x,3))*sin(std::pow(x,2));
    };

    Op::Set(sps, fun);
    Vt x = 2.3;

    Poi p(x);

    std::cout << " localtion = " << x << std::endl;

    Vt val = Op::InterpolateCoordinate(sps, p[0]);
    std::cout << " res = " << val <<std::endl;
    Vt ext = fun(0.0, p[0], 0.0, 0.0);
    std::cout << " ext = " << ext <<std::endl;

    Output_Scalar("center_phi", *sps);
    Output_PointData("result", val, x, 0.0, 0.0);
    Output_PointData("exact",  ext, x, 0.0, 0.0);
}
