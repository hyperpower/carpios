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
typedef std::shared_ptr<Grid> spGrid;
typedef Index_<DIM> Index;
typedef Poisson_<DIM> Poisson;

int main() {
    Poi pmin(1.0);
    Poi pmax(5.0);
    Index mn(20);
    spGrid grid = spGrid(new Grid(pmin, pmax, mn, 3));

    Poisson poisson(grid);
    // set boundary condition
    Poisson::spBC bc1(new Poisson::BC()), bc2(new Poisson::BC());
    bc1->set_default_1_bc(10);
    bc2->set_default_1_bc(1);
    poisson.add_bc(0, 0, "phi", bc1);
    poisson.add_bc(0, 1, "phi", bc2);

    poisson.run();

    //output
    auto spphi = poisson.get_CS("phi");
    //int res  = Output("center_phi", *spphi);
    int res2 = Output_Scalar("center_phi", *spphi);
}
