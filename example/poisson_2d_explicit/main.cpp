#include "structure/s_grid.hpp"
#include "structure/s_poisson.hpp"
#include "structure/s_operation.hpp"
#include "structure/s_io_file.hpp"
#include <math.h>
#include <iostream>
#include <memory>

using namespace structure;


const short DIM = 2;
typedef Poi_<DIM> Poi;
typedef Grid_<DIM> Grid;
typedef std::shared_ptr<Grid> spGrid;
typedef Index_<DIM> Index;
typedef Poisson_<DIM> Poisson;

int main() {
    Poi pmin(0.0, 0.0);
    Poi pmax(1.0, 1.0);
    Index mn(30, 30);
    spGrid grid = spGrid(new Grid(pmin, pmax, mn, 3));

    Poisson poisson(grid);
    // set boundary condition
    Poisson::spBC bc1(new Poisson::BC()), bc2(new Poisson::BC());
    bc1->set_default_1_bc(10);
    bc2->set_default_1_bc(1);
    poisson.add_bc(0, 0, "phi", bc1);
    poisson.add_bc(0, 1, "phi", bc2);
    poisson.add_bc(0, 2, "phi", bc2);
    poisson.add_bc(0, 3, "phi", bc2);

    poisson.set_time(10000, 0.00001);
    poisson.set_output_time(0, -1, 10,
            Poisson::Event::START | Poisson::Event::AFTER);

   

    //output
    std::string filename = "./result/phi";
    poisson.set_output_cs2file("phi", filename, 0, -1, 100,
              Poisson::Event::START |Poisson::Event::AFTER);

    poisson.run();
}
