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
    Index mn(30,  30);
    spGrid grid = spGrid(new Grid(pmin, pmax, mn, 3));

    Poisson poisson(grid);
    // set boundary condition

    Poisson::spBC bc1(new Poisson::BC()), bc2(new Poisson::BC()), bc3(new Poisson::BC());
    Poisson::BC::Fun fun_y2 = [](Vt t, Vt x, Vt y, Vt z) {
        return y * y;
    };
    Poisson::BC::Fun fun_x3 = [](Vt t, Vt x, Vt y, Vt z) {
        return x * x * x;
    };
    bc1->set_default_1_bc(1);
    bc2->set_default_1_bc(fun_y2);
    bc3->set_default_1_bc(fun_x3);
    poisson.add_bc(0, 0, "phi", bc2);
    poisson.add_bc(0, 1, "phi", bc1);
    poisson.add_bc(0, 2, "phi", bc3);
    poisson.add_bc(0, 3, "phi", bc1);

    Poisson::Function fun_source = [](Vt t, Vt x, Vt y, Vt z) {
        return 20 * std::cos(3 * carpio::PI * x) * std::sin(2 * carpio::PI * y);
    };

    poisson.set_source(fun_source);
    poisson.set_output_solver_residual("residual");
    poisson.set_solver("IC_CGS", 0.0, 1000, 1e-7);

    poisson.run();

    //output
    auto spphi = poisson.get_CS("phi");
    int res = Output("center_phi", *spphi);
}
