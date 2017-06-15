#include "structure/s_grid.hpp"
#include "structure/s_poisson.hpp"
#include "structure/s_operation.hpp"
#include "structure/s_io_file.hpp"
#include <math.h>
#include <iostream>
#include <memory>

using namespace structure;

// AN EFFICIENT QUARTER-SWEEP MODIFIED SOR
// ITERATIVE METHOD FOR SOLVING
// HELMHOLTZ EQUATION
const short DIM = 2;
typedef Poi_<DIM> Poi;
typedef Grid_<DIM> Grid;
typedef std::shared_ptr<Grid> spGrid;
typedef Index_<DIM> Index;
typedef Poisson_<DIM> Poisson;


int main() {
    Poi pmin(0.0, 0.0);
    Poi pmax(1.0, 2 * M_PI);
    Index mn(50, 2 * M_PI * 50);
    spGrid grid = spGrid(new Grid(pmin, pmax, mn, 3));

    Poisson poisson(grid);
    
    Poisson::BC::Fun exact = [](Vt t, Vt x, Vt y, Vt z) {
        return x * x * x * x * std::sin(y);
    };

    Poisson::Function source = [](Vt t, Vt x, Vt y, Vt z) {
        return  std::sin(y) * ( 12 * x * x +  3 * x * x * x * x);
    };
    poisson.set_source(source);
    poisson.set_uniform_alpha(-2);

    Poisson::spBC bc0(new Poisson::BC());
    Poisson::BC::Fun fbc0 = [](Vt t, Vt x, Vt y, Vt z) {
        return x * x *x * x* std::sin(y);
    };
    bc0->set_default_1_bc(fbc0);
    poisson.add_bc(0, 0, "phi", bc0);
    poisson.add_bc(0, 1, "phi", bc0);
    poisson.add_bc(0, 2, "phi", bc0);
    poisson.add_bc(0, 3, "phi", bc0);

    poisson.add_CS("exact", exact);

    //poisson.set_output_solver_residual("residual");
    poisson.set_solver("IC_CGS", 0.0, 1000, 1e-9);

    poisson.run();

    auto spphi = poisson.get_CS("phi");
    Output_Scalar("center_phi", *spphi);
    auto spexa = poisson.get_CS("exact");
    Output_Scalar("center_exact", *spexa);
}
