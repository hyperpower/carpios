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

// The Solution of 2D Helmholtz Equations by
// Modified Explicit Group Iterative Method
// Proc. Int. Conf. on Advances in Computing, Control, and Telecommunication Technologies

Vt exact(Vt x, Vt y){
   return 2* x * x + y * y; 
}

int main() {
    Poi pmin(0.0, 0.0);
    Poi pmax(1.0, 1.0);
    Index mn(40, 40);
    spGrid grid = spGrid(new Grid(pmin, pmax, mn, 3));

    Poisson poisson(grid);
    // set boundary condition
    Vt alpha = -10;
    
    Poisson::BC::Fun exact = [](Vt t, Vt x, Vt y, Vt z) {
        return 2 * x * x + y * y;
    };

    Poisson::Function source = [alpha](Vt t, Vt x, Vt y, Vt z) {
        return 6 - alpha *( 2 * x * x +  y * y);
    };
    poisson.set_source(source);
    poisson.set_uniform_alpha(alpha);

    Poisson::spBC bc0(new Poisson::BC());
    Poisson::BC::Fun fbc0 = [](Vt t, Vt x, Vt y, Vt z) {
        return 2 * x * x + y * y;
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
