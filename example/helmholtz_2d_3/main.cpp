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

// On a class of preconditioners for solving
// the Helmholtz equation
// Y.A. Erlangga et al. / Applied Numerical Mathematics 50 (2004) 409–425


int main() {
    int n = 40;
    Vt k  = 2;
    Poi pmin(0.0, 0.0);
    Poi pmax(1.0, 1.0);
    Index mn(n,n);
    spGrid grid = spGrid(new Grid(pmin, pmax, mn, 3));

    Poisson poisson(grid);
    
    Poisson::BC::Fun exact = [](Vt t, Vt x, Vt y, Vt z) {
        return std::sin(carpio::PI * x) * std::sin(2 * carpio::PI * y);
    };

    Poisson::Function source = [k](Vt t, Vt x, Vt y, Vt z) {
        return (k * k - 5 * carpio::PI * carpio::PI) 
                   * std::sin(carpio::PI * x)
                   * std::sin(2 * carpio::PI * y);
    };
    poisson.set_source(source);
    poisson.set_uniform_alpha( k * k);

    Poisson::spBC bc0(new Poisson::BC());
    bc0->set_default_1_bc(0);
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
