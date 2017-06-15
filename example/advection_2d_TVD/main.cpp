#include "structure/s_grid.hpp"
#include "structure/s_poisson.hpp"
#include "structure/s_advection.hpp"
#include "structure/s_io_plotly.hpp"
#include "structure/s_operation.hpp"
#include "structure/s_io_file.hpp"
#include <math.h>
#include <memory>

using namespace structure; 

const short DIM = 2;
typedef Poi_<DIM> Poi;
typedef Grid_<DIM> Grid;
typedef std::shared_ptr<Grid_<DIM> > spGrid;
typedef Index_<DIM> Index;
typedef Advection_<DIM> Eq;

int main(int argc, char** argv) {
    // main sheme
    std::string scheme(argv[1]);
    Vt dt = 0.002;
    Vt v  = 1;
    int step = 1000;

    Poi pmin(0.0, 0.0, 0.0);
    Poi pmax(1, 1, 0.0);
    Index mn(50, 50, 2);
    spGrid spgrid(new Grid(pmin, pmax, mn, 2));
    Eq eq(spgrid);

    // set time
    eq.set_time(step, dt);

    eq.initial_CS("u", v);
    eq.initial_CS("v", v);
    
    eq.initial_CS("phi", 0.0);

    Eq::Function fun_bc_phi = [](Vt t , Vt x, Vt y, Vt z) {
        return (y >= 0.0 && y <= 0.3) ? 1 : 0;
    };

    Eq::spBC bc0(new Eq::BC()), bc1(new Eq::BC());
    bc0->set_default_1_bc(0);
    bc1->set_default_1_bc(fun_bc_phi);
    eq.add_bc(0, 0, "phi", bc1);
    eq.add_bc(0, 1, "phi", bc0);
    eq.add_bc(0, 2, "phi", bc0);
    eq.add_bc(0, 3, "phi", bc0);

    Eq::Function fun_exact = [](Vt t , Vt x, Vt y, Vt z) {
        return (x >= y - 0.3 && x <= y) ? 1 : 0;
    };

    eq.add_CS("exact", fun_exact);

    eq.set_output_time(0, -1, 10, Eq::Event::START | Eq::Event::AFTER);
    eq.set_output_error("exact", "phi", 0, -1, 10, Eq::Event::START | Eq::Event::AFTER, "./result/error_"+ scheme);
    eq.set_center_scalar("exact", fun_exact, 0, -1, 10, Eq::Event::START | Eq::Event::BEFORE);
    std::string filename = "./result/exact";
    eq.set_output_cs2file("exact", filename, 0, -1, 10, Eq::Event::START | Eq::Event::AFTER);
    filename = "./result/phi_" + scheme;
    eq.set_output_cs2file("phi", filename, 0, -1, 10, Eq::Event::START | Eq::Event::AFTER);
    std::cout<<"scheme "<< scheme << std::endl;
    eq.set_advection_scheme(scheme);
    eq.run();
}
