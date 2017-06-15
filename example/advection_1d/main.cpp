#include "structure/s_grid.hpp"
#include "structure/s_poisson.hpp"
#include "structure/s_advection.hpp"
#include "structure/s_io_plotly.hpp"
#include "structure/s_operation.hpp"
#include "structure/s_io_file.hpp"
#include <math.h>
#include <memory>

using namespace structure; 

const short DIM = 1;
typedef Poi_<DIM> Poi;
typedef Grid_<DIM> Grid;
typedef std::shared_ptr<Grid_<DIM> > spGrid;
typedef Index_<DIM> Index;
typedef Advection_<DIM> Eq;

int main(int argc, char** argv) {
    // main sheme
    std::string scheme(argv[1]);
    Vt dt = 0.1;
    Vt v  = 1;
    int step = 500;
    Vt xs  = 50;
    Vt xe  = xs + step * dt * v;
    std::cout<< "x end = " << xe << std::endl;


    Poi pmin(0.0, 0.0, 0.0);
    Poi pmax(200, 200, 0.0);
    Index mn(200, 200, 2);
    spGrid spgrid(new Grid(pmin, pmax, mn, 2));
    Eq eq(spgrid);

    // set time
    eq.set_time(step, dt);

    eq.initial_CS("u", 1);
    eq.initial_CS("v", 1);
    Eq::Function fun_init_phi = [xs](Vt t , Vt x, Vt y, Vt z) {
        return (x >= xs - 10 && x <= xs + 10) ? 1 : 0;
    };
    eq.initial_CS("phi", fun_init_phi);

    Eq::spBC bc0(new Eq::BC()), bc1(new Eq::BC());
    bc0->set_default_1_bc(0);
    bc1->set_default_1_bc(1);
    eq.add_bc(0, 0, "phi", bc0);
    eq.add_bc(0, 1, "phi", bc0);
    //eq.add_bc(0, 2, "phi", bc0);
    //eq.add_bc(0, 3, "phi", bc0);

    Eq::Function fun_exact = [xs,v](Vt t , Vt x, Vt y, Vt z) {
        Vt xe  = xs + t * v;
        return (x >= xe - 10 && x <= xe + 10) ? 1 : 0;
    };

    eq.add_CS("exact", fun_exact);

    //eq.set_output_time(0, -1, 1, Eq::Event::START | Eq::Event::AFTER);
    eq.set_output_error("exact", "phi", 0, -1, 1, Eq::Event::START | Eq::Event::AFTER, "./result/error_"+ scheme);
    eq.set_center_scalar("exact", fun_exact, 0, -1, 1, Eq::Event::START | Eq::Event::BEFORE);
    std::string filename = "./result/exact";
    eq.set_output_cs2file("exact", filename, 0, -1, 1, Eq::Event::START | Eq::Event::AFTER);
    filename = "./result/phi_" + scheme;
    eq.set_output_cs2file("phi", filename, 0, -1, 1, Eq::Event::START | Eq::Event::AFTER);
    std::cout<<"scheme "<< scheme << std::endl;
    eq.set_advection_scheme(scheme);
    eq.run();
}
