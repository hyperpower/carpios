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

    Vt v  = 1;
    Vt l  = 100;
    St n  = 50;
    Vt dx = l / n;
    Vt CFL= 0.25;
    Vt dt = CFL * dx / v;
    Vt dis = 50; 
    Vt t   = std::sqrt( 2* dis * dis)  / v;
    int step = t / dt;
    std::cout << "step = " << step <<std::endl;
    // step = 1;
    // CFL = v * dt / dx
    //

    Poi pmin(0.0, 0.0, 0.0);
    Poi pmax(l, l, 0.0);
    Index mn(n, n, 2);
    spGrid spgrid(new Grid(pmin, pmax, mn, 2));
    Eq eq(spgrid);

    // set time
    eq.set_time(step, dt);

    Vt ux, uy;
    ux = std::sqrt(v * v / 2);
    uy = std::sqrt(v * v / 2);
    Vt ss = 25;
    eq.initial_CS("u", ux);
    eq.initial_CS("v", uy);

    Eq::Function fun_initial_phi = [ss](Vt t , Vt x, Vt y, Vt z) {
        return ( (y - ss) * ( y - ss) + (x - ss) * ( x - ss) < 100) ? 1: 0;
    };
    
    eq.initial_CS("phi", fun_initial_phi);


    Eq::spBC bc0(new Eq::BC());
    bc0->set_default_1_bc(0);
    eq.add_bc(0, 0, "phi", bc0);
    eq.add_bc(0, 1, "phi", bc0);
    eq.add_bc(0, 2, "phi", bc0);
    eq.add_bc(0, 3, "phi", bc0);

    Eq::Function fun_exact = [ss, ux, uy](Vt t , Vt x, Vt y, Vt z) {
        Vt nx = ss + ux * t;
        Vt ny = ss + uy * t;
        return ( (y - ny) * ( y - ny) + (x - nx) * ( x - nx) < 100) ? 1: 0;
    };

    eq.add_CS("exact", fun_exact);

    Vt outnum  = 30;
    int outstep = std::round(step / outnum);
    eq.set_output_time(0, -1, outstep, Eq::Event::START | Eq::Event::AFTER);
    eq.set_output_error("exact", "phi", 0, -1, outstep, Eq::Event::START | Eq::Event::AFTER, "./result/error_"+ scheme);
    eq.set_center_scalar("exact", fun_exact, 0, -1, outstep, Eq::Event::START | Eq::Event::BEFORE);
    std::string filename = "./result/exact";
    eq.set_output_cs2file("exact", filename, 0, -1, outstep, Eq::Event::START | Eq::Event::AFTER);
    filename = "./result/phi_" + scheme;
    eq.set_output_cs2file("phi", filename, 0, -1, outstep, Eq::Event::START | Eq::Event::AFTER);
    std::cout<<"scheme "<< scheme << std::endl;
    eq.set_advection_scheme(scheme);
    eq.run();
}
