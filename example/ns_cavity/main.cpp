#include "structure/s_grid.hpp"
#include "structure/s_ns_explicit.hpp"
#include "structure/s_ns_kim.hpp"

#include "structure/s_io_plotly.hpp"
#include "structure/s_io_file.hpp"
#include <math.h>
#include <memory>

using namespace structure;

const short DIM = 2;
typedef Poi_<DIM> Poi;
typedef Grid_<DIM> Grid;
typedef std::shared_ptr<Grid_<DIM> > spGrid;
typedef Index_<DIM> Index;

typedef NS_explicit_<DIM> NS;

int main(){
    std::cout<< "--------- ns cavity ---------- \n";
    int nx = 30;
    int ny = 30;
    std::cout<< "nx = " << nx << " ny = " << ny << std::endl;

    Poi pmin(-0.5, -0.5, 0.0);
    Poi pmax(0.5, 0.5, 0.0);
    Index mn(nx, ny, 1);
    Grid grid(pmin, pmax, mn, 1);
    spGrid spgrid(new Grid(pmin, pmax, mn, 1));
    NS ns(spgrid);

    // set time
    ns.set_time(1000, 1e-1);


    NS::spBC bc0(new NS::BC()), bc1(new NS::BC());
    bc0->set_default_1_bc(0);
    bc1->set_default_1_bc(1);
    ns.add_bc(0, 0, "u", bc0);
    ns.add_bc(0, 1, "u", bc0);
    ns.add_bc(0, 2, "u", bc0);
    ns.add_bc(0, 3, "u", bc1);
    ns.add_bc(0, 0, "v", bc0);
    ns.add_bc(0, 1, "v", bc0);
    ns.add_bc(0, 2, "v", bc0);
    ns.add_bc(0, 3, "v", bc0);

    //typename NS::Function fun = [](Vt t, Vt x, Vt y, Vt z) {
    //  return  x*x + y*y - 4;
    //};
    Vt Re = 1000;
    ns.set_uniform_rho(Re);
    ns.set_diffusion_number(0.5, Re);

    ns.set_output_time(0, -1, 5,
                NS::Event::START |NS::Event::AFTER);

    std::string filename = "./result/u";
    ns.set_output_cs2file("u", filename, 0, -1, 50,
              NS::Event::START |NS::Event::AFTER);
    filename = "./result/v";
    ns.set_output_cs2file("v", filename, 0, -1, 50,
              NS::Event::START |NS::Event::AFTER);
    filename = "./result/p";
    ns.set_output_cs2file("p", filename, 0, -1, 50,
              NS::Event::START |NS::Event::AFTER);


    //ns.show_events();
    ns.set_projection_solver("IC_CGS", 100, 1e-6);

    ns.run();
    //ns.set_CS("v", fun);
    //ns.apply_bc("v");

    std::cout<< "--------- end ---------- \n";
}


