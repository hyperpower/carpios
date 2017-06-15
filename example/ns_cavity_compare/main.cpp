#include "structure/s_grid.hpp"
#include "structure/s_ns_explicit.hpp"
#include "structure/s_io_plotly.hpp"
#include "structure/s_io_file.hpp"
#include <math.h>
#include <memory>
#include <stdlib.h>

using namespace structure;

const short DIM = 2;
typedef Poi_<DIM> Poi;
typedef Grid_<DIM> Grid;
typedef std::shared_ptr<Grid_<DIM> > spGrid;
typedef Index_<DIM> Index;

typedef NS_explicit_<DIM> NS;

int main(int argc, char** argv){
    std::cout<< "--------- ns cavity ---------- \n";
    std::string strre(argv[1]);
    std::string strmesh(argv[2]);
    Vt Re = std::strtod(argv[1],nullptr);
    double n  = std::strtod(argv[2],nullptr);
    std::string scheme(argv[3]);
    std::cout<< "Re =      " << Re << std::endl;
    std::cout<< "uniform , " << scheme << std::endl;
    int nx = int(n);
    int ny = int(n);
    std::cout<< "nx = " << nx << " ny = " << ny << std::endl;

    Poi pmin(-0.5, -0.5, 0.0);
    Poi pmax(0.5, 0.5, 0.0);
    Index mn(nx, ny, 1);
    spGrid spgrid(new Grid(pmin, pmax, mn, 2));
    NS ns(spgrid);

    ns.set_advection_scheme(scheme);

    // set time
    ns.set_time(100000, 1e-1);


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
    ns.set_uniform_rho(Re);
    ns.set_diffusion_number(0.25, Re);

    ns.set_output_time(0, -1, 100,
                NS::Event::START |NS::Event::AFTER);

    ns.set_stop_cs("u", 1e-6,1e-6,1e-6, 
                    0, -1, 100, NS::Event::AFTER);

    std::string path = "./res_"+strre+"_"+strmesh+"_"+scheme;
    ns.set_output_cs2file("u", path + "/u", 0, -1, 50,
              NS::Event::START |NS::Event::AFTER);
    ns.set_output_cs2file("v", path + "/v", 0, -1, 50,
              NS::Event::START |NS::Event::AFTER);
    ns.set_output_cs2file("p", path + "/p", 0, -1, 50,
              NS::Event::START |NS::Event::AFTER);


    //ns.show_events();
    ns.set_projection_solver("IC_CGS", 100, 1e-6);

    ns.run();
    //ns.set_CS("v", fun);
    //ns.apply_bc("v");

    std::cout<< "--------- end ---------- \n";
}


