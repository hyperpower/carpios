#include "structure/s_grid.hpp"
#include "structure/s_ns_kim.hpp"
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

typedef NS_kim_<DIM> NS;

int main(int argc, char** argv){
    std::cout<< "--------- ns cavity ---------- \n";
    std::string strmesh(argv[1]);
    double tn  = std::strtod(argv[1],nullptr);
    double n = 32;
    int nx = int(n);
    int ny = int(n);
    std::cout<< "nx = " << nx << " ny = " << ny << std::endl;

    Poi pmin(0.0, 0.0, 0.0);
    Poi pmax(carpio::PI, carpio::PI, 0.0);
    Index mn(nx, ny, 1);
    spGrid spgrid(new Grid(pmin, pmax, mn, 2));
    NS ns(spgrid);

    //ns.set_advection_scheme();

    // set time
    Vt te = 2;
    ns.set_time(te * tn, 1/tn);

    NS::BC::Fun exact_u = [](Vt t, Vt x, Vt y, Vt z) {
        //-math.cos(x) * math.sin(y) * math.exp(-2.0 * t)
        return -std::cos(x) * std::sin(y) * std::exp(-2.0 * t);
    };

    NS::BC::Fun exact_v = [](Vt t, Vt x, Vt y, Vt z) {
        return std::sin(x) * std::cos(y) * std::exp(-2.0 * t);
    };

    NS::BC::Fun exact_p = [](Vt t, Vt x, Vt y, Vt z) {
        return -0.25 * (std::cos(2 * x) + std::cos(2 * y)) * std::exp(-4.0 * t);
    };



    NS::spBC bc0(new NS::BC()), bc1(new NS::BC());
    bc0->set_default_1_bc(exact_u);
    bc1->set_default_1_bc(exact_v);
    ns.add_bc(0, 0, "u", bc0);
    ns.add_bc(0, 1, "u", bc0);
    ns.add_bc(0, 2, "u", bc0);
    ns.add_bc(0, 3, "u", bc0);
    ns.add_bc(0, 0, "v", bc1);
    ns.add_bc(0, 1, "v", bc1);
    ns.add_bc(0, 2, "v", bc1);
    ns.add_bc(0, 3, "v", bc1);

    //typename NS::Function fun = [](Vt t, Vt x, Vt y, Vt z) {
    //  return  x*x + y*y - 4;
    //};
    ns.set_uniform_rho(1);
    //ns.set_diffusion_number(0.25, Re);
    

    Vt outstep = 10;
    std::string path = "./res_" + strmesh;
    ns.add_CS("erru", 0);
    ns.add_CS("errv", 0);
    ns.add_CS("errp", 0);
    ns.set_center_scalar("erru", exact_u, 0, -1, outstep, NS::Event::BEFORE);
    ns.set_center_scalar("errv", exact_v, 0, -1, outstep, NS::Event::BEFORE);
    ns.set_center_scalar("errp", exact_p, 0, -1, outstep, NS::Event::BEFORE);
    ns.set_output_error("u",   "erru",    0, -1, outstep, NS::Event::AFTER, path + "/erru");
    ns.set_output_error("v",   "errv",    0, -1, outstep, NS::Event::AFTER, path + "/errv");
    ns.set_output_error("p",   "errp",    0, -1, outstep, NS::Event::AFTER, path + "/errp");

    ns.set_output_time(0, -1, outstep,
                NS::Event::START |NS::Event::AFTER);
    
    ns.set_output_cs2file("u", path + "/u", 0, -1, outstep,
              NS::Event::START |NS::Event::AFTER);
    ns.set_output_cs2file("v", path + "/v", 0, -1, outstep,
              NS::Event::START |NS::Event::AFTER);
    ns.set_output_cs2file("p", path + "/p", 0, -1, outstep,
              NS::Event::START |NS::Event::AFTER);


    //ns.show_events();
    ns.set_projection_solver("IC_CGS", 100, 1e-6);

    ns.run();
    //ns.set_CS("v", fun);
    //ns.apply_bc("v");

    std::cout<< "--------- end ---------- \n";
}


