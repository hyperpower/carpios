#include "structure/s_grid.hpp"
#include "structure/s_poisson.hpp"
#include "structure/s_advection.hpp"
#include "structure/s_io_plotly.hpp"
#include "structure/s_operation.hpp"
#include "structure/s_io_file.hpp"
#include "utility/random.h"
#include <math.h>
#include <memory>

using namespace structure; 

const short DIM = 1;
typedef Poi_<DIM> Poi;
typedef Grid_<DIM> Grid;
typedef std::shared_ptr<Grid_<DIM> > spGrid;
typedef Index_<DIM> Index;
typedef Advection_<DIM> Eq;

Arr random_arr(Vt len, Arr::size_type n) {
    Arr arr(n);
    Vt dx = len / double(n);
    for (Arr::size_type i = 0; i < arr.size(); i++) {
        arr[i] = dx;
    }
    for (Arr::size_type i = 0; i < arr.size() - 1; i++) {
        Vt d2 = arr[i] + arr[i + 1];
        Vt dx1 = carpio::Random::nextFloat(0.7 * d2 / 2.0, 1.3 * d2 / 2.0);
        //Vt dx1 = d2 / 2. - 0.5;
        arr[i] = dx1;
        arr[i + 1] = d2 - dx1;
        ASSERT(arr[i] + arr[i + 1] - d2 < 1e-6);
    }
    Vt nlen =0;
    for (Arr::size_type i = 0; i < arr.size(); i++) {
        nlen += arr[i];
    }
    ASSERT(nlen - len < 1e-6);
    return arr;
}

int main_non(int argc, char** argv) {
    // main sheme
    std::string scheme(argv[1]);
    std::cout<< "Non uniform , " << scheme << std::endl;
    Vt dt = 0.1;
    Vt v  = 1;
    int step = 500;
    Vt xs  = 50;
    Vt xe  = xs + step * dt * v;
   


    Poi pmin(0.0, 0.0, 0.0);
    Arr arrx = random_arr(200, 200);
    spGrid spgrid(new Grid(pmin, 2, arrx));
    Eq eq(spgrid);

    // set time
    eq.set_time(step, dt);

    eq.initial_CS("u", 1);

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

    int outstep = 10;
    //eq.set_output_time(0, -1, 1, Eq::Event::START | Eq::Event::AFTER);
    eq.set_output_error("exact", "phi", 0, -1, outstep, 
                        Eq::Event::START | Eq::Event::AFTER,
                        "./result/nerror_"+ scheme);
    eq.set_center_scalar("exact", fun_exact, 0, -1, outstep, Eq::Event::START | Eq::Event::BEFORE);
    std::string filename = "./result/exact";
    eq.set_output_cs2file("exact", filename, 0, -1, outstep, Eq::Event::START | Eq::Event::AFTER);
    filename = "./result/nphi_" + scheme;
    eq.set_output_cs2file("phi", filename, 0, -1, outstep, Eq::Event::START | Eq::Event::AFTER);
    std::cout<<"scheme "<< scheme << std::endl;
    eq.set_advection_scheme(scheme);
    eq.run();
}

int main_uni(int argc, char** argv) {
    // main sheme
    std::string scheme(argv[1]);
    Vt dt = 0.1;
    Vt v  = 1;
    int step = 500;
    Vt xs  = 50;
    Vt xe  = xs + step * dt * v;

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

    int outstep = 10;
    //eq.set_output_time(0, -1, 1, Eq::Event::START | Eq::Event::AFTER);
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

int main(int argc, char** argv) {
    //
    main_non(argc, argv);
    main_uni(argc, argv);
}
