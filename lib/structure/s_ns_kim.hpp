#ifndef _S_NS_KIM_HPP
#define _S_NS_KIM_HPP

#include "s_data.hpp"
#include "s_define.hpp"
#include "s_stencil.hpp"
#include "s_equation.hpp"
#include "s_solver.hpp"
#include "s_vector.hpp"
#include "s_io_plotly.hpp"
#include "s_ns.hpp"

namespace structure {

// Kim and Moin (JCP 59 (1985), 8-23)!

template<St DIM>
class NS_kim_: public NS_<DIM>
{
public:
    static const St Dim = DIM;
    typedef NS_<DIM> Base;
    typedef Equation_<DIM> Equation;
    typedef Grid_<DIM> Grid;
    typedef std::shared_ptr<Grid> spGrid;
    typedef Index_<Dim> Index;

    typedef Scalar_<DIM> Scalar;
    typedef std::shared_ptr<Scalar> spScalar;
    typedef VectorCenter_<DIM> VectorCenter;
    typedef std::shared_ptr<VectorCenter> spVectorCenter;
    typedef VectorFace_<DIM> VectorFace;
    typedef std::shared_ptr<VectorFace> spVectorFace;

    typedef Stencil_<DIM> Stencil;
    typedef std::shared_ptr<Stencil> spStencil;
    typedef typename Stencil::Expression Expression;
    typedef typename Stencil::spExpression spExpression;

    // function define ------------------------------
    typedef std::function<Vt(Vt, Vt, Vt, Vt)> Function;
    typedef std::unordered_map<std::string, Function> Functions;
    // Vaules ---------------------------------------
    typedef std::unordered_map<std::string, Vt> Values;

    // Variables ------------------------------------
    typedef std::map<std::string, spScalar> CenterScalars;
    // Event ----------------------------------------
    typedef Event_<DIM> Event;
    typedef std::shared_ptr<Event> spEvent;
    typedef std::unordered_map<std::string, spEvent> Events;
    typedef EventFlag_<Dim> Flag;

    typedef Ghost_<Dim> Ghost;
    typedef std::shared_ptr<Ghost> spGhost;
    typedef typename Ghost::Fun_index Fun_index;

    // time -----------------------------------------
    typedef Time_<DIM> Time;
    typedef std::shared_ptr<Time> spTime;
    // solver ---------------------------------------
    typedef Solver_<DIM> Solver;
    typedef std::shared_ptr<Solver_<DIM>> spSolver;
    typedef std::map<std::string, spSolver> Solvers;
    // Operation ------------------------------------
    typedef Poisson_<DIM> Poisson;
    typedef std::shared_ptr<Poisson_<DIM>> spPoisson;
    typedef Operation_<DIM> Operation;

    typedef carpio::Any Any;

protected:
    spVectorCenter _veomc;  // velocity previous step
    spVectorFace   _veomf;  // velocity previous step face
    spVectorCenter _veosc;  // veo star on center
    spVectorFace   _veosf;  // veo star on face

    spVectorCenter _sscr;    // save the source term
    spScalar _smu;           // the new mu

    spPoisson _predict;

    std::array<std::string, 3> _nvs;
    std::array<std::string, 3> _nvm;
    std::array<std::string, 3> _nvsf;
    std::array<std::string, 3> _nvmf;

public:
    NS_kim_(spGrid spg) :
        Base(spg) {
        this->_default_setup();
        //
        Function fun = [](Vt, Vt, Vt, Vt) {return 0;};
        this->_nvs  =  { "us",  "vs",  "ws"  };
        this->_nvm  =  { "um",  "vm",  "wm"  };
        this->_nvsf =  { "usf", "vsf", "wsf" };
        this->_nvmf =  { "umf", "vmf", "wmf" };
        FOR_EACH_DIM {
            this->_functions[_nvs[d]] = fun;
            this->_functions[_nvsf[d]]  = fun;
            this->_functions[_nvm[d]] = fun;
            this->_functions[_nvmf[d]] = fun;
        }
        //
        this->_predict = spPoisson(new Poisson(this->_grid));
        this->_functions["phi"] = fun;
    }

    int initial() {
        this->_ns_initial();
        // initial field veo star on center and face
        this->_veosc = this->_new_veoc("us",  "vs",  "ws");
        this->_veomc = this->_new_veoc("um",  "vm",  "wm");
        this->_sscr = this->_new_veoc("xscr",  "yscr",  "zscr");
        this->_veosf = this->_new_veof("usf", "vsf", "wsf");
        this->_veomf = this->_new_veof("umf", "vmf", "wmf");

        spGrid spg = this->_grid;
        // phi -------------------------------------
        spScalar spphi = spScalar(new Scalar(spg));
        this->_css["phi"] = spphi;
        this->set_CS("phi", this->_functions["phi"]);
        // copy boundary condition ----------------
        FOR_EACH_DIM {
            this->_bi->copy(this->_nv[d], this->_nvs[d]);
        }
        return 1;
    }

    int run_one_step(St step) {
        // debug
        int stepo = 0;
        // copy velocity
        Operation::Copy(this->_veoc, this->_veomc);
        Operation::Copy(this->_veof, this->_veomf);

        // restrict time
        //
        //this->restrict_time();
        //std::cout << "step =" << step << "\n";

        //

        // veo star
        //if (step == 0) {
        //    _veo_star_center_explicit();
        //} else {
            _veo_star_center();
        //}
        if (step == 0){
        Plotly plt;
        	//plt.add(SCALARCENTER(*(ns.get_CS("us"))));
        	plt.add(VECTORCENTER(*(this->_veosc)));
        	//plt.add(SCALARCENTER(*(ns.get_CS("u"))));
        	//plt.add(SCALARCENTER(*(ns.get_CS("v"))));
        	plt.add(WireFrame2(*(this->get_grid())));
        	//	auto actor1 = ScalarCenter(*(ns.get_CS("phi")));
        	plt.size(800, 800);
        	//plt.plot();
        }
        // interpolate veo star to face
        Operation::InterpolateC2F(this->_veosc, this->_veosf);

        // div veo face star to center
        spScalar spveos = Operation::Div(this->_veosf);

        // uniform rho
        ASSERT(this->has_event("uniform_rho"));
        Vt rho = this->_values["uniform_rho"];
        //Operation::Multiply(spveos, rho / this->_time->dt());

        // solve poisson equation
        this->_projection->set_depend(this->_css["phi"], nullptr, spveos);

        this->_projection->run();
        this->apply_bc("phi");

        // VectorFieldFace gradP = Grad(PHI);
        spVectorFace gradphi = Operation::Grad(this->_css["phi"]);

        Operation::Copy(this->_css["phi"], this->_css["p"]);
        Operation::Multiply(this->_css["p"], 1.0 / this->_time->dt());
        spScalar lapphi = nullptr;

        if (this->has_event("uniform_mu")){
        	Vt mu = this->_values["uniform_mu"];
        	lapphi = Operation::NablaMuNabla(mu, *(this->_css["phi"]));
        }else{
        	lapphi = Operation::NablaMuNabla(*(this->_css["mu"]), *(this->_css["phi"]));
        }
        Operation::Multiply(lapphi, 0.5/rho);
        Operation::Minus(this->_css["p"], lapphi);

        // correct veo on face
        Operation::Minus(this->_veosf, gradphi);

        // correct veo on center
        spVectorCenter gradphic = Operation::InterpolateF2C(gradphi);
        Operation::Minus(this->_veosc, gradphic);
        // copy velocity
        Operation::Copy(this->_veosc, this->_veoc);
        Operation::Copy(this->_veosf, this->_veof);

        return -1;
    }

    // debug function
    spVectorFace _veo_star_f() {
        return this->_veosf;
    }
    spVectorCenter _veo_star_c() {
        return this->_veosc;
    }
    spVectorFace _veo_f() {
        return this->_veof;
    }
    spVectorCenter _veo_c() {
        return this->_veoc;
    }

protected:
    void _veo_star_center() {
        FOR_EACH_DIM {
            this->apply_bc(this->_nv[d]);
        }

        spScalar spmu = nullptr, sprho = nullptr;
        if (!this->has_event("uniform_mu")) {
            spmu = this->_css["mu"];
        }
        if (!this->has_event("uniform_rho")) {
            sprho = this->_css["rho"];
        }

        Vt dt = this->_time->dt();
        Grid& grid = *(this->_grid);
        VectorFace& vf  = *(this->_veof);
        VectorFace& vmf = *(this->_veomf);
        VectorCenter& sscr = *(this->_sscr);

        FOR_EACH_DIM {
            Scalar& s  = *(this->_css[this->_nv[d]]);
            Scalar& ss = *(this->_css[this->_nvs[d]]);
            Scalar& sm = *(this->_css[this->_nvm[d]]);
            Vt rho   = SMALL;
            Vt smu   = 1.0;
            for (typename Grid::Ijk ijk = grid.begin_ijk(); !ijk.is_end();
                 ++ijk) {
                Index& index = ijk.current();
                //index.show();

                Vt advn   = __advection_term(vf,  s,  index);
                Vt advnm1 = __advection_term(vmf, sm, index);
                Vt adv    = 1.5 * advn - 0.5 * advnm1;

                Vt dif    = __diffusion_term(spmu, s, index);

                if (sprho == nullptr) {
                    rho = this->_values["uniform_rho"];
                } else {
                    rho = (*sprho)(index);
                }
                Vt src  = __source_term(index);
                Vt veos = 0.0;
                if ( this->_time->current_step() == 0){
                	veos = (advn - (0.5 * dif + src) / rho) * dt - s(index);
                }else{
                	veos = (adv  - (0.5 * dif + src) / rho) * dt - s(index);
                }
                sscr[d](index) = veos;
                // -------------------------------
                if (spmu != nullptr) {
                     SHOULD_NOT_REACH;
                }else{
                	smu = this->_values["uniform_mu"];
                	smu = dt * 0.5 / rho;
                }
            }
            // solve
            this->_predict->set_depend(this->_css[this->_nvs[d]], spmu, this->_sscr->get_spdata(d), smu);
            this->_predict->set_uniform_alpha(1);
            this->_predict->run();
        }
        FOR_EACH_DIM {
            this->apply_bc(this->_nvs[d]);
        }
    }


    Vt __advection_term(const VectorFace& vf, const Scalar& s, const Index& index) {
        return this->fun_v_dot_nabla(vf, s, index);
    }

    Vt __diffusion_term(spScalar spmu, const Scalar& s, const Index& index) {
        if (spmu == nullptr) {
            Vt mu = this->_values["uniform_mu"];
            return Operation::NablaMuNabla(mu, s, index);
        } else {
            return Operation::NablaMuNabla(*spmu, s, index);
        }
    }

    Vt __source_term(const Index& index) {
        return 0;
    }

}
;

}

#endif
