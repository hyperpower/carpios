#ifndef _S_NS_EXPLICIT_HPP
#define _S_NS_EXPLICIT_HPP

#include "s_data.hpp"
#include "s_define.hpp"
#include "s_stencil.hpp"
#include "s_equation.hpp"
#include "s_solver.hpp"
#include "s_vector.hpp"
#include "s_io_plotly.hpp"
#include "s_ns.hpp"

namespace structure {

template<St DIM>
class NS_explicit_: public NS_<DIM>
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
    typedef Operation_<DIM> Operation;

    typedef carpio::Any Any;

protected:
    spVectorCenter _veosc;
    spVectorFace _veosf;

public:
    NS_explicit_(spGrid spg) :
        Base(spg) {
        this->_default_setup();
        Function fun = [](Vt, Vt, Vt, Vt) {return 0;};
        std::string veo_name[] = { "us", "vs", "ws" };
        std::string vf_name[] = { "usf", "vsf", "wsf" };
        FOR_EACH_DIM {
            this->_functions[veo_name[d]] = fun;
            this->_functions[vf_name[d]] = fun;
        }
        //
        this->_functions["phi"] = fun;

    }

    int initial() {
        this->_ns_initial();
        // initial field veo star on center and face
        std::string veo_name[] = { "u", "v", "w" };
        std::string veos_name[] = { "us", "vs", "ws" };
        std::string veosf_name[] = { "usf", "vsf", "wsf" };
        spScalar vc[] = { nullptr, nullptr, nullptr };
        spScalar vf[] = { nullptr, nullptr, nullptr };
        spGrid spg = this->_grid;
        // phi -------------------------------------
        spScalar spphi = spScalar(new Scalar(spg));
        this->_css["phi"] = spphi;
        this->set_CS("phi", this->_functions["phi"]);
        // copy boundary condition ----------------
        FOR_EACH_DIM {
            this->_bi->copy(veo_name[d], veos_name[d]);
        }
        FOR_EACH_DIM {
            vc[d] = spScalar(new Scalar(spg));
            this->_css[veos_name[d]] = vc[d];
            this->set_CS(veos_name[d], this->_functions[veos_name[d]]);
            vf[d] = spScalar(new Scalar(spg));
            this->_css[veosf_name[d]] = vf[d];
            this->set_CS(veosf_name[d], this->_functions[veosf_name[d]]);
        }
        _veosc = spVectorCenter(new VectorCenter(vc[0], vc[1], vc[2]));
        _veosf = spVectorFace(new VectorFace(vf[0], vf[1], vf[2]));
        return 1;
    }

    int run_one_step(St step) {
        // debug
        int stepo = 0;
        // restrict time
        //
        this->restrict_time();
        //std::cout << "step =" << step << "\n";
        //
        // veo star
        _veo_star_center();

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
        // correct veo on face
        //Operation::Multiply(gradphi, this->_time->dt() / rho);
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
        // std::cout << "    veo star center\n";
        //UC.GetData(dim) -= (V_dot_Nabla(UF, UC.GetData(dim), scheme) - Laplacian(UC.GetData(dim))/Re)*timer->GetTimeStep();
        const std::string veo_name[] = { "u", "v", "w" };
        const std::string veos_name[] = { "us", "vs", "ws" };
        FOR_EACH_DIM {
            this->apply_bc(veo_name[d]);
        }
        spVectorFace spvf = this->_veof;
        spScalar spmu = nullptr, sprho = nullptr;
        if (!this->has_event("uniform_mu")) {
            spmu = this->_css["mu"];
        }
        if (!this->has_event("uniform_rho")) {
            sprho = this->_css["rho"];
        }
        Vt dt = this->_time->dt();
        Grid& grid = *(this->_grid);
//#pragma omp parallel for schedule(dynamic,1)
        FOR_EACH_DIM {
            spScalar sps = this->_css[veo_name[d]];
            spScalar spss = this->_css[veos_name[d]];
            for (typename Grid::Ijk ijk = grid.begin_ijk(); !ijk.is_end();
                 ++ijk) {
                Index& index = ijk.current();
                //index.show();
                Vt adv = __advection_term(spvf, sps, index);
                Vt dif = __diffusion_term(spmu, sps, index);
                Vt rho = SMALL;
                if (sprho == nullptr) {
                    rho = this->_values["uniform_rho"];
                } else {
                    rho = (*sprho)(index);
                }
                Vt src = __source_term(index);
                Vt veos = (-adv + (dif + src) / rho) * dt + (*sps)(index);
                (*spss)(index) = veos;
                //
#ifdef __DEBUG__
                if (index == index2) {
                    std::cout << " dim =  " << d << std::endl;
                    std::cout << " name = " << veo_name[d] << std::endl;
                    std::cout << " rho = " << rho;
                    std::cout << " veos = " << veos;
                    std::cout << " adv " << adv << " dif " << dif << std::endl;
                    std::cout << " uo  " << (*sps)(index);
                    std::cout << " un  " << (*spss)(index);
                    std::cout << std::endl;
                }
#endif
            }

        }
        FOR_EACH_DIM {
            this->apply_bc(veos_name[d]);
        }
    }

    Vt __advection_term(spVectorFace spvf, spScalar sps, const Index& index) {
        return this->fun_v_dot_nabla(*spvf, *sps, index);
    }

    Vt __diffusion_term(spScalar spmu, spScalar sps, const Index& index) {
        if (spmu == nullptr) {
            Vt mu = this->_values["uniform_mu"];
            return Operation::NablaMuNabla(mu, *sps, index);
        } else {
            return Operation::NablaMuNabla(*spmu, *sps, index);
        }
    }

    Vt __source_term(const Index& index) {
        return 0;
    }

}
;

}

#endif
