#include "../../lib/algebra/matrix_small.hpp"
#include "../../lib/algebra/algebra.hpp"
#include "../../lib/io/mmio.h"

#include <math.h>
#include <string>
#include <sstream>
#include <vector>
#include <memory>


using namespace carpio;

void solve_one(const std::string& name, 
               MatrixSCR_<Float>& A,
               double o){

    typedef Solver_<Float> Solver;
    int max_iter = 100000;
    double tol = 1e-6;
    Solver* ps;
    typedef Solver_Jacobi_<Float> Jac;
    Jac jac(max_iter, tol);
    typedef Solver_SOR_<Float> Sor;
    Sor sor(max_iter, tol, o);
    if(name == "Jacobi"){
        ps = &jac;
    }
    if(name == "SOR"){
        ps = &sor;
    }
    

    // assign x and b -------------------------------------
    ArrayListV<Float> b(A.iLen());
    b.assign(1);
    ArrayListV<Float> x(A.iLen());
    x.assign(1);

    ps->solve(A, x, b);

    ArrayListV<Float> residual = ps->get_residual_array();

    std::stringstream sst;
    if(name == "SOR"){
        sst << name << "_" << std::setprecision(2) << std::fixed <<o;
    } else {
        sst << name;
    }
    std::string fh = sst.str();

    std::cout << "Solver name : " << fh << "\n";

    std::cout << "  max iter  = " << ps->max_iter();
    std::cout << "  tolerance = " << ps->tolerance() << std::endl;
    std::cout << "  num iter  = " << ps->num_iter();
    std::cout << "  residual  = " << ps->residual() << std::endl;

    mm_write_array(fh + "_x.txt", x);
    mm_write_array(fh + "_residual.txt", residual);

}

int main() {
    std::string workdir = "./";
    MatrixSCO_<Float> mf;
    std::string fn_matrix = "685_bus";
    // read matrix ----------------------------------------
    mm_read_mtx_sparse(workdir + fn_matrix + ".mtx", mf);
    std::cout << "None zeros : " << mf.NumNonzeros() << std::endl;
    MatrixSCR_<Float> mfr(mf);

    solve_one("Jacobi", mfr, 1.0);

    std::vector<double> vomega = {1.0, 1.2, 1.5, 1.7, 1.9};
    for(auto& val : vomega){
        solve_one("SOR", mfr, val);
    }
}
