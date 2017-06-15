/*
 * _test_solver.hpp
 *
 *  Created on: Jun 2, 2017
 *      Author: zhou
 */

#ifndef _TEST_SOLVER_HPP_
#define _TEST_SOLVER_HPP_

#include "gtest/gtest.h"
#include "../../lib/algebra/matrix_small.hpp"
#include "../../lib/algebra/algebra.hpp"
#include "../../lib/io/mmio.h"

#include <math.h>
#include <string>
#include <memory>

namespace carpio {

void read_matrix(std::string filename) {
	int ret_code;
	MM_typecode matcode;
	FILE *f;
	int M, N, nz;
	int i, *I, *J;
	double *val;

	//string filename = "mat1.txt";
	f = fopen(filename.c_str(), "r");
	if (f == nullptr) {
		std::cout << "matrix file not found!\n";
		exit(1);
	}
	if (mm_read_banner(f, &matcode) != 0) {
		printf("Could not process Matrix Market banner.\n");
		exit(1);
	}
	std::cout << mm_typecode_to_str(matcode) << std::endl;
	/*  This is how one can screen matrix types if their application */
	/*  only supports a subset of the Matrix Market data types.      */
	//cout << mm_is_sparse(matcode) << endl;
	//cout << mm_is_dense(matcode) << endl;
	if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)) {
		printf("Sorry, this application does not support ");
		printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
		exit(1);
	}

	/* find out size of sparse matrix .... */
	if ((ret_code = mm_read_mtx_array_size(f, &M, &N)) != 0)
		exit(1);

	/* reseve memory for matrices */

	I = (int *) malloc(nz * sizeof(int));
	J = (int *) malloc(nz * sizeof(int));
	val = (double *) malloc(nz * sizeof(double));

	/* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
	/*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
	/*   (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

	for (i = 0; i < M * N; i++) {
		fscanf(f, "%lf\n", &val[i]);
		I[i]--; /* adjust from 1-based to 0-based */
		J[i]--;
	}

	if (f != stdin)
		fclose(f);

	/************************/
	/* now write out matrix */
	/************************/

	mm_write_banner(stdout, matcode);
	mm_write_mtx_array_size(stdout, M, N);
	for (i = 0; i < M * N; i++)
		fprintf(stdout, "%10.5f\n", val[i]);

}

TEST(Matrix_s, read) {
	std:;string workdir = "./test/input_files/";
	MatrixSCO_<Float> mf;
	std:;string fn_matrix = "685_bus";
	mm_read_mtx_sparse(workdir + fn_matrix + ".mtx", mf);
	std::cout << mf.NumNonzeros() << endl;
	std::cout << "matrix information  oo ==============\n";
	std::cout << "i = " << mf.iLen() << "   j = " << mf.jLen() << endl;
	//mf.show(0);
	MatrixSCR_<Float> mfr(mf);
	std::cout << "matrix information  cr ==============\n";
	std::cout << "i = " << mfr.iLen() << "   j = " << mfr.jLen() << endl;

	std::cout << "matrix information  ==============\n";
	int i = 0, j = 0;
	std::cout << "value at (" << i << ", " << j << ")";
	std::cout << "  oo " << mf(i, j) << "  cr " << mfr(i, j) << endl;
	ASSERT_EQ(mf(i, j), mfr(i, j));

	//arrayListV<Float> ax = mf*x;
	//ax.show();
	//cout.precision(4);
	//cout << mf(20, 22) << "   " << mfr(20, 22) << " \n";

	//set up ========
	//int max_iter = 1000;
	//Float tol = 1e-6;
	//ListT_<Float> lr;  //list residual
	//solver =======================
	//IC_BiCGSTAB(mfr, x, b, max_iter, tol, lr);
	//cout << "max iter = " << max_iter << endl;
	//cout << "tol      = " << tol << endl;
	//gnuplot_show_ylog(lr);
	//cout << "solver jacobi " << endl;
	//x.assign(1);
	//lr.clear();
	//max_iter = 100000;
	//tol = 1e-6;
	//Dia_BiCGSTAB(mfr, x, b, max_iter, tol, lr);
	//cout << "max iter = " << max_iter << endl;
	//cout << "tol      = " << tol << endl;
	//gnuplot_show_ylog(lr);
}

TEST(Matrix_s, solve) {
	string workdir = "./test/input_files/";
	MatrixSCO_<Float> mf;
	string fn_matrix = "685_bus";
	// read matrix ----------------------------------------
	mm_read_mtx_sparse(workdir + fn_matrix + ".mtx", mf);
	std::cout << "None zeros : " << mf.NumNonzeros() << endl;
	MatrixSCR_<Float> mfr(mf);

	// assign x and b -------------------------------------
	ArrayListV<Float> b(mfr.iLen());
	b.assign(1);
	ArrayListV<Float> x(mfr.iLen());
	x.assign(1);

	// solve ----------------------------------------------
	int max_iter = 100000;
	Float tol = 1e-6;
	std::list<Float> lr;  //list residual
	//
	//typedef Solver_SOR_<Float> Solver;
	typedef Solver_Jacobi_<Float> Solver;
	Solver solver(max_iter, tol);

	solver.solve(mfr, x, b);

	//Jacobi(mfr, x, b, max_iter, tol, lr);
	//
	std::cout << "max iter = " << solver.num_iter() << endl;
	std::cout << "tol      = " << solver.residual() << endl;
	//gnuplot_show_ylog(lr);
	// out put --------------------------------------------
	mm_write_array("x.txt", x);

}

TEST(matrix, write_sparse) {
	MM_typecode matcode;
	const int nz = 4;
	const int M = 10;
	const int N = 10;
	int I[nz] = { 0, 4, 2, 8 };
	int J[nz] = { 3, 8, 7, 5 };
	double val[nz] = { 1.1, 2.2, 3.2, 4.4 };
	int i;

	mm_initialize_typecode(&matcode);
	mm_set_matrix(&matcode);
	mm_set_sparse(&matcode);
	mm_set_real(&matcode);

	mm_write_banner(stdout, matcode);
	mm_write_mtx_crd_size(stdout, M, N, nz);

	/* NOTE: matrix market files use 1-based indices, i.e. first element
	 of a vector has index 1, not 0.  */

	for (i = 0; i < nz; i++)
		fprintf(stdout, "%d %d %10.3g\n", I[i] + 1, J[i] + 1, val[i]);
}

}

#endif /* _TEST_SOLVER_HPP_ */
