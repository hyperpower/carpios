/* 
 *   Matrix Market I/O library for ANSI C
 *
 *   See http://math.nist.gov/MatrixMarket for details.
 *
 *
 */

#ifndef MM_IO_H
#define MM_IO_H

#define MM_MAX_LINE_LENGTH 1025
#define MatrixMarketBanner "%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH 64

typedef char MM_typecode[4];

char *mm_typecode_to_str(MM_typecode matcode);

int mm_read_banner(FILE *f, MM_typecode *matcode);
int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz);
int mm_read_mtx_array_size(FILE *f, int *M, int *N);

int mm_write_banner(FILE *f, MM_typecode matcode);
int mm_write_mtx_crd_size(FILE *f, int M, int N, int nz);
int mm_write_mtx_array_size(FILE *f, int M, int N);

/********************* MM_typecode query fucntions ***************************/

#define mm_is_matrix(typecode)	((typecode)[0]=='M')

#define mm_is_sparse(typecode)	((typecode)[1]=='C')
#define mm_is_coordinate(typecode)((typecode)[1]=='C')
#define mm_is_dense(typecode)	((typecode)[1]=='A')
#define mm_is_array(typecode)	((typecode)[1]=='A')

#define mm_is_complex(typecode)	((typecode)[2]=='C')
#define mm_is_real(typecode)		((typecode)[2]=='R')
#define mm_is_pattern(typecode)	((typecode)[2]=='P')
#define mm_is_integer(typecode) ((typecode)[2]=='I')

#define mm_is_symmetric(typecode)((typecode)[3]=='S')
#define mm_is_general(typecode)	((typecode)[3]=='G')
#define mm_is_skew(typecode)	((typecode)[3]=='K')
#define mm_is_hermitian(typecode)((typecode)[3]=='H')

int mm_is_valid(MM_typecode matcode); /* too complex for a macro */

/********************* MM_typecode modify fucntions ***************************/

#define mm_set_matrix(typecode)	((*typecode)[0]='M')       //M
#define mm_set_coordinate(typecode)	((*typecode)[1]='C')   //C
#define mm_set_array(typecode)	((*typecode)[1]='A')       //A
#define mm_set_dense(typecode)	mm_set_array(typecode)     //A
#define mm_set_sparse(typecode)	mm_set_coordinate(typecode)//C

#define mm_set_complex(typecode)((*typecode)[2]='C')
#define mm_set_real(typecode)	((*typecode)[2]='R')
#define mm_set_pattern(typecode)((*typecode)[2]='P')
#define mm_set_integer(typecode)((*typecode)[2]='I')

#define mm_set_symmetric(typecode)((*typecode)[3]='S')
#define mm_set_general(typecode)((*typecode)[3]='G')
#define mm_set_skew(typecode)	((*typecode)[3]='K')
#define mm_set_hermitian(typecode)((*typecode)[3]='H')

#define mm_clear_typecode(typecode) ((*typecode)[0]=(*typecode)[1]= \
									(*typecode)[2]=' ',(*typecode)[3]='G')

#define mm_initialize_typecode(typecode) mm_clear_typecode(typecode)

/********************* Matrix Market error codes ***************************/

#define MM_COULD_NOT_READ_FILE	11
#define MM_PREMATURE_EOF		12
#define MM_NOT_MTX				13
#define MM_NO_HEADER			14
#define MM_UNSUPPORTED_TYPE		15
#define MM_LINE_TOO_LONG		16
#define MM_COULD_NOT_WRITE_FILE	17

/******************** Matrix Market internal definitions ********************

 MM_matrix_typecode: 4-character sequence

 ojbect 		sparse/   	data        storage
 dense     	type        scheme

 string position:	 [0]        [1]			[2]         [3]

 Matrix typecode:  M(atrix)    C(oord)		R(eal)   	G(eneral)
                   A(array)	   C(omplex)    H(ermitian)
                   P(attern)   S(ymmetric)
                   I(nteger)   K(kew)

 ***********************************************************************/

#define MM_MTX_STR		  "matrix"
#define MM_ARRAY_STR	  "array"
#define MM_DENSE_STR	  "array"
#define MM_COORDINATE_STR "coordinate" 
#define MM_SPARSE_STR	  "coordinate"
#define MM_COMPLEX_STR	  "complex"
#define MM_REAL_STR		  "real"
#define MM_INT_STR		  "integer"
#define MM_GENERAL_STR    "general"
#define MM_SYMM_STR		  "symmetric"
#define MM_HERM_STR		  "hermitian"
#define MM_SKEW_STR		  "skew-symmetric"
#define MM_PATTERN_STR    "pattern"

/*  high level routines */

int mm_write_mtx_crd(char fname[], int M, int N, int nz, int I[], int J[],
		double val[], MM_typecode matcode);
int mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, int I[], int J[],
		double val[], MM_typecode matcode);
int mm_read_mtx_crd_entry(FILE *f, int *I, int *J, double *real, double *img,
		MM_typecode matcode);

int mm_read_unsymmetric_sparse(const char *fname, int *M_, int *N_, int *nz_,
		double **val_, int **I_, int **J_);

// High level function ======================
#include "../algebra/matrix.hpp"
#include "io_define.hpp"
#include "../algebra/matrix_SparCoord.hpp"
#include <string>

namespace carpio {

template<class VALUE>
void mm_read_mtx_dense(std::string filename, MatrixV<VALUE>& m) {
	FILE* f;
	if (file_access_check(filename, 4)) { //enable read
		f = fopen(filename.c_str(), "r");
	} else {
		std::cerr << " !> Can't read file " << filename << "\n";
		exit(-1);
	}

	MM_typecode matcode;
	int M;
	int N;

	if (mm_read_banner(f, &matcode) != 0) {
		std::cerr << " !> Could not process Matrix Market banner.\n";
		exit(-1);
	}

	if (!mm_is_matrix(matcode) && !mm_is_dense(matcode) && !mm_is_real(matcode)) {
		std::cerr << " !> This function is for dense matrix\n";
		std::cerr << " !> but, the input file is "
				<< mm_typecode_to_str(matcode) << "\n";
		exit(-1);
	}
	// read dense matrix size
	mm_read_mtx_array_size(f, &M, &N);
	// check size
	if (M <= 0 || N <= 0) {
		std::cerr << " !> Matrix size error: iLen=" << M << "jLen=" << N
				<< "\n";
		exit(-1);
	}

	//reconstruct
	m.reconstruct(M, N);

	//scan
	for (int j = 0; j < N; j++) {
		for (int i = 0; i < M; i++) {
			fscanf(f, "%lf\n", &m[i][j]);
		}
	}
}

template<class VALUE>
void mm_read_array(std::string filename, ArrayListV<VALUE>& arr) {
	FILE* f;
	if (file_access_check(filename, 4)) { //enable read
		f = fopen(filename.c_str(), "r");
	} else {
		std::cerr << " !> Can't read file " << filename << "\n";
		exit(-1);
	}

	MM_typecode matcode;
	int M;
	int N;

	if (mm_read_banner(f, &matcode) != 0) {
		std::cerr << " !> Could not process Matrix Market banner.\n";
		exit(-1);
	}

	if (!mm_is_matrix(matcode) && !mm_is_dense(matcode) && !mm_is_real(matcode)) {
		std::cerr << " !> This function is for Array\n";
		std::cerr << " !> but, the input file is "
				<< mm_typecode_to_str(matcode) << "\n";
		exit(-1);
	}
	// read dense matrix size
	mm_read_mtx_array_size(f, &M, &N);
	// check size
	if (M <= 0 || N != 1) {
		std::cerr << " !> Array size error: iLen=" << M << "jLen=" << N << "\n";
		exit(-1);
	}

	//reconstruct
	arr.reconstruct(M);

	//scan
	for (int i = 0; i < M; i++) {
		fscanf(f, "%lf\n", &arr[i]);
	}
}

template<class VALUE>
void mm_read_mtx_sparse(std::string filename, MatrixSCO_<VALUE>& m) {
	FILE* f;
	if (file_access_check(filename, 4)) { //enable read
		f = fopen(filename.c_str(), "r");
	} else {
		std::cerr << " !> Can't read file " << filename << "\n";
		exit(-1);
	}

	MM_typecode matcode;
	int M;
	int N;
	int nz;

	if (mm_read_banner(f, &matcode) != 0) {
		std::cerr << " !> Could not process Matrix Market banner.\n";
		exit(-1);
	}

	if (!mm_is_matrix(matcode) && !mm_is_sparse(matcode) && !mm_is_real(matcode)) {
		std::cerr << " !> This function is for MatrixSCO \n";
		std::cerr << " !> but, the input file is "
				<< mm_typecode_to_str(matcode) << "\n";
		exit(-1);
	}
	// read dense matrix size
	mm_read_mtx_crd_size(f, &M, &N, &nz);
	// check size
	if (M <= 0 || N <= 0) {
		std::cerr << " !> Matrix size error: iLen=" << M << "jLen=" << N
				<< "\n";
		exit(-1);
	}

	//scan
	if (mm_is_symmetric(matcode)) {
		//reconstruct
		m.newsize(M, N, nz);
		int count_ndia = 0;
		for (int i = 0; i < nz; i++) {
			int rind, cind;
			double value;
			fscanf(f, "%d %d %le\n", &rind, &cind, &value);
			m.row_ind(i) = rind - 1;
			m.col_ind(i) = cind - 1;
			m.val(i) = typename MatrixSCO_<VALUE>::vt(value);
			if (rind != cind) {
				count_ndia++;
			}
		}
		m.resize(M, N, nz + count_ndia);
		int ii = 0;
		for (int i = 0; i < nz; i++) {
			if (m.row_ind(i) != m.col_ind(i)) {
				m.row_ind(nz + ii) = m.col_ind(i);
				m.col_ind(nz + ii) = m.row_ind(i);
				m.val(nz + ii) = m.val(i);
				++ii;
			}
		}
	} else {
		m.newsize(M, N, nz);
		for (int i = 0; i < nz; i++) {
			int rind, cind;
			double val = 0;
			fscanf(f, "%d %d %le\n", &rind, &cind, &val);
			m.val(i) = typename MatrixSCO_<VALUE>::vt(val);
			m.row_ind(i) = rind - 1;
			m.col_ind(i) = cind - 1;
		}
	}
}

template<class VALUE>
void mm_write_mtx_sparse(std::string filename, MatrixSCR_<VALUE>& m) {

	FILE* f = fopen(filename.c_str(), "w");

	MatrixSCO_<VALUE> mco(m);
	MM_typecode matcode;

	mm_initialize_typecode(&matcode);
	mm_set_sparse(&matcode);
	mm_set_coordinate(&matcode);
	mm_set_real(&matcode);

	mm_write_banner(f, matcode);

	mm_write_mtx_crd_size(f, m.size1(), m.size2(), m.NumNonzeros());

	/* NOTE: matrix market files use 1-based indices, i.e. first element
	 of a vector has index 1, not 0.  */

	for (St i = 0; i < m.NumNonzeros(); i++)
		fprintf(f, "%d %d %20.10e\n", int(mco.row_ind(i) + 1), int(mco.col_ind(i) + 1), mco.val(i));

}

template<class VALUE>
void mm_write_array(std::string filename, ArrayListV<VALUE>& m) {

	FILE* f = fopen(filename.c_str(), "w");

	MM_typecode matcode;

	mm_initialize_typecode(&matcode);
	mm_set_matrix(&matcode);
	mm_set_array(&matcode);
	mm_set_real(&matcode);
	mm_set_general(&matcode);

	mm_write_banner(f, matcode);

	fprintf(f, "%d %d\n", m.size(), 1);
	//mm_write_mtx_crd_size(f, m.size(), 1, m.size());

	/* NOTE: matrix market files use 1-based indices, i.e. first element
	 of a vector has index 1, not 0.  */

	for (St i = 0; i < m.size(); i++)
		fprintf(f, "%20.10e\n", m[i]);

}

}
// ==========================================

#endif
