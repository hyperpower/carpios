#ifndef SCALAR_H_
#define SCALAR_H_

#include "../carpio_define.hpp"
#include "domain/domain.hpp"

#include <functional>
#include <math.h>

namespace carpio {
/*
 * New center data on leaf
 * len : the length of array
 */
//void NewCenterDataOnLeaf(Grid_2D& grid, st len);
//void NewCenterDataOnLeaf(Grid_3D& grid, st len);
/*
 * The input indicates the coordinate location
 * return the value on the location
 */
typedef Vt (*scalar_pfun)(Cvt x, Cvt y, Cvt z);

/*
 * set value on the leaf of grid
 *
 * This function will change the value on the grid
 */
void SetScalarOnCenterLeaf( // 2D tree
		Grid_2D& grid,                    //pQuadTree
		St idx,                           //data index
		scalar_pfun pf  //data plus
		);
}

#endif
