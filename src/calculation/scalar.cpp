#include "scalar.h"

namespace carpio {
void _TraversalLeafIterator(Grid_2D& grid, std::function<void(pNode_2D)>& fun) {
	for (Grid_2D::iterator_leaf iter = grid.begin_leaf();
			iter != grid.end_leaf(); ++iter) {
		Grid_2D::pNode pn = iter.get_pointer();
		fun(pn);
	}
}
void _TraversalLeafIterator(Grid_3D& grid, std::function<void(pNode_3D)>& fun) {
	for (Grid_3D::iterator_leaf iter = grid.begin_leaf();
			iter != grid.end_leaf(); ++iter) {
		Grid_3D::pNode pn = iter.get_pointer();
		fun(pn);
	}
}
void SetScalarOnCenterLeaf( // 2D tree
		Grid_2D& grid,                    //pQuadTree
		St idx,                           //data index
		scalar_pfun pf //data plus
		) {
	std::function<void(pNode_2D)> fun = [&pf, &idx](pNode_2D pn) {
		// for 2D z is unsless
			pn->cd(idx) = pf(pn->cp(_X_), pn->cp(_Y_), 0.0);
		};
	_TraversalLeafIterator(grid, fun);
}

void NewCenterDataOnLeaf(Grid_2D& grid, St len) {
	std::function<void(pNode_2D)> fun = [&len](pNode_2D pn) {
			if(pn->data!=nullptr) {   //resize
				pn->data->resize_center(len);
			} else {                  //new
				pn->data = new Grid_2D::Data(len, 0, 0, nullptr);
			}
		};
	_TraversalLeafIterator(grid, fun);
}

void NewCenterDataOnLeaf(Grid_3D& grid, St len) {
	std::function<void(pNode_3D)> fun = [&len](pNode_3D pn) {
			if(pn->data!=nullptr) {   //resize
				pn->data->resize_center(len);
			} else {                  //new
				pn->data = new Grid_3D::Data(len, 0, 0, nullptr);
			}
		};
	_TraversalLeafIterator(grid, fun);
}

}
