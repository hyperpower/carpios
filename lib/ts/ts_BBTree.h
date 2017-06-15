/************************
 //  \file   BBTree.h
 //  \brief
 // 
 //  \author czhou
 //  \date   17 juin 2015 
 ***********************/
#ifndef TS_BBTREE_H_
#define TS_BBTREE_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"
#include "ts_segment.h"
#include <functional>

namespace TS {

//enum BBNodeType {
//	_SEPERATE = 0, _OVERLAP = 1
//};

template<typename BOX>
class BBNode {
protected:
	typedef BBNode<BOX> _Self;
	typedef BBNode<BOX>* pNode;
	typedef st size_type;
public:
	BOX _box;
//	BBNodeType type;

	pNode lchild;
	pNode rchild;
	pNode father;

	BBNode() {
		_box = BOX();
		lchild = nullptr;
		rchild = nullptr;
		father = nullptr;
		//	type = _SEPERATE;
	}

	BBNode(const BOX &data) {
		_box = BOX(data);
		lchild = nullptr;
		rchild = nullptr;
		father = nullptr;
		//	type = _SEPERATE;
	}
	BBNode(const BOX &data, pNode lc, pNode rc, pNode fc) {
		_box = BOX(data);
		lchild = lc;
		rchild = rc;
		father = fc;
		//	type = _SEPERATE;
	}
	pNode leftmost() {
		pNode ptn = this;
		while (ptn->lchild != nullptr)
			ptn = ptn->lchild;
		return ptn;
	}
	pNode rightmost() {
		pNode ptn = this;
		while (ptn->rchild != nullptr)
			ptn = ptn->rchild;
		return ptn;
	}

	size_type _height(const BBNode<BOX>* cur) const {
		if (cur == nullptr)
			return 0;
		else
			return 1 + MAX(_height(cur->lchild), _height(cur->rchild));
	}
	size_type height() const {
		return 0 + MAX(_height(this->lchild), _height(this->rchild));
	}

	size_type _level(const BBNode<BOX>* cur) const {
		if (cur == nullptr)
			return 0;
		else
			return 1 + _level(cur->father);
	}
	size_type level() const {
		// root is root->lchild
		return 0 + _level(father) - 1;
	}

	bool is_leaf() const {
		_return_val_if_fail(this->lchild == nullptr, false);
		_return_val_if_fail(this->rchild == nullptr, false);
		return true;
	}
	BOX& box(){
		return this->_box;
	}
	const BOX& box() const{
		return this->_box;
	}
};

//This is the end of Class BinaryTreeNode
//-----------------------------------------------

template<class _Tp, class _Ref, class _Ptr>
class _BinaryTree_iterator {
public:
	typedef st size_type;
	typedef st difference_type;
	//typedef bidirectional_iterator_tag iterator_category;

	//typedef _BinaryTree_iterator<_Tp, _Tp&, _Tp*> iterator;
	//typedef _BinaryTree_iterator<_Tp, const _Tp&, const _Tp*> const_iterator;
	typedef _BinaryTree_iterator<_Tp, _Ref, _Ptr> _Self;

	typedef _Tp value_type;
	typedef _Ptr pointer;
	typedef _Ref reference;

	_Tp* _ptr;

	_BinaryTree_iterator() {
		_ptr = nullptr;
	}
	_BinaryTree_iterator(_Tp* _x) {
		this->_ptr = _x;
	}
	_BinaryTree_iterator(const _Self& _x) {
		this->_ptr = _x._ptr;
	}

	void _incr() {
		_Tp* fn;
		if (_ptr != nullptr) {
			if (_ptr->rchild != nullptr) {
				_ptr = _ptr->rchild->leftmost();
				return;
			}
			fn = _ptr->father;
			while (fn && _ptr == fn->rchild) {
				_ptr = fn;
				fn = fn->father;
			}
			_ptr = fn;
		}
	}

	void _decr() {
		_Tp* fn;
		if (_ptr != NULL) {
			if (_ptr->lchild) {
				_ptr = _ptr->lchild->rightmost();
				return;
			}
			fn = _ptr->father;
			while (fn && _ptr == fn->lchild) {
				_ptr = fn;
				fn = _ptr->father;
			}
		}
	}

	bool operator==(const _Self& _x) const {
		return _ptr == _x._ptr;
	}
	bool operator!=(const _Self& _x) const {
		return _ptr != _x._ptr;
	}

	reference operator*() const {
		return (*_ptr);
	}

	pointer operator->() const {
		return &(operator*());
	}

	_Self& operator++() {
		this->_incr();
		return *this;
	}

	_Self operator++(int) {
		_Self __tmp = *this;
		this->_incr();
		return __tmp;
	}

	_Self& operator--() {
		this->_decr();
		return *this;
	}

	_Self operator--(int) {
		_Self __tmp = *this;
		this->_decr();
		return __tmp;
	}

	bool isExist() {
		return _ptr != nullptr;
	}
};

//===============================================
template<typename BOX>
class BBTree {
public:
	typedef BBNode<BOX> Node;
	typedef BBNode<BOX>* pNode;
	typedef const BBNode<BOX>* const_pNode;
	typedef Node& reference;
	typedef const Node& const_reference;
	typedef Node* pointer;
	typedef const Node* const_pointer;
	typedef _BinaryTree_iterator<Node, Node&, Node*> iterator;
	typedef _BinaryTree_iterator<Node, const Node&, const Node*> const_iterator;

	typedef st difference_type;
	typedef st size_type;

	typedef std::function<void(pNode)> Fun;
	typedef std::function<void(pNode, bool&)> Fun_flag;
	typedef std::function<void(const_pNode, bool&)> Func_flag;
	;
protected:
	typedef void (*pFun_BBTree)(pNode, utPointer);
	typedef void (*pFun_BBTree_flag)(Int&, pNode, utPointer);

	pNode _root;

	void _preorder(pNode, pFun_BBTree, utPointer);
	void _preorder_flag(Int&, pNode, pFun_BBTree_flag, utPointer);

	void _preorder_flag(bool&, pNode, Fun_flag);
	void _preorder_flag(bool&, const_pNode, Func_flag) const;

	void _inorder(pNode, pFun_BBTree, utPointer);
	void _inorder(pNode, Fun);
	void _inorder(pNode, Fun) const;
	void _postorder(pNode, pFun_BBTree, utPointer);

	void _destory(pNode& cur) {
		if (cur != nullptr) {
			_destory(cur->lchild);
			_destory(cur->rchild);
			delete cur;
			cur = nullptr;
		}
	}
	void _copy(pNode& cur, const pNode& original) {
		if (cur == nullptr) {
			cur = new Node(original->m_value);
		}
		if (original->lchild != nullptr) {
			cur->lchild = new Node(original->lchild->m_value);
			cur->lchild->father = cur;
			_copy(cur->lchild, original->lchild);
		}
		if (original->rchild != nullptr) {
			cur->rchild = new Node(original->rchild->m_value);
			cur->rchild->father = cur;
			_copy(cur->rchild, original->rchild);
		}
	}
	void _reverse(pNode cur) {
		if (cur != nullptr) {
			pNode temp = cur->lchild;
			cur->lchild = cur->rchild;
			cur->rchild = temp;
			_reverse(cur->lchild);
			_reverse(cur->rchild);
		}
	}
	size_type _height(pNode cur) const {
		if (cur == nullptr)
			return 0;
		else
			return 1 + MAX(_height(cur->lchild), _height(cur->rchild));
	}
	size_type _size(const pNode& cur) const {
		if (cur == nullptr) {
			return 0;
		} else {
			return _size(cur->lchild) + _size(cur->rchild) + 1;
		}
	}
public:
	iterator root() {
		return _root->lchild;
	}
	const_iterator root() const {
		return _root->lchild;
	}
	//constructor================================
	BBTree();
	BBTree(const BBTree<BOX>&);

	BBTree(const Set<BOX>&);
	//
	void bottomup_build(List<pNode>& lc, List<pNode>& lf, Int& flag);
	//destructor ================================
	~BBTree();
	//operator= =================================
	BBTree<BOX>& operator=(const BBTree<BOX>&);
	//Traversal =================================
	void PreOrder_flag(Int&, pFun_BBTree_flag, utPointer);
	void PreOrder_flag(Fun_flag);
	void PreOrder_flag(Func_flag) const;
//protected:
//	template<class Args>
//	void _preorder(pNode cur, std::function<void(pNode, Args)> fun,
//			Args &args) {
//		if (cur != nullptr) {
//			fun(cur, args);
//			_preorder(cur->lchild, fun, args);
//			_preorder(cur->rchild, fun, args);
//		}
//	}
public:
	void PreOrder(pFun_BBTree, utPointer);
	void InOrder(pFun_BBTree, utPointer);
	void InOrder(Fun);
	void InOrder(Fun) const;
	void PostOrder(pFun_BBTree, utPointer);
	/*
	 *  iterator===================================
	 */
	iterator begin() {
		pNode pnode = _root->lchild;
		if (pnode == nullptr) {
			return _root;
		}
		if (pnode->lchild) {
			while (pnode->lchild)
				pnode = pnode->lchild;
		}
		return pnode;
	}
	const_iterator begin() const {
		pNode pnode = _root->lchild;
		if (pnode == nullptr) {
			return _root;
		}
		if (pnode->lchild) {
			while (pnode->lchild)
				pnode = pnode->lchild;
		}
		return pnode;
	}
	iterator end() {
		return _root;
	}
	const_iterator end() const {
		return _root;
	}
	//===========================================

	bool empty() const {
		return _root->lchild == nullptr;
	}
	size_type size() const {
		return 0 + _size(_root->lchild);
	}
	void reverse() {
		_reverse(_root->lchild);
	}
	void clear() {
		_destory(_root->lchild);
	}
	size_type height() const {
		return 0 + _height(_root->lchild);
	}

	//show ======================================
protected:
	void _printtree(BBNode<BOX> *cur, int layer);
public:
	void print_tree();
	void output_vtk_height(const String& fn, Int hei);

};

template<typename BOX>
BBTree<BOX>::BBTree() {
	_root = nullptr;
}

template<typename BOX>
BBTree<BOX>::BBTree(const BBTree<BOX>& a) {
	this->_root = nullptr;
	_copy(this->_root, a._root);
}

template<typename BOX>
BBTree<BOX>::~BBTree() {
	_destory(_root);
}

template<typename BOX>
BBTree<BOX>& BBTree<BOX>::operator=(const BBTree<BOX>& original) {
	_destory(this->_root);
	this->_root = nullptr;
	_copy(this->_root, original._root);
	return *this;
}

template<typename BOX>
void BBTree<BOX>::_preorder(pNode cur, pFun_BBTree visit, utPointer utp) {
	if (cur != nullptr) {
		(*visit)(cur->m_value, utp);
		_preorder(cur->lchild, visit, utp);
		_preorder(cur->rchild, visit, utp);
	}
}

template<class BOX>
void BBTree<BOX>::PreOrder(pFun_BBTree visit, utPointer utp) {
	_preorder(_root, visit, utp);
}

template<typename BOX>
void BBTree<BOX>::_preorder_flag(Int& flag, pNode cur, pFun_BBTree_flag visit,
		utPointer utp) {
	if (cur != nullptr) {
		(*visit)(flag, cur, utp);
		if (flag == -1) { //-------------------------
			return;
		}
		_preorder_flag(flag, cur->lchild, visit, utp);
		_preorder_flag(flag, cur->rchild, visit, utp);
	}
}
template<typename BOX>
void BBTree<BOX>::_preorder_flag(bool& flag, pNode cur, Fun_flag visit) {
	if (cur != nullptr) {
		visit(cur, flag);
		if (!flag) { //-------------------------
			return;
		}
		_preorder_flag(flag, cur->lchild, visit);
		_preorder_flag(flag, cur->rchild, visit);
	}
}

template<typename BOX>
void BBTree<BOX>::_preorder_flag(bool& flag, const_pNode cur,
		Func_flag visit) const {
	if (cur != nullptr) {
		visit(cur, flag);
		if (!flag) { //-------------------------
			return;
		}
		_preorder_flag(flag, cur->lchild, visit);
		_preorder_flag(flag, cur->rchild, visit);
	}
}

template<class BOX>
void BBTree<BOX>::PreOrder_flag(Int& flag, pFun_BBTree_flag visit,
		utPointer utp) {
	_preorder_flag(flag, _root->lchild , visit, utp);
}
template<class BOX>
void BBTree<BOX>::PreOrder_flag(Fun_flag visit) {
	bool flag = true;
	_preorder_flag(flag, _root->lchild , visit);
}

template<class BOX>
void BBTree<BOX>::PreOrder_flag(Func_flag visit) const {
	bool flag = true;
	_preorder_flag(flag, _root->lchild , visit);
}

template<typename BOX>
void BBTree<BOX>::_postorder(pNode cur, pFun_BBTree visit, utPointer utp) {
	if (cur != nullptr) {
		_postorder(cur->lchild, visit, utp);
		_postorder(cur->rchild, visit, utp);
		(*visit)(cur->m_value, utp);
	}
}

template<class BOX>
void BBTree<BOX>::PostOrder(pFun_BBTree visit, utPointer utp) {
	_postorder(_root, visit, utp);
}

template<typename BOX>
void BBTree<BOX>::_inorder(pNode cur, pFun_BBTree visit, utPointer utp) {
	if (cur != nullptr) {
		_inorder(cur->lchild, visit, utp);
		(*visit)(cur, utp);
		_inorder(cur->rchild, visit, utp);
	}
}
template<typename BOX>
void BBTree<BOX>::_inorder(pNode cur, Fun visit) {
	if (cur != nullptr) {
		_inorder(cur->lchild, visit);
		visit(cur);
		_inorder(cur->rchild, visit);
	}
}
template<typename BOX>
void BBTree<BOX>::_inorder(pNode cur, Fun visit) const {
	if (cur != nullptr) {
		_inorder(cur->lchild, visit);
		visit(cur);
		_inorder(cur->rchild, visit);
	}
}

template<class BOX>
void BBTree<BOX>::InOrder(pFun_BBTree visit, utPointer utp) {
	_inorder(_root->lchild , visit, utp);
}
template<class BOX>
void BBTree<BOX>::InOrder(Fun visit) {
	_inorder(_root->lchild , visit);
}
template<class BOX>
void BBTree<BOX>::InOrder(Fun visit) const {
	_inorder(_root->lchild , visit);
}

template<typename BOX>
void BBTree<BOX>::_printtree(BBNode<BOX> *cur, int layer) {
	int i;
	if (cur == NULL) {
		return;
	}
	_printtree(cur->rchild, layer + 1);
	for (i = 0; i < layer; i++) {
		if (layer == 1) {
			std::cout << "  R>";
		} else {
			std::cout << "   ";
		}
	}
	std::cout << "O" << '\n';
	_printtree(cur->lchild, layer + 1);
}

template<typename BOX>
void BBTree<BOX>::print_tree() {
	_return_if_fail(_root != nullptr);
	_printtree(this->_root, 1);
}

template<class BOX>
BBTree<BOX>::BBTree(const Set<BOX>& set_box) {
	_root = nullptr;
	List<pNode> lc;
	for (auto iter = set_box.begin(); iter != set_box.end(); ++iter) {
		lc.push_back(new BBNode<BOX>((*iter)));
	}
	List<pNode> lf;
	Int flag = 1;
	List<pNode>* lcr = &lc;
	List<pNode>* lfr = &lf;
	//int idx = 0;
	while (!lcr->empty()) {
		bottomup_build((*lcr), (*lfr), flag);
		List<pNode>* tmpr = lcr;
		lcr = lfr;
		lfr = tmpr;
	}
}

template<class BOX>
void BBTree<BOX>::bottomup_build(List<pNode>& lc, List<pNode>& lf, Int& flag) {
	lf.clear();
	if (lc.size() == 1) {
		_root = new BBNode<BOX>();
		_root->lchild = (*(lc.begin()));
		_root->lchild->father = _root;
		return;
	}
	if (flag == 1) {
		for (auto iter = lc.begin(); iter != lc.end(); ++iter, ++iter) {
			auto iterp = iter;
			iterp++;
			if (iterp != lc.end()) {
				BOX _box((*iter)->_box, (*iterp)->_box);
				BBNode<BOX>* fnode = new BBNode<BOX>(_box, *iter, *iterp,
						nullptr);
				lf.push_back(fnode);
				(*iter)->father = fnode;
				(*iterp)->father = fnode;
			} else {
				flag = 0;
				lf.push_back((*iter));
				return;
			}
		}
	} else {
		for (auto iter = lc.rbegin(); iter != lc.rend(); ++iter, ++iter) {
			auto iterp = iter;
			iterp++;
			if (iterp != lc.rend()) {
				BOX _box((*iter)->_box, (*iterp)->_box);
				BBNode<BOX>* fnode = new BBNode<BOX>(_box, *iter, *iterp,
						nullptr);
				lf.push_back(fnode);
				(*iter)->father = fnode;
				(*iterp)->father = fnode;
			} else {
				flag = 1;
				lf.push_back((*iter));
				return;
			}
		}
	}
}
// cb call back
template<class NODE>
void cb_output_vtk_height(NODE* pn, utPointer utp) {

	Array<utPointer, 2> &arr = (*((Array<utPointer, 2>*) utp));
	Int& hei = (*((Int*) arr[0]));
	List<NODE*>* lbn = (List<NODE*>*) arr[1];
	if (pn->height() == hei) {
		lbn->push_back(pn);
	}
}

template<class BOX>
void BBTree<BOX>::output_vtk_height(const String& fn, Int hei) {
	_return_if_fail(hei >= 0);
	List<pNode> lbn;
	Array<utPointer, 2> arr;
	arr[0] = &hei;
	arr[1] = &lbn;
	InOrder(cb_output_vtk_height, &arr);
	FILE *data;
	data = fopen(fn.c_str(), "w");
	if (data == NULL) {
		std::cerr << "!> Open file error! " << fn << " \n";
		exit(-1);
	}
	uInt num_b = lbn.size();
	uInt NUM_VERTEXES = BOX::NumVertexes;
	fprintf(data, "# vtk DataFile Version 3.0\n");
	fprintf(data, "Gird output\n");
	fprintf(data, "ASCII\n");
	fprintf(data, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(data, "POINTS %d float\n", NUM_VERTEXES * num_b);
	for (auto iter = lbn.begin(); iter != lbn.end(); ++iter) {
		for (uInt i = 0; i < NUM_VERTEXES; i++) {
			fprintf(data, "%f %f %f \n",
					(*iter)->_box.get_point(_ORDER_VTK[i][0], _ORDER_VTK[i][1],
							_ORDER_VTK[i][2]).x(),
					(*iter)->_box.get_point(_ORDER_VTK[i][0], _ORDER_VTK[i][1],
							_ORDER_VTK[i][2]).y(),
					BOX::Dim == 3 ?
							(*iter)->_box.get_point(_ORDER_VTK[i][0],
									_ORDER_VTK[i][1], _ORDER_VTK[i][2]).z() :
							0);
		}
	}
	fprintf(data, "\n");
	fprintf(data, "CELLS %d %d \n", num_b, num_b * (NUM_VERTEXES + 1));
	for (int i = 0; i < NUM_VERTEXES * num_b; ++i) {
		if (i % NUM_VERTEXES == 0) {
			fprintf(data, "%d ", NUM_VERTEXES);
		}
		fprintf(data, "%d ", i);
		if (i % NUM_VERTEXES == 7) {
			fprintf(data, "\n");
		}
	}
	fprintf(data, "\n\n");
	fprintf(data, "CELL_TYPES %d\n", num_b);
	for (int i = 0; i < num_b; ++i) {
		fprintf(data, "%d \n", BOX::Dim == 3 ? 11 : 8);
	}
	//VTK_VOXEL (=11)  -- 3D
	//VTK_PIXEL (=8)   -- 2D
	fclose(data);

}

}

#endif /* BBTREE_H_ */
