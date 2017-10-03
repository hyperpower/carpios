/************************
 //  \file   ts_surface.h
 //  \brief
 //
 //  \author czhou
 //  \date   15 juin 2015
 ***********************/
#ifndef _SURFACE_HPP_
#define _SURFACE_HPP_

#include "_edge.hpp"
#include "_vertex.hpp"
#include "_triface.hpp"
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <algorithm>
#include <memory>

namespace carpio {

template<class TYPE, St DIM>
class TriSurface_ {
public:
	static const St Dim = DIM;
	typedef TriSurface_<TYPE, DIM> Self;
	typedef St size_type;
	typedef TYPE vt;
	typedef TriFace_<TYPE, DIM, Self> Fac;
	typedef Edge_<TYPE, DIM, Fac> Edg;
	typedef Point_<TYPE, DIM> Poi;
	typedef Poi* pPoi;
	typedef Vertex_<TYPE, DIM, Edg> Ver;
	typedef Ver* pVer;
	typedef Edg* pEdg;
	typedef Fac* pFac;
	typedef const Ver* const_pVer;
	typedef const Edg* const_pEdg;
	typedef const Fac* const_pFac;
	typedef TriSurface_<TYPE, DIM> Sur;
	typedef Sur* pSur;
	typedef std::list<pVer> list_pVer;
	typedef std::list<pFac> list_pFac;
	typedef std::list<pSur> list_pSur;
	typedef std::function<void(Fac&)> Fun_Fac;

	typedef typename std::list<pFac>::iterator iterator;
	typedef typename std::list<pFac>::const_iterator const_iterator;
	public:
	std::set<pFac> faces;

	Any _any_data;
	//std::set<pEdg> c_edge;
	//std::set<pVer> c_vertex;
public:
	TriSurface_() {
	}
	TriSurface_(const std::string& filename) {
		this->load_gts_file(filename);
	}
	void load_gts_file(const std::string& filename);

	~TriSurface_() {
		//faces.clear();
		clear();
	}

	void copy_to(Self& other) {
		other.clear();
		for (auto& pf : faces) {
			other.insert(pf);
		}
	}

	void insert(pFac f) {
		auto r = std::find(faces.begin(), faces.end(), f);
		if (r == faces.end()) { //not found
			faces.insert(f);
			f->attach(this);
		} else {
			return;
		}
	}

	void erase(pFac f) {
		_IF_TRUE_RETRUN(f == nullptr);
		std::list<pFac> ldfac;
		std::list<pEdg> ldedg;
		std::list<pVer> ldver;
		auto iter = std::find(faces.begin(), faces.end(), f);
		if (iter != faces.end()) { // found
			pFac pf = (*iter);
			if (pf->has_one_parent_surface(this)) {
				ldfac.push_back(pf);
			}
			pf->detach(this);
			this->faces.erase(pf);
			if (pf->is_unattached()) {
				for (int i = 0; i < 3; i++) {
					pEdg e = pf->edge(i);
					if (e->has_one_face(pf)) {
						ldedg.push_back(e);
					}
					e->detach(pf);
					///
					if (e->is_unattached()) {
						for (int iv = 0; iv < 2; iv++) {
							pVer v = e->vertex(iv);
							if (v->has_one_edge(e)) {
								ldver.push_back(v);
							}
							v->detach(e);
						}
					}
				}
			}
		} else {
			return;
		}
		for (auto& fac : ldfac) {
			delete fac;
		}
		for (auto& edg : ldedg) {
			delete edg;
		}
		for (auto& ver : ldver) {
			delete ver;
		}

	}

	void clear() {
		std::list<pFac> ldfac;
		std::list<pEdg> ldedg;
		std::list<pVer> ldver;
		for (typename std::set<pFac>::iterator iter = faces.begin();
				iter != faces.end(); ++iter) {
			pFac pf = (*iter);
			if (pf->has_one_parent_surface(this)) {
				ldfac.push_back(pf);
			}
			pf->detach(this);
			if (pf->is_unattached()) {
				for (int ie = 0; ie < 3; ie++) {
					pEdg e = pf->edge(ie);
					if (e->has_one_face(pf)) {
						ldedg.push_back(e);
					}
					e->detach(pf);
					///
					for (int iv = 0; iv < 2; iv++) {
						pVer v = e->vertex(iv);
						if (v->has_one_edge(e)) {
							ldver.push_back(v);
						}
						v->detach(e);
					}
				}
			}
		}
		for (auto& fac : ldfac) {
			this->faces.erase(fac);
			delete fac;
		}
		for (auto& edg : ldedg) {
			delete edg;
		}
		for (auto& ver : ldver) {
			delete ver;
		}

	}

	void clear2() {
//		for (auto iter = c_vertex.begin(); iter != c_vertex.end(); ++iter) {
//			delete (*iter);
//		}
//		c_vertex.clear();
//		for (auto iter = c_edge.begin(); iter != c_edge.end(); ++iter) {
//			delete (*iter);
//		}
//		c_edge.clear();
		for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
			delete (*iter);
		}
		faces.clear();
	}
	size_type size_edge() const {
		return count_edge();   // can be improved
	}
	size_type size_vertex() const {
		return count_vertex(); // can be improved
	}
	size_type count_edge() const {
		size_type count_edge = 0;
		std::function<void(Sur::Edg&)> fun2 = [&count_edge](Sur::Edg&) {
			count_edge++;
		};
		this->foreach_edge(fun2);
		return count_edge;
	}
	size_type count_vertex() const {
		size_type count = 0;
		std::function<void(Sur::Ver&)> fun = [&count](Sur::Ver&) {
			count++;
		};
		this->foreach_vertex(fun);
		return count;
	}
	size_type size_face() const {
		return faces.size();
	}

	bool empty() const {
		return faces.empty();
	}
	void foreach_face(std::function<void(Fac&)> fun) {
		for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
			pFac pf = (*iter);
			fun(*pf);
		}
	}

	void foreach_face(std::function<void(Fac&)> fun) const {
		for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
			pFac pf = (*iter);
			fun(*pf);
		}
	}

	void walk_each_face(std::function<void(Fac&)> fun) {
		auto iterb = faces.begin();
		pFac spf = (*iterb);
		std::list<pFac> qold, qnew;
		qnew.push_back(spf);
		while (!qnew.empty()) {
			pFac v = qnew.front();
			auto r = std::find(qold.begin(), qold.end(), v);
			if (r == qold.end()) { //not found
				fun(*v);
				qold.push_back(v);
				// found its neighbor
				list_pFac ln = v->get_neighbor_faces(this);
				for (auto i = ln.begin(); i != ln.end(); i++) {
					if (std::find(qnew.begin(), qnew.end(), *i) == qnew.end()) {
						qnew.push_back(*i);
					}
				}
			}
			qnew.pop_front();
		}

	}

	void walk_each_face(std::function<void(Fac&)> fun) const {
		auto iterb = faces.begin();
		pFac spf = (*iterb);
		std::list<pFac> qold, qnew;
		qnew.push_back(spf);
		while (!qnew.empty()) {
			pFac v = qnew.front();
			auto r = std::find(qold.begin(), qold.end(), v);
			if (r == qold.end()) { //not found
				fun(*v);
				qold.push_back(v);
				// found its neighbor
				list_pFac ln = v->get_neighbor_faces(this);
				for (auto i = ln.begin(); i != ln.end(); i++) {
					if (std::find(qnew.begin(), qnew.end(), *i) == qnew.end()) {
						qnew.push_back(*i);
					}
				}
			}
			qnew.pop_front();
		}
	}

	bool is_compatible() const {
		bool res = true;
		auto iterb = faces.begin();
		pFac spf = (*iterb);
		std::list<pFac> qold, qnew;
		qnew.push_back(spf);
		while (!qnew.empty()) {
			pFac v = qnew.front();
			auto r = std::find(qold.begin(), qold.end(), v);
			if (r == qold.end()) { //not found
				qold.push_back(v);
				// found its neighbor
				list_pFac ln = v->get_neighbor_faces(this);
				for (auto i = ln.begin(); i != ln.end(); i++) {
					if (std::find(qnew.begin(), qnew.end(), *i) == qnew.end()) {
						pFac nei = *i;
						pEdg come = Fac::GetCommon_pEdg(v, nei);
						res = Fac::AreCompatible(v, nei, come);
						if (res == false) {
							return res;
						}
						qnew.push_back(nei);
					}
				}
			}
			qnew.pop_front();
		}
		return res;
	}

	void revert() {
		for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
			pFac pf = (*iter);
			pf->revert();
		}
	}

	void foreach_edge(std::function<void(Edg&)> fun) {
		std::set<pEdg> tmp;
		for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
			pFac pf = (*iter);
			for (St i = 0; i < 3; i++) {
				pEdg edg = pf->operator[](i);
				if (tmp.find(edg) == tmp.end()) {
					fun(*edg);
				}
				tmp.insert(edg);
			}
		}
	}
	void foreach_edge(std::function<void(Edg&)> fun) const {
		std::set<pEdg> tmp;
		for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
			pFac pf = (*iter);
			for (St i = 0; i < 3; i++) {
				pEdg edg = pf->operator[](i);
				if (tmp.find(edg) == tmp.end()) {
					fun(*edg);
				}
				tmp.insert(edg);
			}
		}
	}

	void foreach_vertex(std::function<void(Ver&)> fun) {
		std::set<pVer> tmp;
		for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
			pFac pf = (*iter);
			for (St i = 0; i < 3; i++) {
				pVer ver = pf->vertex(i);
				if (tmp.find(ver) == tmp.end()) {
					fun(*ver);
				}
				tmp.insert(ver);
			}
		}
	}
	void foreach_vertex(std::function<void(Ver&)> fun) const {
		std::set<pVer> tmp;
		for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
			pFac pf = (*iter);
			for (St i = 0; i < 3; i++) {
				pVer ver = pf->vertex(i);
				if (tmp.find(ver) == tmp.end()) {
					fun(*ver);
				}
				tmp.insert(ver);
			}
		}
	}

	void transfer(vt dx, vt dy, vt dz = 0) {
		std::function<void(Sur::Ver&)> fun = [&dx, &dy, &dz](Sur::Ver& ver) {
			ver[0] += dx;
			ver[1] += dy;
			if (Ver::Dim == 3) {
				ver[2] += dz;
			}
		};
		foreach_vertex(fun);
	}

	bool is_orientable() const {
		bool res = true;
		std::function<void(Sur::Edg&)> fun = [&res, this](Sur::Edg& edg) {
			pFac f1 = nullptr, f2 = nullptr;
			int i = 0;
			for (pFac pf : edg.faces) {
				//pFac pf = dynamic_cast<pFac>(pt);
				if (pf != nullptr && pf->has_parent_surface(this)) {
					if (f1 == nullptr) {
						f1 = pf;
					} else if (f2 == nullptr) {
						f2 = pf;
					} else {
						res = false;
					}
				}
				i++;
			}
			if (f1 != nullptr && f2 != nullptr && !Fac::AreCompatible(f1, f2, &edg)) {
				res = false;
				//std::cout<<"size "<<pe->faces.size()<<endl;
				//f2->show();
				//std::cout<<"compatible\n";
			}
		};
		foreach_edge(fun);
		return res;
	}

	void add_face(const pFac& spf) {
		pFac f = spf;
		f->attach_surface(this);
		faces.insert(spf);
	}

// iterator
	typename std::set<pFac>::iterator begin() {
		return this->faces.begin();
	}
	typename std::set<pFac>::const_iterator begin() const {
		return this->faces.begin();
	}
	typename std::set<pFac>::iterator end() {
		return this->faces.end();
	}
	typename std::set<pFac>::const_iterator end() const {
		return this->faces.end();
	}
// connection -------------------------------
// check ------------------------------------

// show -------------------------------------
	void show_vertex() const;
	void show_edge() const;
	void show_face() const;

// output -----------------------------------
	void output_vtk(const std::string& fn) const;

	static void Copy(pSur src, pSur dst) {
		src->copy_to(*dst);
	}
	static void Copy(Sur& src, Sur& dst) {
		src.copy_to(dst);
	}

	/// get all the faces to vertex
	static void FacesConnectTo(pVer v, pSur sur, std::set<pFac>& res) {
		for (auto iter = v->begin_edge(); iter != v->end_edge(); ++iter) {
			pEdg pe = *iter;
			for (auto itface = pe->begin_face(); itface != pe->end_face();
					++itface) {
				pFac pf = *itface;
				if (pf->has_parent_surface(sur)) {
					res.insert(pf);
				}
			}
		}
	}


	template<class Container>
	static std::set<pFac> FacesConnectTo(const Container& con, pSur sur) {
		typename Container::value_type dummy;
		return _FacesConnectTo(con, sur, dummy);
	}

	template<class Container>
	static std::set<pFac> _FacesConnectTo(
			Container& con,
			pSur sur,
			pVer dummy) {
		std::set<pFac> res;
		for (auto& pv : con) {
			FacesConnectTo(pv, sur, res);
		}
		return res;
	}

	static pFac Neighbor(
			pFac f,
			pEdg e,
			Sur& surface) {
		for (auto itf = e->begin_face(); itf != e->end_face(); ++itf) {
			pFac pf = *itf;
			if (pf != f && pf->has_parent_surface(&surface)) {
				return pf;
			}
		}
		return nullptr;
	}

}
;

template<class TYPE, St DIM>
void TriSurface_<TYPE, DIM>::load_gts_file(const std::string& filename) {
	this->faces.clear();
	std::ifstream fs;
	fs.open(filename.c_str(), std::fstream::in); //read
	if (!fs.is_open()) {
		std::cerr << "!> Can't find file. " << filename.c_str() << " \n";
		exit(-1);
	}
	uInt n = 0;
	uInt nv, ne, nf;
	std::vector<pVer> v_vertex;
	std::vector<Int> v_vc;
	std::vector<pEdg> v_edge;
	std::vector<Int> v_ec;
	std::vector<pFac> v_face;
	//pSur psur = this->get_this();
	while (!fs.eof()) {
		std::string sline;
		getline(fs, sline, '\n');
		if (sline.length() != 0) {
			if (n == 0) {
				std::istringstream istr(sline);
				istr >> nv >> ne >> nf;
			} else if (n <= nv) {
				std::istringstream istr(sline);
				TYPE x, y, z;
				istr >> x >> y >> z;
				pVer pver = new Ver(x, y, z);
				v_vertex.push_back(pver);
				v_vc.push_back(0);
			} else if (n <= nv + ne) {
				std::istringstream istr(sline);
				uInt i_v1, i_v2;
				istr >> i_v1 >> i_v2;
				pEdg pedg = new Edg(v_vertex[i_v1 - 1], v_vertex[i_v2 - 1]);
				v_vc[i_v1 - 1]++;
				v_vc[i_v2 - 1]++;
				v_edge.push_back(pedg);
				v_ec.push_back(0);
			} else if (n <= nv + ne + nf) {
				std::istringstream istr(sline);
				uInt i_e1, i_e2, i_e3;
				istr >> i_e1 >> i_e2 >> i_e3;
				pFac pfac = new Fac(v_edge[i_e1 - 1], v_edge[i_e2 - 1],
						v_edge[i_e3 - 1], this);
				v_ec[i_e1 - 1]++;
				v_ec[i_e2 - 1]++;
				v_ec[i_e3 - 1]++;
				v_face.push_back(pfac);
			}
		}
		n++;
	}
	//insert pface to set
	for (auto iter = v_face.begin(); iter != v_face.end(); ++iter) {
		auto ptr = (*iter);
		this->insert(ptr);
	}
}
template<class TYPE, St DIM>
void TriSurface_<TYPE, DIM>::show_vertex() const {
	std::cout << " size vertex = " << this->size_vertex() << "\n";
	std::function<void(Sur::Ver&)> fun = [](Sur::Ver& ver) {
		ver.show();
	};
	this->foreach_vertex(fun);
}
template<class TYPE, St DIM>
void TriSurface_<TYPE, DIM>::show_edge() const {
	std::cout << " size edge = " << size_edge() << "\n";
	std::function<void(Sur::Edg&)> fun = [](Sur::Edg& edg) {
		edg.show();
	};
	this->foreach_edge(fun);
}

template<class TYPE, St DIM>
void TriSurface_<TYPE, DIM>::show_face() const {
	std::cout << " size face = " << faces.size() << "\n";
	for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
		(*iter)->show();
	}
}

template<class TYPE, St DIM>
void TriSurface_<TYPE, DIM>::output_vtk(const std::string& fn) const {
//	FILE* fptr = fopen(fn.c_str(), "w"); //write
//	if (fptr == NULL) {
//		std::cerr << "!> Open file error! " << fn << " \n";
//		exit(-1);
//	}
//	fprintf(fptr, "# vtk DataFile Version 2.0\n"
//			"Generated by LarusTS\n"
//			"ASCII\n"
//			"DATASET POLYDATA\n"
//			"POINTS %lu float\n", size_vertex());
//	std::map<pVer, uInt> m_veridx;
//	uInt count = 0;
//	// not ok
//	for (auto iter = c_vertex.begin(); iter != c_vertex.end(); ++iter) {
//		auto pt = (*iter);
//		fprintf(fptr, "%f %f %f \n", pt->x(), pt->y(), pt->z());
//		m_veridx.insert(std::pair<pVer, uInt>(pt, count));
//		count++;
//	}
//	fprintf(fptr, "POLYGONS %lu %lu\n", faces.size(), faces.size() * 4);
//	for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
//		auto pt = (*iter);
//		uInt i1, i2, i3;
//		i1 = m_veridx.find(pt->e1->v1)->second;
//		i2 = m_veridx.find(pt->e1->v2)->second;
//		if (pt->e2->v2 != pt->e1->v2 && pt->e2->v2 != pt->e1->v1) {
//			i3 = m_veridx.find(pt->e2->v2)->second;
//		} else {
//			i3 = m_veridx.find(pt->e2->v1)->second;
//		}
//		fprintf(fptr, "3 %u %u %u\n", i1, i2, i3);
//	}
//	fclose(fptr);
}

}

#endif /* _SURFACE_H_ */
