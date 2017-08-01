
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
#include <list>
#include <algorithm>

namespace carpio{

template<class TYPE, St DIM>
class TriSurface_ {
public:
	static const St Dim = DIM;
	typedef TriSurface_<TYPE, DIM> self_class;
	typedef St size_type;
	typedef TYPE vt;
	typedef Point_<TYPE, DIM> Poi;
	typedef Poi* pPoi;
	typedef Vertex_<TYPE, DIM> Ver;
	typedef Ver* pVer;
	typedef Edge_<TYPE, DIM> Edg;
	typedef Edg* pEdg;
	typedef TriFace_<TYPE, DIM> Fac;
	typedef Fac* pFac;
	typedef TriSurface_<TYPE, DIM> Sur;
	typedef std::shared_ptr<Sur> pSur;
	typedef std::list<pVer> list_pVer;
	typedef std::list<pFac> list_pFac;
	typedef std::list<pSur> list_pSur;
	typedef std::function<void(Fac&)> Fun_Fac;
public:
	std::set<pFac> faces;
	std::set<pVer> c_vertex;
	std::set<pEdg> c_edge;
public:
	TriSurface_() {
	}
	TriSurface_(const std::string& filename) {
		this->load_gts_file(filename);
	}
	void load_gts_file(const std::string& filename);

	~TriSurface_() {
		faces.clear();
		//clear();
	}

	void clear() {
		//for (auto iter = c_vertex.begin(); iter != c_vertex.end(); ++iter) {
		//delete (*iter);
		//}
		//c_vertex.clear();
		//for (auto iter = c_edge.begin(); iter != c_edge.end(); ++iter) {
		//delete (*iter);
		//}
		//c_edge.clear();
		//for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
		//delete (*iter);
		//}
		//faces.clear();
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
						pEdg come = Tri::GetCommon_pEdge(v, nei);
						res = Tri::AreCompatible(v, nei, come);
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
				pVer ver = pf->get_vertex(i);
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
				pVer ver = pf->get_vertex(i);
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
			if (f1 != nullptr && f2 != nullptr && !AreCompatible(f1.get(), f2.get(), &edg)) {
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
	typename std::set<pFac>::iterator begin_face() {
		return this->faces.begin();
	}
	typename std::set<pFac>::const_iterator begin_face() const {
		return this->faces.begin();
	}
	typename std::set<pFac>::iterator end_face() {
		return this->faces.end();
	}
	typename std::set<pFac>::const_iterator end_face() const {
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

};

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
				pVer pver(new Ver(x, y, z));
				v_vertex.push_back(pver);
				v_vc.push_back(0);
			} else if (n <= nv + ne) {
				std::istringstream istr(sline);
				uInt i_v1, i_v2;
				istr >> i_v1 >> i_v2;
				pEdg pedg(new Edg(v_vertex[i_v1 - 1], v_vertex[i_v2 - 1]));
				pedg->attach();
				v_vc[i_v1 - 1]++;
				v_vc[i_v2 - 1]++;
				v_edge.push_back(pedg);
				v_ec.push_back(0);
			} else if (n <= nv + ne + nf) {
				std::istringstream istr(sline);
				uInt i_e1, i_e2, i_e3;
				istr >> i_e1 >> i_e2 >> i_e3;
				pFac pfac(
						new Fac(v_edge[i_e1 - 1], v_edge[i_e2 - 1],
								v_edge[i_e3 - 1], this));
				pfac->attach();
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
		faces.insert(ptr);
	}
	//for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
	//	auto ptr = (*iter);
	//	ptr->show();
	//}

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
	FILE* fptr = fopen(fn.c_str(), "w"); //write
	if (fptr == NULL) {
		std::cerr << "!> Open file error! " << fn << " \n";
		exit(-1);
	}
	fprintf(fptr, "# vtk DataFile Version 2.0\n"
			"Generated by LarusTS\n"
			"ASCII\n"
			"DATASET POLYDATA\n"
			"POINTS %lu float\n", size_vertex());
	std::map<pVer, uInt> m_veridx;
	uInt count = 0;
	// not ok
	for (auto iter = c_vertex.begin(); iter != c_vertex.end(); ++iter) {
		auto pt = (*iter);
		fprintf(fptr, "%f %f %f \n", pt->x(), pt->y(), pt->z());
		m_veridx.insert(std::pair<pVer, uInt>(pt, count));
		count++;
	}
	fprintf(fptr, "POLYGONS %lu %lu\n", faces.size(), faces.size() * 4);
	for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
		auto pt = (*iter);
		uInt i1, i2, i3;
		i1 = m_veridx.find(pt->e1->v1)->second;
		i2 = m_veridx.find(pt->e1->v2)->second;
		if (pt->e2->v2 != pt->e1->v2 && pt->e2->v2 != pt->e1->v1) {
			i3 = m_veridx.find(pt->e2->v2)->second;
		} else {
			i3 = m_veridx.find(pt->e2->v1)->second;
		}
		fprintf(fptr, "3 %u %u %u\n", i1, i2, i3);
	}
	fclose(fptr);
}

}

#endif /* _SURFACE_H_ */
