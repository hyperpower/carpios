#ifndef _PLOTLY_H_
#define _PLOTLY_H_

#include "io_define.hpp"
#include "algebra/array_list.hpp"
#include <Python.h>
#include <map>
#include <cmath>
#include <memory>
#include <algorithm>

namespace carpio {

class Plotly_actor {
public:
	typedef PyObject* pPO;
	typedef std::map<std::string, pPO> Map;
	typedef ArrayListV<double> Arrd;
	typedef std::list<double> Listd;
protected:
	Map _map;
	pPO _module;
public:
	Plotly_actor();
	virtual ~Plotly_actor();

	pPO get_p_python_object() {
		pPO dict = PyDict_New();
		for (Map::const_iterator iter = this->_map.begin();
				iter != this->_map.end(); iter++) {
			const std::string& skey = iter->first;
			pPO key = Py_BuildValue("s", skey.c_str());
			pPO val = iter->second;
			PyDict_SetItem(dict, key, val);
		}
		pPO args = Py_BuildValue("(O)", dict);
		return PyObject_CallObject(this->_module, args);
	}

	pPO get_module(std::string name) {
		pPO moduleName = PyUnicode_FromString("plotly.graph_objs");
		pPO imn = PyImport_Import(moduleName);
		pPO module = PyObject_GetAttrString(imn, name.c_str());
		ASSERT(!(!module || !PyCallable_Check(module)));
		return module;
	}

	void set_name(const std::string& name) {
		pPO pv = Py_BuildValue("s", name.c_str());
		this->_map["name"] = pv;
		//this->_map["showlegend"] = Py_BuildValue("p", 1);
	}

	void set_opacity(const double& val) {
		ASSERT(val >= 0 && val <= 1);
		pPO pv = Py_BuildValue("d", val);
		this->_map["opacity"] = pv;
	}

	void set_mode(const std::string& mode) {
		pPO pv = Py_BuildValue("s", mode.c_str());
		this->_map["mode"] = pv;
	}

protected:

	template <class Container>
	pPO _to_list(const Container& arr, int jump = 0) const {
		bool checkt = IsIterable<Container>::value;
		ASSERT(checkt);
		pPO pl = PyList_New(0);
		int i = 0;
		for (Listd::value_type v : arr) {
			pPO p = Py_BuildValue("d", v);
			if (jump > 0) {
				if (i % jump == 0 && i > 0) {
					PyList_Append(pl, Py_None);
				}
			}
			PyList_Append(pl, p);
			i++;
		}
		return pl;
	}

	pPO _to_list_list(const Listd& arr,
			const typename Listd::size_type& ndim1) {
		ASSERT(ndim1 > 0);
		pPO pl = PyList_New(0);
		int i = 0;
		pPO prow = PyList_New(0);
		for (Listd::value_type v : arr) {
			if (i % ndim1 != 0 || i == 0) {
				pPO p = PyFloat_FromDouble(v);
				PyList_Append(prow, p);
			} else {
				PyList_Append(pl, prow);
				prow = PyList_New(0);
				pPO p = PyFloat_FromDouble(v);
				PyList_Append(prow, p);
			}
			i++;
		}
		return pl;
	}

};

class Plotly_actor_scatter: public Plotly_actor {
public:
	typedef PyObject* pPO;
	typedef ArrayListV<double> Arrd;
	typedef std::list<double> Listd;
public:
	Plotly_actor_scatter(const Arrd& x, const Arrd& y, int jump = 0) {
		//Py_Initialize();
		_module = this->get_module("Scatter");
		pPO px = this->_to_list(x, jump);
		pPO py = this->_to_list(y, jump);
		this->_map["x"] = px;
		this->_map["y"] = py;
		//_data = nullptr;
	}
	Plotly_actor_scatter(const Listd& x, const Listd& y, int jump = 0) {
		//Py_Initialize();
		_module = this->get_module("Scatter");
		pPO px = this->_to_list(x, jump);
		pPO py = this->_to_list(y, jump);
		this->_map["x"] = px;
		this->_map["y"] = py;
		//_data = nullptr;
	}

	void set_colorscale_range(const double& minv, const double & maxv) {
		// make sure the minv and maxv is the real min and max, otherwise the function will not work.
		pPO dict = this->_map["marker"];
		pPO cs = Py_BuildValue("s", "color");
		pPO lv = PyDict_GetItem(dict, cs);
		pPO minp = Py_BuildValue("d", minv);
		PyList_Append(lv, minp);
		pPO maxp = Py_BuildValue("d", maxv);
		PyList_Append(lv, maxp);
		//PyDict_SetItem(dict, cs, lv);
		//this->_map["marker"] = dict;
	}

	void set_colorscale_name(const std::string& name) {
		pPO dict = this->_map["marker"];
		pPO key = Py_BuildValue("s", "colorscale");
		pPO val = Py_BuildValue("s", name.c_str());
		PyDict_SetItem(dict, key, val);
	}

	void set_colorscale(
			const Listd& d, //
			const int& size = 10, const std::string name = "Viridis",
			const double& minv = 0.0, const double& maxv = 0.0) {
		typedef typename Listd::value_type vt;
		pPO dict = PyDict_New();
		pPO key = Py_BuildValue("s", "color");
		//vt min = *std::min_element(d.begin(), d.end());
		//vt max = *std::max_element(d.begin(), d.end());
		//vt scale = 1 / min;
		//Listd nd;
		//std::for_each(d.begin(), d.end(), [&nd](vt i) {
		//std::cout<<-i<<" \n";
		//	nd.push_back(-i * 10000);
		//});
		pPO val = this->_to_list(d);
		if (!(minv == maxv && minv == 0.0)) {
			pPO minp = Py_BuildValue("d", minv);
			PyList_Append(val, minp);
			pPO maxp = Py_BuildValue("d", maxv);
			PyList_Append(val, maxp);
		}
		//
		PyDict_SetItem(dict, key, val);
		pPO cs = Py_BuildValue("s", "colorscale");
		pPO val2 = Py_BuildValue("s", name.c_str());
		PyDict_SetItem(dict, cs, val2);
		pPO key3 = Py_BuildValue("s", "size");
		pPO val3 = Py_BuildValue("i", size);
		PyDict_SetItem(dict, key3, val3);
		//pPO key4 = Py_BuildValue("s", "colorscale");
		//pPO val4 = Py_BuildValue("i", 1);
		//std::cout<<"val4 "<<PyObject_IsTrue(val4)<<"\n";
		//PyDict_SetItem(dict, key4, val4);
		this->_map["marker"] = dict;
	}

	void set_add_val(const Listd& d) {
		pPO val = this->_to_list(d);
		this->_map["text"] = val;
	}

};

class Plotly_actor_scatter3d: public Plotly_actor {
public:
	typedef PyObject* pPO;
	typedef ArrayListV<double> Arrd;
	typedef std::list<double> Listd;
public:
	template<class Container>
	Plotly_actor_scatter3d(
			const Container& x,
			const Container& y,
			const Container& z,
			int jump = 0) {
		bool checkt = IsIterable<Container>::value;
		//Py_Initialize();
		_module = this->get_module("Scatter3d");
		pPO px = this->_to_list(x, jump);
		pPO py = this->_to_list(y, jump);
		pPO pz = this->_to_list(z, jump);
		this->_map["x"] = px;
		this->_map["y"] = py;
		this->_map["z"] = pz;
	}

	void set_colorscale(const Listd& d, const int& size = 10,
			const std::string name = "Viridis") {
		typedef typename Listd::value_type vt;
		pPO dict = PyDict_New();
		pPO key = Py_BuildValue("s", "color");
		pPO val = this->_to_list(d);
		PyDict_SetItem(dict, key, val);
		pPO key2 = Py_BuildValue("s", "colorscale");
		pPO val2 = Py_BuildValue("s", name.c_str());
		PyDict_SetItem(dict, key2, val2);
		pPO key3 = Py_BuildValue("s", "size");
		pPO val3 = Py_BuildValue("i", size);
		PyDict_SetItem(dict, key3, val3);
		pPO key4 = Py_BuildValue("s", "showscale");
		pPO val4 = Py_BuildValue("i", size);
		PyDict_SetItem(dict, key4, val4);
		this->_map["marker"] = dict;
	}

	void set_add_val(const Listd& d) {
		pPO val = this->_to_list(d);
		this->_map["text"] = val;
	}

};

class Plotly_actor_mesh3d: public Plotly_actor {
public:
	typedef PyObject* pPO;
	typedef ArrayListV<double> Arrd;
public:
	template<class Container>
	Plotly_actor_mesh3d(const Container& x, const Container& y, const Container& z) {
		bool checkt = IsIterable<Container>::value;
		ASSERT(checkt);
		_module = this->get_module("Mesh3d");
		pPO px = this->_to_list(x);
		pPO py = this->_to_list(y);
		pPO pz = this->_to_list(z);
		this->_map["x"] = px;
		this->_map["y"] = py;
		this->_map["z"] = pz;
		//_data = nullptr;
	}

	void set_ijk(const Listd& x, const Listd& y, const Listd& z) {
		pPO px = this->_to_list(x);
		pPO py = this->_to_list(y);
		pPO pz = this->_to_list(z);
		this->_map["i"] = px;
		this->_map["j"] = py;
		this->_map["k"] = pz;
	}
};

class Plotly_actor_heatmap: public Plotly_actor {
public:
	typedef PyObject* pPO;
	typedef ArrayListV<double> Arrd;
	typedef std::list<double> Listd;
public:
	// x, y indicate the location of the faces
	Plotly_actor_heatmap(const Listd& x, const Listd& y, const Listd& z) {
		//Py_Initialize();
		_module = this->get_module("Heatmap");
		pPO px = this->_to_list(x);
		pPO py = this->_to_list(y);
		pPO pz = this->_to_list_list(z, x.size() - 1);
		this->_map["x"] = px;
		this->_map["y"] = py;
		this->_map["z"] = pz;
	}

	void colorscale(const std::string& cs = "Vidridis") {
		pPO val = Py_BuildValue("s", cs.c_str());
		this->_map["colorscale"] = val;
	}

};

class Plotly {
public:
	typedef PyObject* pPO;
	typedef std::shared_ptr<Plotly_actor> spPA;
	typedef std::map<std::string, pPO> Map;
public:
	Plotly() {
		_init();
	}

	~Plotly() {
		//Py_Finalize();
	}

	std::string version() const {
		pPO version = PyObject_GetAttrString(_modulePlotly, "__version__");
		//std::string a(PyString_AsString(version));
		std::string a("ooooo");
		return a;
	}

	void add(spPA spa) {
		this->_actors.push_back(spa);
	}

	void plot() const {
		pPO data = PyList_New(0);
		//pPO list_actors;
		for (spPA pa : _actors) {
			pPO p = pa->get_p_python_object();
			PyList_Append(data, p);
		}
		if (PyList_Size(data) == 0) {
			return;
		}
		//
		pPO layout = _get_layout();
		pPO args;
			args = Py_BuildValue("({s:O, s:O})", "data", data, "layout",
					layout);
		pPO pRet = PyObject_CallObject(_funPlot, args);
	}

	void title(const std::string& t) {
		pPO val = Py_BuildValue("s", t.c_str());
		_map_layout["title"] = val;
	}

	void width(const int& w) {
		ASSERT(w > 10);
		pPO val = Py_BuildValue("i", w);
		_map_layout["width"] = val;
	}
	void height(const int& w) {
		ASSERT(w > 10);
		pPO val = Py_BuildValue("i", w);
		_map_layout["height"] = val;
	}
	void size(const int& w, const int& h) {
		//this->disable_autosize();
		this->width(w);
		this->height(h);
	}

	//void set_output_file(const std::string& fn) {
	//	this->_filename = Py_BuildValue("s", fn.c_str());
	//}

protected:
// data
	pPO _modulePlotly;
	pPO _funPlot;
	//pPO _filename;

	Map _map_layout;

	std::list<spPA> _actors;

	pPO _get_layout() const {
		pPO dict = PyDict_New();
		for (Map::const_iterator iter = this->_map_layout.begin();
				iter != this->_map_layout.end(); iter++) {
			const std::string& skey = iter->first;
			pPO key = Py_BuildValue("s", skey.c_str());
			pPO val = iter->second;
			PyDict_SetItem(dict, key, val);
		}
		pPO args = Py_BuildValue("(O)", dict);
		pPO _module = _get_module("Layout");
		return PyObject_CallObject(_module, args);
	}

	pPO _get_module(std::string name) const {
		pPO moduleName = PyUnicode_FromString("plotly.graph_objs");
		pPO imn = PyImport_Import(moduleName);
		pPO module = PyObject_GetAttrString(imn, name.c_str());
		ASSERT(!(!module || !PyCallable_Check(module)));
		return module;
	}

	void _init() {
		Py_Initialize();
		pPO moduleName = PyUnicode_FromString("plotly");
		_modulePlotly = PyImport_Import(moduleName);

		moduleName = PyUnicode_FromString("plotly.offline");
		pPO offline = PyImport_Import(moduleName);
		_funPlot = PyObject_GetAttrString(offline, "plot");

		//_filename = nullptr;
		ASSERT(!(!_funPlot || !PyCallable_Check(_funPlot)));
	}

}
;

}

#endif
