#ifndef _S_EVENT_HPP
#define _S_EVENT_HPP

#include "s_define.hpp"
#include <map>
#include <utility>
#include "utility/clock.h"

namespace structure {

/*
 * The event class
 * This is the base class
 */
template<St DIM> class Equation_;

template<St Dim> class Operation_;

template<St DIM>
class Event_ {
public:
	static const St Dim = DIM;
	static const St NumFaces = DIM + DIM;
	static const St NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);

	typedef Equation_<Dim> Equ;
	typedef Equ* pEqu;
	typedef const Equ* const_pEqu;

	static const int START = 0x10;
	static const int END = 0x20;
	static const int BEFORE = 0x100;
	static const int AFTER = 0x200;

protected:
	/*
	 *  _flag = 0      disabled flag
	 *  _flag = -1     forward
	 *  _flag = 1      backward      (default)
	 *  _flag = 2      both run, forward and backward
	 *  _flag = -2     both do not run, forward and backward
	 */
	int _flag;
	int _istart, _iend, _istep;

	bool _has_flag(int f) const {
		return (_flag | f) == _flag ? true : false;
	}

public:
	/**
	 *  @brief Event constructor
	 *
	 *  default constructor is running event after each step
	 */

	Event_(int is = -1, int ie = -1, int istep = -1, int flag = 0) {
		_flag = flag;
		_istart = is;
		_iend = ie;
		_istep = istep;
	}

	virtual int execute(St step, Vt t, int fob, pEqu pd = nullptr) {
		std::cout << "Event : execute \n";
		return -1;
	}

	virtual int flag() const {
		return 0;
	}

	virtual void set_flag(int i) {
		return;
	}

	//virtual DATA data() const {
	//	return 0;
	//}

	//virtual void set_data(DATA v) {
	//	return;
	//}

	virtual ~Event_() {
	}

	/**
	 *  @brief Dose the event need to be executed
	 *
	 *  @param [step] is step right now
	 *  @param [t]    is time right now
	 *  @param [fob]  forward or backward,
	 *                fob = -1 is before the step,
	 *                fob = 1 is after the step
	 */
	bool _do_execute(St step, Vt t, int fob) const {
		if (this->_has_flag(fob)) {
			bool res = ((step - this->_istart) % this->_istep == 0);
			if (this->_iend != -1) {
				res = res && int(step) <= this->_iend;
			}
			if (this->_istart != -1) {
				res = res && int(step) >= this->_istart;
			}
			// self value
			//  &&
			return res;
		} else {
			return false;
		}
	}

};

/*
 * derived class
 */
template<St DIM>
class EventOutput_: public Event_<DIM> {
public:
	typedef Event_<DIM> Event;
protected:
	std::FILE* _f;
	int _flag_output;
public:
	EventOutput_(int is, int ie, int istep, int flag, std::string fn = "") :
			Event(is, ie, istep, flag) {
		if (fn == "") {
			_f = nullptr;
			_flag_output = 0;
		} else {
			_f = std::fopen(fn.c_str(), "w");
			if (!_f) {
				std::perror("File opening failed");
			}
			_flag_output = 1;
		}
	}

	int flag() const {
		return _flag_output;
	}

	std::string name() const {
		return "Output";
	}

	virtual ~EventOutput_() {
		if (_f != nullptr) {
			std::fclose(_f);
		}
	}

};
/*
 template<St DIM>
 class EventTraceCenterValue_: public EventOutput_<DIM> {
 public:
 typedef EventOutput_<DIM> EventOutput;

 typedef Equation_<DIM> Equ;
 typedef Equ* pEqu;
 typedef const Equ* const_pEqu;
 protected:
 ArrayListV<cvt> _pos;
 St _idx;

 std::list<vt> _l_t;     //time
 std::list<St> _l_step;  //step
 std::list<tick_t> _l_v;   //value
 public:
 EventTraceCenterValue_(int is = -1, int ie = -1, int istep = -1, int flag =
 0, std::FILE* f = nullptr) :
 EventOutput(is, ie, istep, flag, f), _pos(3), _idx(0) {
 _pos.zeros();

 }
 void set_location(cvt x, cvt y = 0, cvt z = 0) {
 _pos[0] = x;
 _pos[1] = y;
 _pos[2] = z;
 }
 void set_value_idx(St i) {
 this->_idx = i;
 }
 int execute(St step, vt t, int fob, pDomain pd = nullptr) {
 _record_value(step, t, pd);
 return -1;
 }

 protected:
 void _record_value(const St& step, vt& t, pDomain pd) {
 _l_step.push_back(step);
 _l_t.push_back(t);
 pNode pn = pd->grid().get_pnode(_pos[0], _pos[1], _pos[2]);
 vt value = pn->cda(_idx);
 _l_v.push_back(value);
 std::cout << "Trace Center Value  " << value << "\n";
 }
 };
 */
template<St DIM>
class EventOutputCenterScalar_: public EventOutput_<DIM> {
public:
	typedef EventOutput_<DIM> EventOutput;

	typedef Equation_<DIM> Equ;
	typedef Equ* pEqu;
	typedef const Equ* const_pEqu;

	typedef Scalar_<DIM> Scalar;
	typedef const Scalar_<DIM> const_Scalar;
	typedef const Scalar_<DIM>& const_ref_Scalar;
	typedef std::shared_ptr<Scalar> spScalar;
	typedef std::shared_ptr<const Scalar> spcScalar;
protected:
	std::string _vname;
	std::string _fname;

public:
	EventOutputCenterScalar_(const std::string& vn, const std::string& fn,
			int is = -1, int ie = -1, int istep = -1, int flag = 0) :
			EventOutput(is, ie, istep, flag, ""), _vname(vn), _fname(fn) {
	}

	int execute(St step, Vt t, int fob, pEqu pd = nullptr) {
		_record_value(step, t, pd);
		return -1;
	}

protected:
	void _record_value(const St& step, Vt& t, pEqu pd) {
		spScalar cs = pd->get_CS(this->_vname);
		//std::stringstream sst;
		//std::setprecision(3);
		//sst << this->_fname << "_" << step << "_" << t;
		std::string name = fmt::format("{:s}_{:d}_{:.4e}", this->_fname, step,
				t);
		Output_Scalar(name, *cs);
		//sst.str("");
	}
};

template<St DIM>
class EventOutputPointScalar_: public EventOutput_<DIM> {
public:
	typedef EventOutput_<DIM> EventOutput;

	typedef Equation_<DIM> Equ;
	typedef Equ* pEqu;
	typedef const Equ* const_pEqu;

	typedef Scalar_<DIM> Scalar;
	typedef const Scalar_<DIM> const_Scalar;
	typedef const Scalar_<DIM>& const_ref_Scalar;
	typedef std::shared_ptr<Scalar> spScalar;
	typedef std::shared_ptr<const Scalar> spcScalar;
protected:
	std::string _vname;
	std::string _fname_out;
	std::string _fname_in;

public:
	EventOutputPointScalar_(const std::string& vn, const std::string& fnin,
			const std::string& fnout, int is = -1, int ie = -1, int istep = -1,
			int flag = 0) :
			EventOutput(is, ie, istep, flag, ""), _vname(vn), _fname_in(fnin), _fname_out(
					fnout) {
	}

	int execute(St step, Vt t, int fob, pEqu pd = nullptr) {
		_record_value(step, t, pd);
		return -1;
	}

protected:
	void _record_value(const St& step, Vt& t, pEqu pd) {
		spScalar cs = pd->get_CS(this->_vname);
		//std::stringstream sst;
		//std::setprecision(3);
		//sst << this->_fname << "_" << step << "_" << t;
		std::string name = fmt::format("{:s}_{:d}_{:.4e}", this->_fname, step,
				t);
		//Vt res = Operation::InterpolateCoordinate(*cs,x, y,z);
		//Output_Scalar(name, *cs);
		//sst.str("");
	}
};

template<St DIM>
class EventSetCenterScalar_: public Event_<DIM> {
public:
	typedef Event_<DIM> Event;
	typedef Equation_<DIM> Equ;
	typedef Equ* pEqu;
	typedef const Equ* const_pEqu;

	typedef typename Equ::Function Function;

	typedef Scalar_<DIM> Scalar;
	typedef const Scalar_<DIM> const_Scalar;
	typedef const Scalar_<DIM>& const_ref_Scalar;
	typedef std::shared_ptr<Scalar> spScalar;
	typedef std::shared_ptr<const Scalar> spcScalar;
protected:
	std::string _vname;
	Function _fun;
public:
	EventSetCenterScalar_(const std::string& vn, Function fun, int is = -1,
			int ie = -1, int istep = -1, int flag = 0) :
			Event(is, ie, istep, flag), _vname(vn), _fun(fun) {
	}

	int execute(St step, Vt t, int fob, pEqu pd = nullptr) {
		_set_value(step, t, pd);
		return -1;
	}

protected:
	void _set_value(const St& step, Vt& t, pEqu pd) {
		pd->set_CS(_vname, _fun);
	}
};

template<St DIM>
class EventOutputTime_: public EventOutput_<DIM> {
public:
	typedef EventOutput_<DIM> EventOutput;

	typedef Equation_<DIM> Equ;
	typedef Equ* pEqu;
	typedef const Equ* const_pEqu;
	typedef carpio::tick_t tick_t;
protected:
	//St _count_exe;

	std::list<Vt> l_t;     //time
	std::list<St> l_step;  //step
	std::list<tick_t> l_cpu;   //cpu time
	std::list<tick_t> l_wall;  //wall time

public:
	EventOutputTime_(int is = -1, int ie = -1, int istep = -1, int flag = 0,
			std::string fn = "") :
			EventOutput(is, ie, istep, flag, fn) {
		//this->_count_exe = 0;
	}

	int execute(St step, Vt t, int fob, pEqu pd = nullptr) {
		//_count_exe++;
		_record_time(step, t);
		_output(step, t);
		return -1;
	}

	~EventOutputTime_() {
	}

protected:
	tick_t _d_last(const std::list<tick_t>& list) const {
		tick_t res = 0;
		if (list.size() < 2) {
			return res;
		}
		auto last = list.rbegin();
		tick_t cpu_l = (*last);
		last++;
		tick_t cpu_ll = (*last);
		return (cpu_l - cpu_ll);
	}

	void _output(const St& step, const Vt& t) {
		Vt d_cpu = carpio::Clock::TicksToSecondsD(_d_last(l_cpu));
		Vt d_wall = carpio::Clock::TicksToSecondsD(_d_last(l_wall));
		if (this->_flag_output == 0) {
			fmt::print(
					"step: {:>8d} t: {:>10.5e} cpu: {:>8.5f} wall: {:>8.5f}\n",
					step, t, d_cpu, d_wall);
		} else {
			fmt::print(this->_f,
					"step: {:>8d} t: {:>10.5e} cpu: {:>8.5f} wall: {:>8.5f}\n",
					step, t, d_cpu, d_wall);
		}
	}

	void _record_time(const St& step, Vt& t) {
		l_t.push_back(step);
		l_t.push_back(t);
		l_cpu.push_back(carpio::Clock::SystemTime());
		l_wall.push_back(carpio::Clock::Tick());
	}

	std::string name() const {
		return "OutputTime";
	}
};

template<St DIM>
class EventOutputError_: public EventOutput_<DIM> {
public:
	typedef EventOutput_<DIM> EventOutput;

	typedef Equation_<DIM> Equ;
	typedef Equ* pEqu;
	typedef const Equ* const_pEqu;

	typedef Scalar_<DIM> Scalar;
	typedef const Scalar_<DIM> const_Scalar;
	typedef const Scalar_<DIM>& const_ref_Scalar;
	typedef std::shared_ptr<Scalar> spScalar;
	typedef std::shared_ptr<const Scalar> spcScalar;

	typedef Operation_<DIM> Operation;
protected:
	std::string _name_exa;
	std::string _name_val;

	std::list<Vt> l_t;     //time
	std::list<St> l_step;  //step
	std::list<Vt> l_e1;   //
	std::list<Vt> l_e2;  //
	std::list<Vt> l_ei;  //

public:
	EventOutputError_(std::string name_exa, std::string name_val, int is = -1,
			int ie = -1, int istep = -1, int flag = 0, std::string fn = "") :
			EventOutput(is, ie, istep, flag, fn) {
		this->_name_exa = name_exa;
		this->_name_val = name_val;
	}

	int execute(St step, Vt t, int fob, pEqu pe = nullptr) {
		Vt e1 = _cal_error1(pe);
		Vt e2 = _cal_error2(pe);
		Vt ei = _cal_errori(pe);
		_output(step, t, e1, e2, ei);
		_record_error(step, t, e1, e2, ei);
		return -1;
	}

	~EventOutputError_() {

	}

protected:
	Vt _cal_error1(pEqu pe) {
		Scalar& exa = *(pe->get_CS(this->_name_exa));
		Scalar& val = *(pe->get_CS(this->_name_val));
		return Operation::Error1(exa, val);
	}
	Vt _cal_error2(pEqu pe) {
		Scalar& exa = *(pe->get_CS(this->_name_exa));
		Scalar& val = *(pe->get_CS(this->_name_val));
		return Operation::Error2(exa, val);
	}
	Vt _cal_errori(pEqu pe) {
		Scalar& exa = *(pe->get_CS(this->_name_exa));
		Scalar& val = *(pe->get_CS(this->_name_val));
		return Operation::Errori(exa, val);
	}

	void _output(const St& step, const Vt& t, Vt e1, Vt e2, Vt ei) {
		if (this->_flag_output == 0) {
			fmt::print(
					"step: {:>8d} t: {:>8.5e} e1: {:>8.5e} e2: {:>8.5e} ei: {:>8.5e}\n",
					step, t, e1, e2, ei);
		} else {
			fmt::print(this->_f,
					"step: {:>8d} t: {:>8.5e} e1: {:>8.5e} e2: {:>8.5e} ei: {:>8.5e}\n",
					step, t, e1, e2, ei);
		}
	}

	void _record_error(const St& step, Vt& t, Vt e1, Vt e2, Vt ei) {
		l_t.push_back(step);
		l_t.push_back(t);
		l_e1.push_back(e1);
		l_e2.push_back(e2);
		l_ei.push_back(ei);
	}

	std::string name() const {
		return "Error";
	}

};

template<St DIM>
class EventStopScalarError_: public Event_<DIM> {
public:
	static const St Dim = DIM;
	typedef Event_<DIM> Event;
	typedef Equation_<Dim> Equ;
	typedef Equ* pEqu;
	typedef const Equ* const_pEqu;

	typedef Scalar_<DIM> Scalar;
	typedef const Scalar_<DIM> const_Scalar;
	typedef const Scalar_<DIM>& const_ref_Scalar;
	typedef std::shared_ptr<Scalar> spScalar;
	typedef std::shared_ptr<const Scalar> spcScalar;

	typedef Operation_<DIM> Operation;
protected:
	std::string _vname;
	spScalar _prev;
	Vt _e[3];
public:
	EventStopScalarError_(int is, int ie, int istep, int flag,
			std::string vname, Vt e1 = -1, Vt e2 = -1, Vt e3 = -1) :
			Event(is, ie, istep, flag) {
		_prev = nullptr;
		_e[0] = e1;
		_e[1] = e2;
		_e[2] = e3;
		this->_vname = vname;

	}

	int execute(St step, Vt t, int fob, pEqu pd = nullptr) {
		bool f[3] = { false, false, false };
		_check_error(step, t, pd, f);
		if (f[0] || f[1] || f[2]) {
			_set_stop_flag(pd);
			//std::cout << "stop --------------------------------------\n";
		}
		_save_scalar(step, t, pd);

		return -1;
	}

	std::string name() const {
		return "EventStopScalarError";
	}

	virtual ~EventStopScalarError_() {
	}
protected:
	void _save_scalar(const St& step, Vt& t, pEqu pd) {
		spScalar cs = pd->get_CS(this->_vname);
		if (this->_prev == nullptr) {
			this->_prev = spScalar(new Scalar(*cs));
		} else {
			Operation::Copy(*cs, *(this->_prev));
		}
	}

	void _set_stop_flag(pEqu pd) {
		pd->add_event_flag("_STOP_");
	}

	void _check_error(const St& step, Vt& t, pEqu pd, bool* f) {
		// if error less than _e is true;
		if (_prev == nullptr) {
			return;
		}
		spScalar cs = pd->get_CS(this->_vname);
		f[0] = _e[0] > 0 ? Operation::Error1(*_prev, *cs) < _e[0] : false;
		f[1] = _e[1] > 0 ? Operation::Error2(*_prev, *cs) < _e[1] : false;
		f[2] = _e[2] > 0 ? Operation::Errori(*_prev, *cs) < _e[2] : false;
		fmt::print("step: {:5d} e1 {:7.3e} e2 {:7.3e} ei:{:7.3e}\n", step,
				Operation::Error1(*_prev, *cs), Operation::Error2(*_prev, *cs),
				Operation::Errori(*_prev, *cs));
	}

};

template<St DIM>
class EventFlag_: public Event_<DIM> {
public:
	static const St Dim = DIM;
	typedef Event_<DIM> Event;
	typedef Equation_<Dim> Equ;
	typedef Equ* pEqu;
	typedef const Equ* const_pEqu;
protected:
	int _flag;
	std::string _name; // the name of flag, it is also the key

public:
	EventFlag_(const std::string& name, int f) :
			Event(-1, -1, -1) {
		this->_flag = f;
		this->_name = name;
	}

	~EventFlag_() {

	}

	bool _do_execute(St step, Vt t, int fob) const {
		return false;
	}

	int execute(St step, Vt t, int fob, pEqu pd = nullptr) {
		return -1;
	}

	int flag() const {
		return _flag;
	}

	void set_flag(int i) {
		_flag = i;
	}

	std::string name() const {
		return _name;
	}
};

template<St DIM, class DATA>
class EventFlagData_: public Event_<DIM> {
public:
	static const St Dim = DIM;
	typedef Event_<DIM> Event;

	typedef DATA Data;
	typedef Equation_<Dim> Equ;
	typedef Equ* pEqu;
	typedef const Equ* const_pEqu;
protected:
	int _flag;
	Data _data;
	std::string _name; // the name of flag, it is also the key

public:
	EventFlagData_(const std::string& name, int f, Data val) :
			Event(-1, -1, -1) {
		this->_flag = f;
		this->_name = name;
		this->_data = val;
	}

	~EventFlagData_() {

	}

	bool _do_execute(St step, Vt t, int fob) const {
		return false;
	}

	int execute(St step, Vt t, int fob, pEqu pd = nullptr) {
		return -1;
	}

	int flag() const {
		return _flag;
	}

	void set_flag(int i) {
		_flag = i;
	}

	const Data& data() const {
		return this->_data;
	}

	void set_data(Data v) {
		this->_data = v;
	}

	std::string name() const {
		return _name;
	}
};

}

#endif
