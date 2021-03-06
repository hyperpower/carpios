/* Copyright 2010 Jukka Jylanki

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License. */
#pragma once

/**
 *  @file  Clock.h
 *  @brief The Clock class. Supplies timing facilities.
 **/

#ifdef __EMSCRIPTEN__
// The native type for high-resolution timing is double, use
// that instead of uint64, which is not natively supported but
// must be emulated, which is slow.
#define _TICK_IS_FLOAT
#include <limits>
#endif


#include "format.h"
#include <iostream>
#include <iterator>
#include <list>
#include "../type_define.hpp"

namespace carpio {

/// A tick is the basic unit of the high-resolution timer. If MATH_TICK_IS_FLOAT is defined,
/// then tick_t is a floating-point type. Otherwise a 64-bit unsigned integer is used instead.
#ifdef _TICK_IS_FLOAT
typedef double tick_t;
const tick_t TICK_INF = std::numeric_limits<double>::infinity();
#else
typedef unsigned long long tick_t;
const tick_t TICK_INF = 0xFFFFFFFFFFFFFFFFULL;
#endif

/**
 *  @brief High-resolution timing and system time.
 *
 *  Gives out timing information in various forms. Use this rather than
 *  any platform-dependent perf-counters or rdtsc or whatever.
 */
class Clock {
public:
	Clock();
	Clock(const std::string& name);

//	~Clock() {}

/// Sleeps the current thread for the given amount of milliseconds.
	static void Sleep(int milliseconds);

	/// @return The current year.
	static int Year();

	/// @return The current month, [1,12]
	static int Month();

	/// @return The current day of month, [1,31]
	static int Day();

	/// @return The current hour of day, [0,23]
	static int Hour();

	/// @return The current clock minute, [0,59]
	static int Min();

	/// @return The current clock second, [0,59]
	static int Sec();

	/// @return The current system time counter in milliseconds.
	static unsigned long SystemTime();

	/// @return The number of ticks since application start.
	static tick_t TimeElapseInTick();

	/// @return The number of milliseconds since application start.
	static inline float TimeElapseInMillisecondsF() {
		return TimespanToMillisecondsF(appStartTime, Tick());
	}
	/// @return The number of milliseconds since application start.
	static inline double TimeElapseInMillisecondsD() {
		return TimespanToMillisecondsD(appStartTime, Tick());
	}

	/// @return The number of seconds since application start.
	static inline float TimeElapseInSecondsF() {
		return TimespanToSecondsF(appStartTime, Tick());
	}
	/// @return The number of seconds since application start.
	static inline double TimeElapseInSecondsD() {
		return TimespanToSecondsD(appStartTime, Tick());
	}

	/// @return The low part of the current tick-time (using whatever high-resolution counter available)
	static unsigned long TickU32();

	/// @return The current tick-time (using whatever high-resolution counter available)
	static tick_t Tick();

	/// @return How many ticks make up a second.
	static tick_t TicksPerSec();

	static inline tick_t TicksPerMillisecond() {
		return TicksPerSec() / 1000;
	}

	/// @returns the number of ticks occurring between the two wallclock times.
	static inline tick_t TicksInBetween(tick_t newTick, tick_t oldTick) {
		return (tick_t) (newTick - oldTick);
	}

	/// @returns true if newTick represents a later wallclock time than oldTick.
	static inline bool IsNewer(tick_t newTick, tick_t oldTick) {
#ifdef _TICK_IS_FLOAT
		return newTick >= oldTick;
#else
		return TicksInBetween(newTick, oldTick) < ((tick_t) (-1) >> 1);
#endif
	}

	static inline float MillisecondsSinceF(tick_t oldTick) {
		return TimespanToMillisecondsF(oldTick, Tick());
	}
	static inline double MillisecondsSinceD(tick_t oldTick) {
		return TimespanToMillisecondsD(oldTick, Tick());
	}

	static inline float SecondsSinceF(tick_t oldTick) {
		return TimespanToSecondsF(oldTick, Tick());
	}
	static inline double SecondsSinceD(tick_t oldTick) {
		return TimespanToSecondsD(oldTick, Tick());
	}

	static inline float TicksToMillisecondsF(tick_t ticks) {
		return ticks * 1000.0 / (float) TicksPerSec();
	}
	static inline double TicksToMillisecondsD(tick_t ticks) {
		return ticks * 1000.0 / (double) TicksPerSec();
	}

	static inline float TicksToSecondsF(tick_t ticks) {
		return ticks / (float) TicksPerSec();
	}
	static inline double TicksToSecondsD(tick_t ticks) {
		return ticks / (double) TicksPerSec();
	}

	static inline float TimespanToMillisecondsF(tick_t oldTick,
			tick_t newTick) {
		return TicksToMillisecondsF(TicksInBetween(newTick, oldTick));
	}
	static inline double TimespanToMillisecondsD(tick_t oldTick,
			tick_t newTick) {
		return TicksToMillisecondsD(TicksInBetween(newTick, oldTick));
	}

	static inline float TimespanToSecondsF(tick_t oldTick, tick_t newTick) {
		return TicksToSecondsF(TicksInBetween(newTick, oldTick));
	}
	static inline double TimespanToSecondsD(tick_t oldTick, tick_t newTick) {
		return TicksToSecondsD(TicksInBetween(newTick, oldTick));
	}

	/*
	 * start a clock
	 */
	void start() {
		this->_start_time_cpu = Clock::SystemTime();
		this->_start_time_wt = Clock::Tick();
	}

	void break_point(const std::string& name = "", const Float& num = 1) {
		_tp_name.push_back(name);
		_tp_cpu.push_back(Clock::SystemTime());
		_tp_wt.push_back(Clock::Tick());
		_tp_num.push_back(num);
	}

	void clear_records() {
		_tp_cpu.clear();
		_tp_wt.clear();
		_tp_name.clear();
		_tp_num.clear();
	}

	St size() const {
		return _tp_name.size();
	}

protected:
	typedef std::list<tick_t>::iterator iterator_t;
	typedef std::list<tick_t>::const_iterator const_iterator_t;

	typedef std::list<Float>::iterator iterator_n;
	typedef std::list<Float>::const_iterator const_iterator_n;

public:
	tick_t start_time_cpu() const {
		return this->_start_time_cpu;
	}

	tick_t start_time_wall() const {
		return this->_start_time_wt;
	}

	iterator_t begin_cpu() {
		return this->_tp_cpu.begin();
	}

	const_iterator_t begin_cpu() const {
		return this->_tp_cpu.begin();
	}

	iterator_t end_cpu() {
		return this->_tp_cpu.end();
	}

	const_iterator_t end_cpu() const {
		return this->_tp_cpu.end();
	}

	iterator_t begin_wall() {
		return this->_tp_wt.begin();
	}

	const_iterator_t begin_wall() const {
		return this->_tp_wt.begin();
	}

	iterator_t end_wall() {
		return this->_tp_wt.end();
	}

	const_iterator_t end_wall() const {
		return this->_tp_wt.end();
	}

	iterator_n begin_num() {
		return this->_tp_num.begin();
	}

	const_iterator_n begin_num() const {
		return this->_tp_num.begin();
	}

	iterator_n end_num() {
		return this->_tp_num.end();
	}

	const_iterator_n end_num() const {
		return this->_tp_num.end();
	}

	std::string name() const{
		return this->_name;
	}

	void show() {
		// get leatest time
		tick_t d_end_cpu, d_end_wt;
		fmt::print("{:<7} {:^10} {:^10} {:^5} {:^5} {:^7} {}\n", " Index",
				"Cpu", "Wall", "Cpu%", "Wall%", "Num", "Name");
		fmt::print("{:-<60}\n", " ");
		if (_tp_cpu.empty()) {
			std::cout << " Empty\n";
			return;
		} else {
			d_end_cpu = Clock::TicksInBetween(*(_tp_cpu.rbegin()),
					_start_time_cpu);
			d_end_wt = Clock::TicksInBetween(*(_tp_wt.rbegin()),
					_start_time_wt);
		}
		std::list<tick_t>::iterator iter_c = _tp_cpu.begin(), iter_wt =
				_tp_wt.begin();
		std::list<std::string>::iterator iter_name = _tp_name.begin();
		iterator_n iter_num = _tp_num.begin();

		// loop --------------------
		int i = 0;
		while (iter_c != _tp_cpu.end()) {
			tick_t dt_cpu, dt_wt;
			if (iter_c == _tp_cpu.begin()) {
				dt_cpu = Clock::TicksInBetween(*iter_c, _start_time_cpu);
				dt_wt = Clock::TicksInBetween(*iter_wt, _start_time_wt);
			} else {
				// c++11
				auto prev_c = std::prev(iter_c, 1);
				auto prev_wt = std::prev(iter_wt, 1);
				dt_cpu = Clock::TicksInBetween(*iter_c, *prev_c);
				dt_wt = Clock::TicksInBetween(*iter_c, *prev_wt);
			}
			// output
			fmt::print(
					"{: ^7} {: >10.5f} {: >10.5f} {: >5.2f} {: >5.2f} {: >4.1e} {}\n",
					i,  // Index
					Clock::TicksToSecondsD(dt_cpu), // d Cpu time
					Clock::TicksToSecondsD(dt_wt),  // d Wall time
					Clock::TicksToSecondsD(dt_cpu)
							/ Clock::TicksToSecondsD(d_end_cpu) * 100,
					Clock::TicksToSecondsD(dt_wt)
							/ Clock::TicksToSecondsD(d_end_wt) * 100,
					(*iter_num), (*iter_name));

			//std::cout << "  " << i << "  " << Clock::TicksToSecondsD(dt_cpu)
			//		<< "   " << Clock::TicksToSecondsD(dt_wt) << "   "
			//		<< Clock::TicksToSecondsD(dt_cpu)
			//				/ Clock::TicksToSecondsD(d_end_cpu) * 100 << "   "
			//		<< Clock::TicksToSecondsD(dt_wt)
			//				/ Clock::TicksToSecondsD(d_end_wt) * 100 << "   "
			//		<< (*iter_name) << "   " << "\n";
			//
			++iter_c;
			++iter_wt;
			++iter_name;
			++iter_num;
			++i;
		}
		fmt::print("{:-<60}<\n", " ");
		fmt::print("{: ^7} {: >10.5f} {: >10.5f}\n", "sum",
				Clock::TicksToSecondsD(d_end_cpu),
				Clock::TicksToSecondsD(d_end_wt));

	}

protected:
	static tick_t appStartTime;      ///< Application startup time in ticks.

	/// Initializes clock tick frequency and marks the application startup time.
	static void InitClockData();

	tick_t _start_time_cpu;
	tick_t _start_time_wt;

	std::string _name;          // clock name
	std::list<tick_t> _tp_cpu;  // time points of cpu time
	std::list<tick_t> _tp_wt;   // time points of wall time
	std::list<std::string> _tp_name;
	std::list<Float> _tp_num;   // store

#ifdef WIN32
	// The following two are actually used as LARGE_INTEGERs, but to avoid having to pull Windows.h in Clock.h, define them
	// as identically sized u64 instead.
	static u64/*LARGE_INTEGER*/ddwTimerFrequency;///< Ticks per second.
#endif
#ifdef __APPLE__
	static tick_t ticksPerSecond;
#endif
};

}
