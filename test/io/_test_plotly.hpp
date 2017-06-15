#ifndef __TEST_PLOTLY_HPP_
#define __TEST_PLOTLY_HPP_

//#include "../test_define.hpp"
#include <io/plotly.h>
#include <io/plotly_actor.h>
#include "Python.h"
#include <iostream>
using namespace std;
namespace carpio {

TEST(PLOTLY, test){
	Plotly p;
	std::cout<< "version : " << p.version()<<std::endl;
	typedef ArrayListV<double> Arr;
	Arr arrx(3);
	Arr arry(3);
	arrx.assign_forward(0,1);
	arry.assign_forward(0,2);
	arry[2] = 5;
	PlotlyActor::spPA actor = PlotlyActor::XY(arrx, arry);
	actor->set_name("line one");
	actor->set_opacity(0.5);
	p.add(actor);
	//p.plot();
}
}
#endif
